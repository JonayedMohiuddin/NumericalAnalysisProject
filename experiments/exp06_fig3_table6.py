"""Figure 3 and Table 6: five point paths {0, 0.05, 0.1, t3, 1} with K = 0.001.

Index 0-4 on the x axis are the homotopy points and 5, 6, ... are NR iterations.
Table 6 gives the time of homotopy + NR relative to plain NR from the case file (NR-MAT).
"""

import argparse

import numpy as np

from dynhomotopy.datasets import SYSTEMS, load_problem
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.solvers import newton_raphson

from .common import benchmark, label, savefig, write_table
from .paper_values import TABLE6
from .plots import BACKGROUND, CASE_LINESTYLE, METHOD_COLOR, METHOD_MARKER, new_figure

K = 1e-3
CASES = ["case109272", "case36964", "case_ACTIVSg2000limit"]
T3_VALUES = [0.20, 0.25, 0.30]
METHODS = ["FE", "RK2", "BE"]


def run(pf, method, t3):
    times = pathway.explicit(0, 0.05, 0.1, t3, 1)
    return solve_hybrid(pf, pf.flat_start(), K, times, method, "NR")


def plot(results, t3, panel):
    fig, ax = new_figure(width=7, height=4.6)
    for name in CASES:
        short = SYSTEMS[name].short
        for method in METHODS:
            r = results[(method, t3, name)]
            norms = r.trajectory.norms + (r.refine.norms[1:] if r.refine else [])
            ax.plot(np.log10(norms), CASE_LINESTYLE[short], color=METHOD_COLOR[method],
                    marker=METHOD_MARKER[method], markersize=5, markeredgecolor=BACKGROUND,
                    label=f"{short}-{method}")
    ax.axvline(4, color="#52514e", lw=0.8, ls=":")
    ax.set(xlabel="Iterations", ylabel=r"$\log_{10}\|g(x)\|_\infty$", ylim=(-12.5, 4.5),
           title=f"Fig. 3({panel}) path with t3 = {t3}")
    ax.legend(ncol=3, fontsize=8, loc="lower left")
    savefig(fig, f"fig3{panel}_pathway_t3_{t3:.2f}")


def main(reps=5):
    results, times = {}, {}
    for name in CASES:
        pf = load_problem(name)
        to_time = {("NR-MAT", None, name): lambda pf=pf: newton_raphson(pf, pf.case_start())}
        for method in METHODS:
            for t3 in T3_VALUES:
                res = run(pf, method, t3)
                results[(method, t3, name)] = res
                times[(method, t3, name)] = None
                if res.converged:
                    to_time[(method, t3, name)] = lambda pf=pf, m=method, t3=t3: run(pf, m, t3)
        times.update(benchmark(to_time, reps))

    for t3, panel in zip(T3_VALUES, "abc"):
        plot(results, t3, panel)

    rows = []
    paper_ref = TABLE6[("NR-MAT", None)]
    for method, t3 in [("NR-MAT", None)] + [(m, t3) for m in METHODS for t3 in T3_VALUES]:
        row = [method, "-" if t3 is None else f"{t3:.2f}"]
        for j, name in enumerate(CASES):
            t, ref = times[(method, t3, name)], times[("NR-MAT", None, name)]
            paper = TABLE6[(method, t3)][j]
            paper = paper if isinstance(paper, str) else f"{100 * paper / paper_ref[j]:.0f}%"
            if t is None:
                ours = "fail"
            else:
                ours = f"{t:.3f} s ({100 * t / ref:.0f}%)"
                if method != "NR-MAT":
                    ours += f" [NR {results[(method, t3, name)].refine.iterations}]"
            row.append(f"{ours} / {paper}")
        rows.append(row)
    write_table("table6", ["method", "t3"] + [label(c) for c in CASES], rows,
                f"Table 6: median CPU time of {reps} runs, % of NR-MAT (ours / paper)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=5)
    main(parser.parse_args().reps)
