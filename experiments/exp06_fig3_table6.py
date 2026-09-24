"""Section 4.4 - Figure 3 and Table 6: five-point pathways {0, 0.05, 0.1, t3, 1}.

dt0 = 0.05, K = 0.001, t3 in {0.20, 0.25, 0.30}; case109272 ("109k"),
case36964 ("36k") and case_ACTIVSg2000limit ("lim"); FE, RK2 or BE followed by NR.
Abscissa: index k of t_k for 0..4, then NR iteration i at 4 + i.
Table 6 reports the median CPU time (homotopy + NR) relative to NR-MAT.
"""

from __future__ import annotations

import argparse

import numpy as np

from common import benchmark, label, savefig, write_table
from dynhomotopy.datasets import SYSTEMS, load_problem
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.solvers import newton_raphson
from paper_values import TABLE6
from plots import CASE_STYLE, METHOD_COLOR, METHOD_MARKER, new_figure

K_44 = 1e-3
CASES = ["case109272", "case36964", "case_ACTIVSg2000limit"]
T3 = [0.20, 0.25, 0.30]
METHODS = ["FE", "RK2", "BE"]


def run(pf, method, t3):
    return solve_hybrid(pf, pf.flat_start(), K_44, pathway.explicit(0, 0.05, 0.1, t3, 1.0), method, "NR")


def main(reps: int = 5):
    results, times = {}, {}
    for name in CASES:
        pf = load_problem(name)
        fns = {("NR-MAT", None, name): lambda pf=pf: newton_raphson(pf, pf.case_start())}
        for method in METHODS:
            for t3 in T3:
                res = run(pf, method, t3)
                results[(method, t3, name)] = res
                times[(method, t3, name)] = None
                if res.converged:
                    fns[(method, t3, name)] = lambda pf=pf, m=method, t=t3: run(pf, m, t)
        times.update(benchmark(fns, reps))
        print("done", name)

    for t3, tag in zip(T3, "abc"):
        fig, ax = new_figure(width=7, height=4.6)
        for name in CASES:
            short = SYSTEMS[name].short
            for method in METHODS:
                r = results[(method, t3, name)]
                norms = r.trajectory.norms + (r.refine.norms[1:] if r.refine else [])
                ax.plot(np.arange(len(norms)), np.log10(norms), CASE_STYLE[short], color=METHOD_COLOR[method],
                        marker=METHOD_MARKER[method], markersize=5, markeredgecolor="#fcfcfb",
                        label=f"{short}-{method}")
        ax.set(xlabel="Iterations", ylabel=r"$\log_{10}\|g(x)\|_\infty$", title=f"Fig. 3({tag}) pathway with t3 = {t3}",
               ylim=(-12.5, 4.5))
        ax.axvline(4, color="#52514e", lw=0.8, ls=(0, (2, 3)))
        ax.legend(ncol=3, fontsize=8, loc="lower left")
        savefig(fig, f"fig3{tag}_pathway_t3_{t3:.2f}")

    # ---- Table 6 -------------------------------------------------------------------
    rows = []
    for method, t3 in [("NR-MAT", None)] + [(m, t) for m in METHODS for t in T3]:
        row = [method, "-" if t3 is None else f"{t3:.2f}"]
        paper = TABLE6[(method, t3)]
        for j, name in enumerate(CASES):
            t, ref = times[(method, t3, name)], times[("NR-MAT", None, name)]
            p = paper[j]
            p_txt = "fail" if isinstance(p, str) else f"{100 * p / TABLE6[('NR-MAT', None)][j]:.0f}%"
            ours = "fail" if t is None else f"{t:.3f} s ({100 * t / ref:.0f}%)"
            iters = "" if method == "NR-MAT" else f" [NR {results[(method, t3, name)].refine.iterations}]" if t else ""
            row.append(f"{ours}{iters} / {p_txt}")
        rows.append(row)
    write_table("table6", ["method", "t3"] + [label(c) for c in CASES], rows,
                f"Table 6: median CPU time over {reps} interleaved runs, % of NR-MAT (ours / paper %)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=5)
    main(ap.parse_args().reps)
