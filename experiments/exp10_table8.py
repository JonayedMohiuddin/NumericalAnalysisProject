"""Table 8: run time of each method on the ill-conditioned cases, relative to NR-MAT.

Only the solve is timed; Ybus is built beforehand. The paper averages 100
MATLAB runs, we take the median of `reps` runs.
"""

import argparse

from dynhomotopy.datasets import load_problem

from .common import benchmark, label, write_table
from .paper_values import TABLE8, TABLE8_CASES
from .setups import methods

ROWS = ["NR-MAT", "GSH-NR", "BE(NR)", "RK2(NR)", "BE(FDXB)", "RK2(FDXB)"]


def paper_percent(method, j):
    value, ref = TABLE8[method][j], TABLE8["NR-MAT"][j]
    return value if isinstance(value, str) else f"{100 * value / ref:.0f}%"


def main(cases=TABLE8_CASES, reps=5):
    times = {}
    for name in cases:
        pf = load_problem(name)
        to_time = {}
        for m, fn in methods(name).items():
            if m not in ROWS:
                continue
            converged, _ = fn(pf)
            times[(m, name)] = None
            if converged:
                to_time[(m, name)] = lambda fn=fn, pf=pf: fn(pf)
        times.update(benchmark(to_time, reps))

    rows = []
    for m in ROWS:
        row = [m]
        for name in cases:
            t, ref = times[(m, name)], times[("NR-MAT", name)]
            ours = "fail" if t is None else f"{t:.3f} s ({100 * t / ref:.0f}%)"
            row.append(f"{ours} / {paper_percent(m, TABLE8_CASES.index(name))}")
        rows.append(row)
    write_table("table8", ["method"] + [label(c) for c in cases], rows,
                f"Table 8: median CPU time of {reps} runs, % of NR-MAT (ours / paper)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=5)
    main(reps=parser.parse_args().reps)
