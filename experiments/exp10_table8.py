"""Section 4.8 - Table 8: CPU time of each method on the ill-conditioned cases.

Times cover the solution stage (homotopy + refiner, or the NR/GSH iterations)
on a pre-built Ybus, as the median of ``--reps`` interleaved runs (the paper
averages 100 MATLAB runs). Percentages are relative to NR-MAT; the paper's
value follows "/".
"""

from __future__ import annotations

import argparse

from common import benchmark, label, write_table
from dynhomotopy.datasets import load_problem
from paper_values import TABLE8, TABLE8_CASES
from setups import methods

ROWS = ["NR-MAT", "GSH-NR", "BE(NR)", "RK2(NR)", "BE(FDXB)", "RK2(FDXB)"]


def paper_pct(method: str, j: int) -> str:
    v, ref = TABLE8[method][j], TABLE8["NR-MAT"][j]
    return "fail" if isinstance(v, str) else f"{100 * v / ref:.0f}%"


def main(cases=TABLE8_CASES, reps: int = 5):
    times = {}
    for name in cases:
        pf = load_problem(name)
        fns, timed = methods(name), {}
        for m in ROWS:
            ok, _ = fns[m](pf)  # warm-up; also decides success
            times[(m, name)] = None
            if ok:
                timed[(m, name)] = lambda fn=fns[m], pf=pf: fn(pf)
        times.update(benchmark(timed, reps))
        print("done", name)
    rows = []
    for m in ROWS:
        row = [m]
        for j, name in enumerate(cases):
            t, ref = times[(m, name)], times[("NR-MAT", name)]
            ours = "fail" if t is None else f"{t:.3f} s ({100 * t / ref:.0f}%)"
            row.append(f"{ours} / {paper_pct(m, TABLE8_CASES.index(name))}")
        rows.append(row)
    write_table("table8", ["method"] + [label(c) for c in cases], rows,
                f"Table 8: median CPU time over {reps} interleaved runs, % of NR-MAT (ours / paper %)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=5)
    main(reps=ap.parse_args().reps)
