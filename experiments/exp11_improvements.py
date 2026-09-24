"""Our changes to the method, each switched on and off.

The three changes that help (multiplier, corrector, adaptive) are run in all
8 on/off combinations. The four that did not help (BE-chord, RK4, the scaled
and the Newton homotopy) are run on their own. Each configuration runs on 34 cases
(the paper's 9 plus 25 more from the same Zenodo records) under 5 settings
of (dt0, K), with the path {0, dt0, 1}. A run is solved if NR then converges
within 10 iterations, as in the paper. Cost is the number of LU
factorizations, which does not depend on the machine.

When the paper's method fails on a case that an improvement solves, we rerun
the paper's method with 30 NR iterations to see whether it was only slow or
really diverged.
"""

import csv
import itertools

import numpy as np

from improvements import Options, solve
from improvements.cases import ALL, load

from .common import TABLES, write_table

SETTINGS = {
    "S1 paper": (0.005, 1e-4),
    "S2 Sec 4.4": (0.05, 1e-3),
    "S3 Sec 4.4": (0.1, 2e-3),
    "S4 weak K": (0.005, 1e-5),
    "S5 strong K": (0.005, 1e-3),
}
COMBINATIONS = [Options(*flags) for flags in itertools.product([False, True], repeat=3)]
COMBINATIONS += [Options(step="BE-chord"), Options(step="RK4"),
                 Options(homotopy="scaled"), Options(homotopy="newton")]
PAPER = Options()


def run(name, setting, options, max_it=10):
    dt0, K = SETTINGS[setting]
    pf = load(name)
    return solve(pf, pf.flat_start(), K, [0.0, dt0, 1.0], options, max_it=max_it)


def run_everything():
    results = {}
    for name in ALL:
        for setting in SETTINGS:
            for options in COMBINATIONS:
                r = run(name, setting, options)
                results[(name, setting, options.name())] = (r.converged, r.factorizations, r.time)
        print(f"{name} done", flush=True)

    with open(TABLES / "improvements_raw.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["case", "setting", "options", "converged", "lu", "time"])
        for (name, setting, options), (ok, lu, t) in results.items():
            writer.writerow([name, setting, options, ok, lu, f"{t:.4f}"])
    return results


def table_solved(results):
    rows = []
    for options in COMBINATIONS:
        counts = [sum(results[(c, s, options.name())][0] for c in ALL) for s in SETTINGS]
        rows.append([options.name()] + counts + [sum(counts)])
    write_table("improvements_solved", ["options"] + list(SETTINGS) + ["total"], rows,
                f"Cases solved out of {len(ALL)} (NR converges within 10 iterations)")


def table_cost(results):
    rows = []
    for options in COMBINATIONS:
        row = [options.name()]
        for s in SETTINGS:
            common = [c for c in ALL if results[(c, s, PAPER.name())][0] and results[(c, s, options.name())][0]]
            ours = np.mean([results[(c, s, options.name())][1] for c in common])
            paper = np.mean([results[(c, s, PAPER.name())][1] for c in common])
            row.append(f"{ours:.2f} vs {paper:.2f}" if common else "-")
        times = [results[(c, "S1 paper", options.name())][2] / results[(c, "S1 paper", PAPER.name())][2]
                 for c in ALL
                 if results[(c, "S1 paper", PAPER.name())][0] and results[(c, "S1 paper", options.name())][0]]
        row.append(f"{np.median(times):.2f}")
        rows.append(row)
    write_table("improvements_cost", ["options"] + list(SETTINGS) + ["S1 median time ratio"], rows,
                "Mean LU factorizations on the cases solved by both this combination and the paper's "
                "method (ours vs paper), and median run time relative to the paper's method")


def table_gains(results):
    """Cases the paper's method fails on but some combination solves."""
    rows = []
    for s in SETTINGS:
        for c in ALL:
            if results[(c, s, PAPER.name())][0]:
                continue
            solved_by = [o.name() for o in COMBINATIONS if results[(c, s, o.name())][0]]
            if not solved_by:
                continue
            slow = run(c, s, PAPER, max_it=30)
            kind = f"slow (converges in {slow.refine.iterations})" if slow.converged else "diverges"
            rows.append([s, c, kind, ", ".join(solved_by)])
    write_table("improvements_gains", ["setting", "case", "paper's method with 30 NR iterations",
                                       "solved within 10 by"], rows,
                "Cases the paper's method fails on that an improvement solves")


def table_cases(results):
    rows = []
    for c in ALL:
        cells = []
        for options in COMBINATIONS:
            ok, lu, _ = results[(c, "S1 paper", options.name())]
            cells.append(str(lu) if ok else "fail")
        rows.append([c] + cells)
    write_table("improvements_cases_S1", ["case"] + [o.name() for o in COMBINATIONS], rows,
                "Setting S1: LU factorizations per case, or fail")


def main():
    results = run_everything()
    table_solved(results)
    table_cost(results)
    table_gains(results)
    table_cases(results)


if __name__ == "__main__":
    main()
