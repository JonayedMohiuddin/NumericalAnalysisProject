"""Our changes to the method, each switched on and off.

The multiplier, corrector and adaptive switches are run in all 8 on/off
combinations. Richardson step control is run with BE, FE, RK2 and RK4 and
together with the multiplier. BE-chord, fixed-step RK4 and the scaled and
Newton homotopies are run on their own.

Each configuration runs on 34 cases (the paper's 9 plus 25 more from the same
Zenodo records) under 5 settings of (dt0, K), with the path {0, dt0, 1}.
A run counts as solved only if NR converges within 10 iterations (as in the
paper) AND reaches the reference operating point of the case
(improvements.cases.reference_voltage). A run that converges to another
root is counted separately.

Cost is the number of LU factorizations. Run times are not reported here,
because every configuration is run only once.

To run only some configurations (the paper's method is always included):
    python -m experiments.exp11_improvements --configs PC "OM + PC" richardson
The tables of such a run are written as improvements_selected_*.
"""

import argparse
import csv
import itertools

import numpy as np

from improvements import Options, solve
from improvements.cases import ALL, load, on_reference

from .common import TABLES, write_table

SETTINGS = {
    "S1 paper": (0.005, 1e-4),
    "S2 Sec 4.4": (0.05, 1e-3),
    "S3 Sec 4.4": (0.1, 2e-3),
    "S4 weak K": (0.005, 1e-5),
    "S5 strong K": (0.005, 1e-3),
}
COMBINATIONS = [Options(*flags) for flags in itertools.product([False, True], repeat=3)]
COMBINATIONS += [Options(richardson=True, step=s) for s in ["BE", "FE", "RK2", "RK4"]]
COMBINATIONS += [Options(richardson=True, multiplier=True)]
COMBINATIONS += [Options(step="BE-chord"), Options(step="RK4"),
                 Options(homotopy="scaled"), Options(homotopy="newton")]
PAPER = Options()

SOLVED, OTHER_ROOT, FAILED = "solved", "other root", "failed"


def run(name, setting, options, max_it=10):
    dt0, K = SETTINGS[setting]
    pf = load(name)
    res = solve(pf, pf.flat_start(), K, [0.0, dt0, 1.0], options, max_it=max_it)
    if not res.converged:
        return FAILED, res
    return (SOLVED if on_reference(pf, res.refine.x) else OTHER_ROOT), res


def run_everything(configs, prefix):
    results = {}
    for name in ALL:
        for setting in SETTINGS:
            for options in configs:
                outcome, res = run(name, setting, options)
                results[(name, setting, options.name())] = (outcome, res.factorizations)
        print(f"{name} done", flush=True)

    with open(TABLES / f"{prefix}_raw.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["case", "setting", "options", "outcome", "lu"])
        for (name, setting, options), (outcome, lu) in results.items():
            writer.writerow([name, setting, options, outcome, lu])
    return results


def count(results, options, setting, outcome):
    return sum(results[(c, setting, options.name())][0] == outcome for c in ALL)


def table_solved(results, configs, prefix):
    rows = []
    for options in configs:
        solved = [count(results, options, s, SOLVED) for s in SETTINGS]
        other = sum(count(results, options, s, OTHER_ROOT) for s in SETTINGS)
        lost = sum(results[(c, s, PAPER.name())][0] == SOLVED and results[(c, s, options.name())][0] != SOLVED
                   for c in ALL for s in SETTINGS)
        rows.append([options.name()] + solved + [sum(solved), other, lost])
    write_table(f"{prefix}_solved",
                ["options"] + list(SETTINGS) + ["total solved", "other root", "lost vs paper"], rows,
                f"Runs that reach the reference operating point, out of {len(ALL)} per setting. "
                "'other root' counts runs that converged to a different solution")


def table_cost(results, configs, prefix):
    rows = []
    for options in configs:
        row = [options.name()]
        for s in SETTINGS:
            common = [c for c in ALL
                      if results[(c, s, PAPER.name())][0] == SOLVED and results[(c, s, options.name())][0] == SOLVED]
            if not common:
                row.append("-")
                continue
            ours = np.mean([results[(c, s, options.name())][1] for c in common])
            paper = np.mean([results[(c, s, PAPER.name())][1] for c in common])
            row.append(f"{ours:.2f} vs {paper:.2f} ({len(common)})")
        rows.append(row)
    write_table(f"{prefix}_cost", ["options"] + list(SETTINGS), rows,
                "Mean LU factorizations on the cases solved by both this configuration and the "
                "paper's method (ours vs paper, number of cases)")


def table_gains(results, configs, prefix):
    """Cases the paper's method does not solve but some configuration does."""
    rows = []
    for s in SETTINGS:
        for c in ALL:
            if results[(c, s, PAPER.name())][0] == SOLVED:
                continue
            solved_by = [o.name() for o in configs if results[(c, s, o.name())][0] == SOLVED]
            if not solved_by:
                continue
            outcome, res = run(c, s, PAPER, max_it=30)
            if outcome == SOLVED:
                paper = f"slow (converges in {res.refine.iterations})"
            else:
                paper = "reaches another root" if outcome == OTHER_ROOT else "diverges"
            rows.append([s, c, paper, "; ".join(solved_by)])
    write_table(f"{prefix}_gains",
                ["setting", "case", "paper's method with 30 NR iterations", "solved within 10 by"], rows,
                "Cases the paper's method does not solve that another configuration solves")


def table_cases(results, configs, prefix):
    rows = []
    for c in ALL:
        cells = []
        for options in configs:
            outcome, lu = results[(c, "S1 paper", options.name())]
            cells.append(str(lu) if outcome == SOLVED else outcome)
        rows.append([c] + cells)
    write_table(f"{prefix}_cases_S1", ["case"] + [o.name() for o in configs], rows,
                "Setting S1: LU factorizations per case, or the outcome if not solved")


def normalise(name):
    return name.replace(" ", "").lower()


def select(names):
    """The configurations named on the command line, with the paper's method first."""
    by_name = {normalise(o.name()): o for o in COMBINATIONS}
    unknown = [n for n in names if normalise(n) not in by_name]
    if unknown:
        known = ", ".join(o.name() for o in COMBINATIONS)
        raise SystemExit(f"unknown configuration {unknown}. Known: {known}")
    chosen = [by_name[normalise(n)] for n in names]
    return [PAPER] + [o for o in chosen if o != PAPER]


def main(names=None):
    configs = select(names) if names else COMBINATIONS
    prefix = "improvements_selected" if names else "improvements"
    results = run_everything(configs, prefix)
    table_solved(results, configs, prefix)
    table_cost(results, configs, prefix)
    table_gains(results, configs, prefix)
    table_cases(results, configs, prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--configs", nargs="+", help='configurations to run, e.g. PC "OM + PC" richardson')
    main(parser.parse_args().configs)
