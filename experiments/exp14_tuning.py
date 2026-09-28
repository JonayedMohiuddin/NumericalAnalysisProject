"""Proposal: tuning delta and dt0 by golden-section search, and predicting delta from the spectrum.

1. For every case, golden-section search over log10(delta) (dt0 = 0.005),
   then over log10(dt0), minimising the total LU factorizations of the
   paper's method (improvements/tuning.py).
2. Two rules for choosing delta without a search are fitted on the paper's
   9 cases and tested on the 25 other cases:
     fixed rule     delta = geometric mean of the tuned deltas
     spectral rule  delta = c * lambda_min(J(x0)), with c the geometric mean
                    of tuned delta / lambda_min
   The spectral rule costs one extra LU (inverse iteration).
3. Both rules, the paper's delta = 0.02 and the per-case search are compared
   by cases solved (reference operating point) and LU factorizations.
"""

import numpy as np

from dynhomotopy.datasets import ALL as PAPER_CASES
from improvements.spectrum import inverse_iteration
from improvements.tuning import tune
from improvements import solve
from improvements.cases import EXTRA, load, on_reference

from .common import write_table

DT0 = 0.005
PAPER_DELTA = 0.02


def evaluate(pf, delta, dt0=DT0):
    """(solved, LU factorizations) of the paper's method with this delta and dt0."""
    res = solve(pf, pf.flat_start(), delta * dt0, [0.0, dt0, 1.0])
    solved = res.converged and on_reference(pf, res.refine.x)
    return solved, res.factorizations


def geometric_mean(values):
    return float(np.exp(np.mean(np.log(values))))


def main():
    per_case = {}
    rows = []
    for name in PAPER_CASES + EXTRA:
        pf = load(name)
        lam_min, _ = inverse_iteration(pf.jacobian(pf.flat_start()).tocsc())
        delta, dt0, evaluations = tune(pf)
        solved, lus = evaluate(pf, delta, dt0)
        per_case[name] = dict(pf=pf, lam_min=lam_min, delta=delta, dt0=dt0)
        rows.append([name, f"{lam_min:.3g}", f"{delta:.3g}", f"{dt0:.3g}", evaluations,
                     "yes" if solved else "no", lus])
        print(f"{name} tuned", flush=True)
    write_table("tuning_per_case", ["case", "lambda_min", "tuned delta", "tuned dt0", "search evaluations",
                                    "solved at tuned point", "LUs at tuned point"], rows,
                "Golden-section tuning of delta, then dt0, for each case")

    fixed_delta = geometric_mean([per_case[c]["delta"] for c in PAPER_CASES])
    c_spectral = geometric_mean([per_case[c]["delta"] / per_case[c]["lam_min"] for c in PAPER_CASES])

    rules = {
        "paper, delta = 0.02": lambda c: (PAPER_DELTA, DT0, 0),
        f"fixed rule, delta = {fixed_delta:.3g}": lambda c: (fixed_delta, DT0, 0),
        f"spectral rule, delta = {c_spectral:.3g} lambda_min": lambda c: (c_spectral * per_case[c]["lam_min"], DT0, 1),
        "golden search per case": lambda c: (per_case[c]["delta"], per_case[c]["dt0"], 0),
    }
    rows = []
    for rule, choose in rules.items():
        row = [rule]
        for group, cases in (("paper's 9 (training)", PAPER_CASES), ("25 others (test)", EXTRA)):
            outcomes = []
            for c in cases:
                delta, dt0, extra_lus = choose(c)
                solved, lus = evaluate(per_case[c]["pf"], delta, dt0)
                outcomes.append((solved, lus + extra_lus))
            solved = [lu for ok, lu in outcomes if ok]
            row += [f"{len(solved)} of {len(cases)}", f"{np.mean(solved):.2f}" if solved else "-"]
        rows.append(row)
    write_table("tuning_rules", ["how delta is chosen", "solved, training", "mean LUs, training",
                                 "solved, test", "mean LUs, test"], rows,
                "Choosing delta: rules fitted on the paper's 9 cases, tested on the other 25. "
                "The golden search row excludes the cost of the search itself")


if __name__ == "__main__":
    main()
