"""Testing modifications of the method that are suggested by the paper or its references.

Every variant is run on the paper's 9 cases plus 25 more cases from the same
Zenodo records, under 5 settings of (dt0, K). A run succeeds if the refining
NR converges to 1e-8 within 10 iterations, the same rule as in the paper.
Cost is counted in LU factorizations (homotopy + NR), which unlike run time
does not depend on the machine.

Raw results go to results/tables/modifications_raw.csv.
"""

import csv
import time

import numpy as np

from dynhomotopy.datasets import ALL, EXTRA, load_problem
from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.solvers import newton_raphson

from .common import TABLES, write_table

CASES = ALL + EXTRA

# name: (dt0, K). The path is {0, dt0, 1} as in Sections 4.5 and 4.7.
SETTINGS = {
    "S1 paper (dt0=0.005, K=1e-4)": (0.005, 1e-4),
    "S2 Sec 4.4 (dt0=0.05, K=1e-3)": (0.05, 1e-3),
    "S3 Sec 4.4 (dt0=0.1, K=2e-3)": (0.1, 2e-3),
    "S4 weak K (dt0=0.005, K=1e-5)": (0.005, 1e-5),
    "S5 strong K (dt0=0.005, K=1e-3)": (0.005, 1e-3),
}

VARIANTS = {
    "BE (paper)": {},
    "BE-chord": dict(method="BE-chord", first_step="BE-chord"),
    "BE-PC": dict(method="BE-PC", first_step="BE-PC"),
    "BE + OM": dict(multiplier=True),
    "BE-PC + OM": dict(method="BE-PC", first_step="BE-PC", multiplier=True),
    "RK4": dict(method="RK4"),
    "adaptive BE": dict(adaptive=True),
    "adaptive BE-PC": dict(adaptive=True, method="BE-PC"),
    "Newton homotopy": dict(homotopy="newton"),
    "scaled K": dict(homotopy="scaled"),
}


def run_all_combinations():
    rows = []
    for name in CASES:
        pf = load_problem(name)
        for label, kwargs in [("NR flat", {}), ("NR flat + OM", dict(multiplier=True))]:
            start = time.perf_counter()
            r = newton_raphson(pf, pf.flat_start(), **kwargs)
            rows.append([name, "-", label, r.converged, r.factorizations, r.iterations,
                         time.perf_counter() - start])
        for setting, (dt0, K) in SETTINGS.items():
            for variant, kwargs in VARIANTS.items():
                r = solve_hybrid(pf, pf.flat_start(), K, [0.0, dt0, 1.0], **kwargs)
                iterations = r.refine.iterations if r.refine else -1
                rows.append([name, setting, variant, r.converged, r.factorizations, iterations, r.time])
        print(f"{name} done", flush=True)

    with open(TABLES / "modifications_raw.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["case", "setting", "variant", "converged", "lu", "nr_iterations", "time"])
        writer.writerows(rows)
    return rows


def load_raw():
    with open(TABLES / "modifications_raw.csv") as fh:
        return [[r["case"], r["setting"], r["variant"], r["converged"] == "True", int(r["lu"]),
                 int(r["nr_iterations"]), float(r["time"])] for r in csv.DictReader(fh)]


def summarise(rows):
    result = {(r[0], r[1], r[2]): r for r in rows}
    settings = list(SETTINGS)
    variants = list(VARIANTS)

    # 1. number of cases solved
    table = []
    for v in variants:
        counts = [sum(result[(c, s, v)][3] for c in CASES) for s in settings]
        table.append([v] + counts + [sum(counts)])
    for v in ["NR flat", "NR flat + OM"]:
        solved = sum(result[(c, "-", v)][3] for c in CASES)
        table.append([v] + [solved] * len(settings) + ["-"])
    write_table("modifications_solved", ["variant"] + [s.split(" (")[0] for s in settings] + ["total"],
                table, f"Cases solved out of {len(CASES)} (NR converges within 10 iterations)")

    # 2. cost on the cases the paper's method already solves
    table = []
    base = "BE (paper)"
    for v in variants:
        row = [v]
        for s in settings:
            both = [c for c in CASES if result[(c, s, base)][3] and result[(c, s, v)][3]]
            if not both:
                row.append("-")
                continue
            lu_v = np.mean([result[(c, s, v)][4] for c in both])
            lu_b = np.mean([result[(c, s, base)][4] for c in both])
            row.append(f"{lu_v:.2f} vs {lu_b:.2f} ({len(both)})")
        table.append(row)
    write_table("modifications_cost", ["variant"] + [s.split(" (")[0] for s in settings], table,
                "Mean LU factorizations on cases solved by both the variant and the paper's BE "
                "(variant vs paper, number of cases)")

    # 3. which cases each variant solves in the paper's setting
    s1 = settings[0]
    table = []
    for c in CASES:
        row = [c, "yes" if result[(c, "-", "NR flat")][3] else "no"]
        row += [f"{result[(c, s1, v)][4]}" if result[(c, s1, v)][3] else "fail" for v in variants]
        table.append(row)
    write_table("modifications_cases_S1", ["case", "NR flat"] + variants, table,
                "Setting S1: total LU factorizations per case, or fail")


def main(rerun=True):
    rows = run_all_combinations() if rerun else load_raw()
    summarise(rows)


if __name__ == "__main__":
    main()
