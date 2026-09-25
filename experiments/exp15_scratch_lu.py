"""Proposal: NR and FDXB with our own Gauss elimination / LU with partial pivoting.

The paper's method is run twice on the smaller systems, once with SciPy's
sparse LU (SuperLU) and once with the dense LU written in extensions/gauss.py.
The iterations and the solution should be the same. The run times show why
a sparse LU is needed for large grids.
"""

import time

import numpy as np

from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.problem import TutorialProblem
from extensions import scratch
from extensions.gauss import gauss_solve, lu_factor, lu_solve
from improvements.cases import load

from .common import write_table

CASES = ["case69limit", "case141limit", "case_ACTIVSg500limit", "case_ACTIVSg2000limit"]
PATH = [0.0, 0.005, 1.0]
K = 1e-4


def compare(pf, x0, refiner, tol=1e-8):
    start = time.perf_counter()
    ref = solve_hybrid(pf, x0, K, PATH, refiner=refiner, tol=tol, max_it=50)
    t_scipy = time.perf_counter() - start

    start = time.perf_counter()
    x1, _ = scratch.backward_euler_path(pf, x0, K, PATH)
    if refiner == "NR":
        ours = scratch.newton_raphson(pf, x1, tol, max_it=50)
    else:
        ours = scratch.fast_decoupled_xb(pf, x1, tol)
    t_ours = time.perf_counter() - start

    diff = np.max(np.abs(ours.x - ref.refine.x))
    return [ref.refine.iterations, ours.iterations, f"{diff:.1e}", f"{t_scipy:.3f}", f"{t_ours:.3f}"]


def gauss_versus_lu(n=400, right_hand_sides=20):
    """Many solves with one matrix, as in FDXB: Gauss redoes the elimination, LU does not."""
    rng = np.random.default_rng(0)
    a = rng.standard_normal((n, n)) + n * np.eye(n)
    bs = rng.standard_normal((right_hand_sides, n))
    start = time.perf_counter()
    for b in bs:
        gauss_solve(a, b)
    t_gauss = time.perf_counter() - start
    start = time.perf_counter()
    lu, perm = lu_factor(a)
    for b in bs:
        lu_solve(lu, perm, b)
    t_lu = time.perf_counter() - start
    return t_gauss, t_lu


def main(cases=CASES):
    rows = []
    tutorial = TutorialProblem()
    rows.append(["tutorial (Sec. 3.3)", "NR", *compare(tutorial, tutorial.x0, "NR", tol=1e-5)])
    for name in cases:
        pf = load(name)
        for refiner in ["NR", "FDXB"]:
            rows.append([name, refiner, *compare(pf, pf.flat_start(), refiner)])
            print(name, refiner, "done", flush=True)
    write_table("scratch_lu", ["case", "refiner", "iterations (SuperLU)", "iterations (our LU)",
                               "max |x difference|", "time SuperLU (s)", "time our LU (s)"], rows,
                "BE path + NR / FDXB with SciPy's sparse LU and with our dense LU")

    t_gauss, t_lu = gauss_versus_lu()
    write_table("gauss_vs_lu", ["method", "time for 20 solves with one 400 x 400 matrix (s)"],
                [["Gauss elimination for every right-hand side", f"{t_gauss:.3f}"],
                 ["one LU, then 20 triangular solves", f"{t_lu:.3f}"]],
                "Why FDXB factorises B' and B'' once")


if __name__ == "__main__":
    main()
