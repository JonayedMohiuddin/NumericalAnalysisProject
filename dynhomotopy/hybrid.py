"""The full method of the paper: dynamic homotopy, then a classical refiner.

    x0 (flat start) --[FE | RK2 | BE along gamma(t)]--> x^(0) = x(t)|_{t=1}
                    --[NR | FDXB]--> x*  with ||g(x*)||_inf < tol
"""

from __future__ import annotations

import time

import numpy as np

from .homotopy import integrate
from .problem import NonlinearProblem
from .results import HybridResult
from .solvers import fast_decoupled_xb, newton_raphson


def solve_hybrid(problem: NonlinearProblem, x0: np.ndarray, K: float, times: list[float],
                 method: str = "BE", refiner: str = "NR", tol: float = 1e-8,
                 max_it: int | None = None, first_step: str | None = "BE",
                 record_states: bool = False) -> HybridResult:
    """Run the dynamic homotopy over ``times`` and refine x(1) with NR or FDXB."""
    start = time.perf_counter()
    traj = integrate(problem, x0, K, times, method, first_step)
    refine = None
    if not traj.failed:
        with np.errstate(all="ignore"):
            if refiner == "NR":
                refine = newton_raphson(problem, traj.x_final, tol, max_it or 10, record_states)
            elif refiner == "FDXB":
                refine = fast_decoupled_xb(problem, traj.x_final, tol, max_it or 100, record_states)
            else:
                raise ValueError(f"unknown refiner {refiner!r}")
    return HybridResult(traj, refine, time.perf_counter() - start)
