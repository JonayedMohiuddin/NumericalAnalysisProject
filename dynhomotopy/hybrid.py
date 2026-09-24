import time

import numpy as np

from .homotopy import integrate
from .results import HybridResult
from .solvers import fast_decoupled_xb, newton_raphson


def solve_hybrid(problem, x0, K, times, method="BE", refiner="NR", tol=1e-8,
                 max_it=None, first_step="BE", record_states=False):
    """The method of the paper: follow the homotopy to t = 1, then refine x(1) with NR or FDXB."""
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
                raise ValueError(f"unknown refiner {refiner}")
    return HybridResult(traj, refine, time.perf_counter() - start)
