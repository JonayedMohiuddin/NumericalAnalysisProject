import time

import numpy as np

from .homotopy import integrate
from .homotopy.dynamic import integrate_adaptive
from .results import HybridResult
from .solvers import fast_decoupled_xb, newton_raphson


def solve_hybrid(problem, x0, K, times, method="BE", refiner="NR", tol=1e-8,
                 max_it=None, first_step="BE", record_states=False,
                 homotopy="fpv", adaptive=False, multiplier=False):
    """The method of the paper: follow the homotopy to t = 1, then refine x(1) with NR or FDXB.

    The defaults are the paper's method. The other options switch on the
    modifications compared in experiments/exp11_modifications.py:
    homotopy="scaled" or "newton", adaptive=True (only times[1] is used, the
    rest of the path is chosen on the fly) and multiplier=True (NR with
    Iwamoto's optimal multiplier).
    """
    start = time.perf_counter()
    if adaptive:
        traj = integrate_adaptive(problem, x0, K, times[1], method, homotopy)
    else:
        traj = integrate(problem, x0, K, times, method, first_step, homotopy)
    refine = None
    if not traj.failed:
        with np.errstate(all="ignore"):
            if refiner == "NR":
                refine = newton_raphson(problem, traj.x_final, tol, max_it or 10, record_states, multiplier)
            elif refiner == "FDXB":
                refine = fast_decoupled_xb(problem, traj.x_final, tol, max_it or 100, record_states)
            else:
                raise ValueError(f"unknown refiner {refiner}")
    return HybridResult(traj, refine, time.perf_counter() - start)
