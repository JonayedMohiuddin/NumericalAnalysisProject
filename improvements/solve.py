import time
from dataclasses import dataclass

import numpy as np

from dynhomotopy.homotopy.integrators import backward_euler
from dynhomotopy.results import HybridResult
from dynhomotopy.solvers import fast_decoupled_xb, newton_raphson

from .corrector import backward_euler_corrected
from .multiplier import newton_raphson_om
from .path import follow_path, follow_path_adaptive


@dataclass(frozen=True)
class Options:
    """Switches for the three improvements. All off gives the paper's BE + NR."""

    multiplier: bool = False  # optimal multiplier in the final NR
    corrector: bool = False   # Newton corrector after every BE step
    adaptive: bool = False    # choose the time points on the fly

    def name(self):
        parts = [p for p, on in (("OM", self.multiplier), ("PC", self.corrector),
                                 ("adaptive", self.adaptive)) if on]
        return " + ".join(parts) or "paper"


def solve(problem, x0, K, times, options=Options(), refiner="NR", tol=1e-8, max_it=None,
          record_states=False):
    """Backward Euler along the homotopy path, then NR or FDXB from x(1).

    With options.adaptive only times[1] (the first step dt0) is used.
    The multiplier only applies when the refiner is NR.
    """
    start = time.perf_counter()
    step = backward_euler_corrected if options.corrector else backward_euler
    if options.adaptive:
        traj = follow_path_adaptive(problem, x0, K, times[1], step)
    else:
        traj = follow_path(problem, x0, K, times, step)

    refine = None
    if not traj.failed:
        with np.errstate(all="ignore"):
            if refiner == "NR":
                newton = newton_raphson_om if options.multiplier else newton_raphson
                refine = newton(problem, traj.x_final, tol, max_it or 10, record_states)
            elif refiner == "FDXB":
                refine = fast_decoupled_xb(problem, traj.x_final, tol, max_it or 100, record_states)
            else:
                raise ValueError(f"unknown refiner {refiner}")
    return HybridResult(traj, refine, time.perf_counter() - start)
