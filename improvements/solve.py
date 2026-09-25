import time
from dataclasses import dataclass

import numpy as np

from dynhomotopy.homotopy.integrators import backward_euler
from dynhomotopy.linalg import LUCounter
from dynhomotopy.results import HybridResult
from dynhomotopy.solvers import fast_decoupled_xb, newton_raphson

from .homotopies import HOMOTOPIES
from .multiplier import newton_raphson_om
from .path import follow_path, follow_path_adaptive, follow_path_richardson
from .steps import EXPLICIT, ORDER, STEPS, with_corrector


@dataclass(frozen=True)
class Options:
    """Switches for each change to the paper's method. The defaults give the paper's BE + NR.

    multiplier  optimal multiplier in the final NR
    corrector   Newton corrector after every step
    adaptive    choose the time points on the fly, halving while ||g|| grows
    richardson  choose the time points from a Richardson error estimate,
                with tolerance richardson_tol (cannot be combined with adaptive)
    step        "BE" (paper), "BE-chord", "FE", "RK2" or "RK4"
    homotopy    "fpv" (paper), "scaled" or "newton"

    results/tables/improvements_*.md shows which of these help.
    """

    multiplier: bool = False
    corrector: bool = False
    adaptive: bool = False
    richardson: bool = False
    richardson_tol: float = 0.5
    step: str = "BE"
    homotopy: str = "fpv"

    def __post_init__(self):
        if self.adaptive and self.richardson:
            raise ValueError("choose either adaptive or richardson, not both")

    def name(self):
        parts = [label for label, on in (("OM", self.multiplier), ("PC", self.corrector),
                                         ("adaptive", self.adaptive),
                                         ("richardson", self.richardson)) if on]
        if self.step != "BE":
            parts.append(self.step)
        if self.homotopy != "fpv":
            parts.append(f"{self.homotopy} homotopy")
        return " + ".join(parts) or "paper"


def step_rules(options):
    """The step rule for the first step and for the rest.

    Explicit rules blow up on the first step (Section 4.2.1), so as in the
    paper that step is always taken with BE.
    """
    step = STEPS[options.step]
    first = backward_euler if options.step in EXPLICIT else step
    if options.corrector:
        step, first = with_corrector(step), with_corrector(first)
    return first, step


def solve(problem, x0, K, times, options=Options(), refiner="NR", tol=1e-8, max_it=None,
          record_states=False):
    """Follow the homotopy path, then refine x(1) with NR or FDXB.

    With options.adaptive or options.richardson only times[1] (the first step dt0) is used.
    The multiplier only applies when the refiner is NR.
    """
    start = time.perf_counter()
    h = HOMOTOPIES[options.homotopy](problem, np.asarray(x0), K, LUCounter())
    first, step = step_rules(options)
    if options.adaptive:
        traj = follow_path_adaptive(h, times[1], step, first)
    elif options.richardson:
        traj = follow_path_richardson(h, times[1], step, first, ORDER[options.step],
                                      options.richardson_tol)
    else:
        traj = follow_path(h, times, step, first)

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
