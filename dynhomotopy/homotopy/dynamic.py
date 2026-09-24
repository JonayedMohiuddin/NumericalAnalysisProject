"""Integrate the dynamic-homotopy IVP along a pathway (Section 3 and 4.3).

Following Section 4.3, the first point x(t1) is always computed with the
implicit rule (BE or, equivalently, the linear approximation (27)), because
explicit schemes blow up there when ||g(x0)|| >> 1 (Section 4.2.1). The
chosen scheme (FE, RK2 or BE) then advances the remaining points.
"""

from __future__ import annotations

import time

import numpy as np

from ..linalg import LUCounter, SingularMatrixError
from ..problem import NonlinearProblem
from ..results import Trajectory
from .fpv import FixedPointHomotopy
from .integrators import SCHEMES


def integrate(problem: NonlinearProblem, x0: np.ndarray, K: float, times: list[float],
              method: str = "BE", first_step: str | None = "BE") -> Trajectory:
    """Discrete pathway {x(t_k)} from the IC x(0) = x0.

    method      scheme for t_k > t_1: "FE", "BE" or "RK2".
    first_step  scheme for [0, t_1]: "BE" (default), "LIN" (eq. 27), or None
                to use ``method`` there too (the Section 4.2.1 experiment).
    """
    start = time.perf_counter()
    counter = LUCounter()
    h = FixedPointHomotopy(problem, x0, K, counter)
    x = h.x0.copy()
    states, norms = [x.copy()], [problem.norm(x)]
    failed = False
    for k in range(len(times) - 1):
        scheme = SCHEMES[first_step if (k == 0 and first_step) else method]
        t, dt = times[k], times[k + 1] - times[k]
        try:
            with np.errstate(all="ignore"):
                x = scheme(h, x, t, dt)
                nrm = problem.norm(x)
        except SingularMatrixError:
            failed = True
            break
        states.append(x.copy())
        norms.append(nrm)
        if not np.isfinite(nrm):
            failed = True
            break
    return Trajectory(method, K, list(times[: len(states)]), states, norms,
                      counter.factorizations, time.perf_counter() - start, failed)
