import time

import numpy as np

from ..linalg import LUCounter, SingularMatrixError
from ..results import Trajectory
from .fpv import FixedPointHomotopy
from .integrators import SCHEMES


def integrate(problem, x0, K, times, method="BE", first_step="BE"):
    start = time.perf_counter()
    counter = LUCounter()
    h = FixedPointHomotopy(problem, np.asarray(x0), K, counter)
    x = h.x0.copy()
    states, norms = [x], [problem.norm(x)]
    failed = False

    for k in range(len(times) - 1):
        name = first_step if (k == 0 and first_step) else method
        t, dt = times[k], times[k + 1] - times[k]
        try:
            with np.errstate(all="ignore"):
                x = SCHEMES[name](h, x, t, dt)
                norm = problem.norm(x)
        except SingularMatrixError:
            failed = True
            break
        states.append(x)
        norms.append(norm)
        if not np.isfinite(norm):
            failed = True
            break

    return Trajectory(method, K, times[:len(states)], states, norms,
                      counter.factorizations, time.perf_counter() - start, failed)
