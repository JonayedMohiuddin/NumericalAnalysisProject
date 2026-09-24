"""Following the homotopy path with a fixed or an adaptive set of time points.

`step` is a function step(h, x, t, dt) -> x_new, either the paper's
backward_euler or backward_euler_corrected.
"""

import time

import numpy as np

from dynhomotopy.homotopy.fpv import FixedPointHomotopy
from dynhomotopy.linalg import LUCounter, SingularMatrixError
from dynhomotopy.results import Trajectory


def try_step(step, h, x, t, dt):
    """Returns (x_new, ||g(x_new)||). The norm is inf if the step broke down."""
    try:
        with np.errstate(all="ignore"):
            x_new = step(h, x, t, dt)
            return x_new, h.problem.norm(x_new)
    except SingularMatrixError:
        return None, np.inf


def follow_path(problem, x0, K, times, step):
    """Take one step between each pair of consecutive time points."""
    start = time.perf_counter()
    h = FixedPointHomotopy(problem, np.asarray(x0), K, LUCounter())
    states, norms = [h.x0.copy()], [problem.norm(h.x0)]
    failed = False
    for t, t_next in zip(times, times[1:]):
        x, norm = try_step(step, h, states[-1], t, t_next - t)
        if not np.isfinite(norm):
            failed = True
            break
        states.append(x)
        norms.append(norm)
    return Trajectory("BE", K, list(times[:len(states)]), states, norms,
                      h.counter.factorizations, time.perf_counter() - start, failed)


def follow_path_adaptive(problem, x0, K, dt0, step, max_halvings=6):
    """Choose the time points while following the path.

    The first step goes to t1 = dt0 as in the paper. After that we try to jump
    straight to t = 1. A step is kept if ||g|| went down, otherwise it is
    halved and tried again. After a kept step the next one is doubled.
    Rejected steps still count as LU factorizations.
    """
    start = time.perf_counter()
    h = FixedPointHomotopy(problem, np.asarray(x0), K, LUCounter())
    times, states, norms = [0.0], [h.x0.copy()], [problem.norm(h.x0)]
    failed = False
    dt = dt0
    while times[-1] < 1 - 1e-12:
        t = times[-1]
        for _ in range(max_halvings + 1):
            dt = min(dt, 1 - t)
            x, norm = try_step(step, h, states[-1], t, dt)
            if norm < norms[-1] or (t == 0 and np.isfinite(norm)):
                break
            dt /= 2
        else:
            failed = True
            break
        times.append(t + dt)
        states.append(x)
        norms.append(norm)
        dt = 1 - dt0 if t == 0 else 2 * dt
    return Trajectory("BE-adaptive", K, times, states, norms,
                      h.counter.factorizations, time.perf_counter() - start, failed)
