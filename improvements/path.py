"""Following the homotopy path h with a fixed or an adaptive set of time points.

`step` is a function step(h, x, t, dt) -> x_new from improvements.steps,
or the paper's backward_euler.
"""

import time

import numpy as np

from dynhomotopy.linalg import SingularMatrixError
from dynhomotopy.results import Trajectory


def try_step(step, h, x, t, dt):
    """Returns (x_new, ||g(x_new)||). The norm is inf if the step broke down."""
    try:
        with np.errstate(all="ignore"):
            x_new = step(h, x, t, dt)
            return x_new, h.problem.norm(x_new)
    except SingularMatrixError:
        return None, np.inf


def follow_path(h, times, step, first_step):
    """Take one step between each pair of consecutive time points."""
    start = time.perf_counter()
    states, norms = [h.x0.copy()], [h.problem.norm(h.x0)]
    failed = False
    for k, (t, t_next) in enumerate(zip(times, times[1:])):
        x, norm = try_step(first_step if k == 0 else step, h, states[-1], t, t_next - t)
        if not np.isfinite(norm):
            failed = True
            break
        states.append(x)
        norms.append(norm)
    return Trajectory("fixed", h.K, list(times[:len(states)]), states, norms,
                      h.counter.factorizations, time.perf_counter() - start, failed)


def follow_path_adaptive(h, dt0, step, first_step, max_halvings=6):
    """Choose the time points while following the path.

    The first step goes to t1 = dt0 as in the paper. After that we try to jump
    straight to t = 1. A step is kept if ||g|| went down, otherwise it is
    halved and tried again. After a kept step the next one is doubled.
    Rejected steps still count as LU factorizations.
    """
    start = time.perf_counter()
    times, states, norms = [0.0], [h.x0.copy()], [h.problem.norm(h.x0)]
    failed = False
    dt = dt0
    while times[-1] < 1 - 1e-12:
        t = times[-1]
        rule = first_step if t == 0 else step
        for _ in range(max_halvings + 1):
            dt = min(dt, 1 - t)
            x, norm = try_step(rule, h, states[-1], t, dt)
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
    return Trajectory("adaptive", h.K, times, states, norms,
                      h.counter.factorizations, time.perf_counter() - start, failed)


def follow_path_richardson(h, dt0, step, first_step, order, tol=0.5, max_lu=150):
    """Choose the time points from a Richardson estimate of the local error.

    Each attempt compares one step of size dt with two steps of size dt / 2.
    For a method of order p their difference, divided by 2^p - 1, estimates
    the error of the two half steps. The step is kept if that error is below
    tol (radians and p.u.), and the kept point is the Richardson extrapolation
    of the two results. The next dt is scaled by (tol / error)^(1 / (p + 1)).
    Every attempt costs three steps. A path that needs more than max_lu LU
    factorizations is stopped and counted as failed.
    """
    start = time.perf_counter()
    times, states, norms = [0.0], [h.x0.copy()], [h.problem.norm(h.x0)]
    x, norm = try_step(first_step, h, h.x0, 0.0, dt0)
    failed = not np.isfinite(norm)
    if not failed:
        times.append(dt0)
        states.append(x)
        norms.append(norm)

    dt = 1 - dt0
    while not failed and times[-1] < 1 - 1e-12:
        if h.counter.factorizations > max_lu:
            failed = True
            break
        t, x = times[-1], states[-1]
        dt = min(dt, 1 - t)
        full, n_full = try_step(step, h, x, t, dt)
        half, n_half = try_step(step, h, x, t, dt / 2)
        both, n_both = (None, np.inf) if half is None else try_step(step, h, half, t + dt / 2, dt / 2)
        if not np.isfinite(n_full + n_half + n_both):
            dt /= 4
            continue
        error = np.linalg.norm(both - full, np.inf) / (2**order - 1)
        if error <= tol:
            x_new = both + (both - full) / (2**order - 1)
            n_new = h.problem.norm(x_new)
            if not np.isfinite(n_new):
                x_new, n_new = both, n_both
            times.append(t + dt)
            states.append(x_new)
            norms.append(n_new)
        dt *= min(4.0, max(0.1, 0.9 * (tol / max(error, 1e-12)) ** (1 / (order + 1))))

    return Trajectory("richardson", h.K, times, states, norms,
                      h.counter.factorizations, time.perf_counter() - start, failed)
