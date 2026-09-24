import time

import numpy as np

from ..linalg import LUCounter, SingularMatrixError
from ..results import Trajectory
from .fpv import FixedPointHomotopy, NewtonHomotopy, jacobian_diagonal_scale
from .integrators import SCHEMES


def make_homotopy(kind, problem, x0, K, counter):
    if kind == "fpv":
        return FixedPointHomotopy(problem, x0, K, counter)
    if kind == "scaled":
        return FixedPointHomotopy(problem, x0, K, counter, jacobian_diagonal_scale(problem, x0))
    if kind == "newton":
        return NewtonHomotopy(problem, x0, counter=counter)
    raise ValueError(f"unknown homotopy {kind}")


def take_step(scheme, h, x, t, dt, problem):
    """One step; returns (x_new, norm), with norm = inf if the step broke down."""
    try:
        with np.errstate(all="ignore"):
            x_new = SCHEMES[scheme](h, x, t, dt)
            return x_new, problem.norm(x_new)
    except SingularMatrixError:
        return None, np.inf


def integrate(problem, x0, K, times, method="BE", first_step="BE", homotopy="fpv"):
    """Follow the homotopy path from x(0) = x0 over the given time points.

    As in Section 4.3, the first step to t1 is taken with BE (explicit schemes
    blow up there) and `method` is used for the rest. Pass first_step=None to
    use `method` for the first step too.
    """
    start = time.perf_counter()
    counter = LUCounter()
    h = make_homotopy(homotopy, problem, np.asarray(x0), K, counter)
    x = h.x0.copy()
    states, norms = [x], [problem.norm(x)]
    failed = False

    for k in range(len(times) - 1):
        scheme = first_step if (k == 0 and first_step) else method
        x, norm = take_step(scheme, h, x, times[k], times[k + 1] - times[k], problem)
        if not np.isfinite(norm):
            failed = True
            break
        states.append(x)
        norms.append(norm)

    return Trajectory(method, K, times[:len(states)], states, norms,
                      counter.factorizations, time.perf_counter() - start, failed)


def integrate_adaptive(problem, x0, K, dt0, method="BE", homotopy="fpv", max_halvings=6):
    """Like integrate(), but the points after t1 are chosen on the fly.

    After t1 = dt0 it first tries to jump straight to t = 1. A step is accepted
    if ||g|| went down, otherwise the step is halved. After an accepted step the
    next step is doubled. Rejected steps still count as LU factorizations.
    """
    start = time.perf_counter()
    counter = LUCounter()
    h = make_homotopy(homotopy, problem, np.asarray(x0), K, counter)
    x = h.x0.copy()
    times, states, norms = [0.0], [x], [problem.norm(x)]

    x, norm = take_step("BE", h, x, 0.0, dt0, problem)
    failed = not np.isfinite(norm)
    if not failed:
        times.append(dt0)
        states.append(x)
        norms.append(norm)

    dt = 1 - dt0
    while not failed and times[-1] < 1 - 1e-12:
        t = times[-1]
        for _ in range(max_halvings + 1):
            dt = min(dt, 1 - t)
            x_new, norm = take_step(method, h, x, t, dt, problem)
            if norm < norms[-1]:
                break
            dt /= 2
        else:
            failed = True
            break
        x = x_new
        times.append(t + dt)
        states.append(x)
        norms.append(norm)
        dt *= 2

    return Trajectory(f"{method}-adaptive", K, times, states, norms,
                      counter.factorizations, time.perf_counter() - start, failed)
