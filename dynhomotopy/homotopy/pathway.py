"""Discrete time instants t_k of the homotopy pathway gamma(t), 0 <= t <= 1 (Section 4.2).

    t_0 = 0,  t_1 = dt0,  t_{k+1} = t_k + dt,  with the last step shortened so t_N = 1.
"""

from __future__ import annotations


def constant_step(dt0: float, dt: float, t2: float | None = None) -> list[float]:
    """[0, dt0, (t2,) dt0 + dt, ...,  1].

    ``t2`` inserts an extra small step, e.g. t2 = 2 dt0 in Tables 2-4.
    Example (Section 4.2): dt0 = 0.005, dt = 1.0 gives [0, 0.005, 1.0]
    (dt is shortened to 0.995 so the pathway ends at t = 1).
    """
    times = [0.0, dt0]
    if t2 is not None:
        times.append(t2)
    while times[-1] < 1.0 - 1e-12:
        times.append(min(round(times[-1] + dt, 12), 1.0))
    return times


def explicit(*times: float) -> list[float]:
    """A pathway given point by point, e.g. explicit(0, 0.05, 0.1, 0.2, 1)."""
    out = [float(t) for t in times]
    if out[0] != 0.0 or abs(out[-1] - 1.0) > 1e-12 or any(b <= a for a, b in zip(out, out[1:])):
        raise ValueError(f"pathway must increase from 0 to 1, got {out}")
    return out
