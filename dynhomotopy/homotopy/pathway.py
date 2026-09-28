"""Time points t_k of the homotopy path, from 0 to 1."""


def constant_step(dt0, dt, t2=None):
    """
    The last step is shortened so the path ends at 1
    constant_step(0.005, 1.0) gives [0, 0.005, 1.0].
    """
    times = [0.0, dt0]
    if t2 is not None:
        times.append(t2)
    while times[-1] < 1 - 1e-12:
        times.append(min(round(times[-1] + dt, 12), 1.0))
    return times


def explicit(*times):
    times = [float(t) for t in times]
    increasing = all(b > a for a, b in zip(times, times[1:]))
    if times[0] != 0 or abs(times[-1] - 1) > 1e-12 or not increasing:
        raise ValueError(f"path must increase from 0 to 1, got {times}")
    return times
