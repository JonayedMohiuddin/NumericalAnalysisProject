"""Backward Euler with a Newton corrector (predictor-corrector, Section 2.1.2 and refs. [36, 37])."""

from dynhomotopy.homotopy.integrators import backward_euler
from dynhomotopy.linalg import solve


def backward_euler_corrected(h, x, t, dt):
    """One BE step (the paper's predictor), then one Newton step on G(x, t + dt) = 0.

    The BE step only follows the path approximately. The corrector pulls the
    point back towards the path, at the cost of one more LU per point.
    """
    x_new = backward_euler(h, x, t, dt)
    return x_new - solve(h.Gx(x_new, t + dt), h.G(x_new, t + dt), h.counter)
