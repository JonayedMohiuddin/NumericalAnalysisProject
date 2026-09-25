"""Step rules for following the path: step(h, x, t, dt) -> x at t + dt.

backward_euler (from dynhomotopy) is the paper's rule. The others are the
alternatives we tested. with_corrector() adds a Newton corrector to any of them.
"""

from dynhomotopy.homotopy.integrators import backward_euler, forward_euler, runge_kutta2
from dynhomotopy.linalg import factorize, solve


def backward_euler_chord(h, x, t, dt, iterations=3):
    """BE with extra fixed point iterations that reuse one LU of Gx.

    This is the refinement the paper suggests after eq. (19). In our tests it
    made the start worse, because the iteration drifts when dt is large.
    """
    lu = factorize(h.Gx(x, t + dt), h.counter)
    x_new = x
    for _ in range(iterations):
        x_new = x - dt * lu(h.Gt(x_new))
    return x_new


def runge_kutta4(h, x, t, dt):
    """Classical RK4. Like FE and RK2 it is explicit, and it fails on most cases."""
    k1 = h.dxdt(x, t)
    k2 = h.dxdt(x + dt / 2 * k1, t + dt / 2)
    k3 = h.dxdt(x + dt / 2 * k2, t + dt / 2)
    k4 = h.dxdt(x + dt * k3, t + dt)
    return x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def with_corrector(step):
    """After `step`, take one Newton step on G(x, t + dt) = 0 (predictor-corrector).

    The corrector pulls the point back towards the path, at the cost of one
    more LU per point.
    """
    def corrected(h, x, t, dt):
        x_new = step(h, x, t, dt)
        return x_new - solve(h.Gx(x_new, t + dt), h.G(x_new, t + dt), h.counter)
    return corrected


STEPS = {
    "BE": backward_euler,
    "BE-chord": backward_euler_chord,
    "FE": forward_euler,
    "RK2": runge_kutta2,
    "RK4": runge_kutta4,
}
EXPLICIT = {"FE", "RK2", "RK4"}
ORDER = {"BE": 1, "BE-chord": 1, "FE": 1, "RK2": 2, "RK4": 4}
