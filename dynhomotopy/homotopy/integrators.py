"""One step of each integration scheme, from (x, t) to t + dt."""

from ..linalg import factorize, solve


def forward_euler(h, x, t, dt):
    # eq. (16)
    return x + dt * h.dxdt(x, t)


def backward_euler(h, x, t, dt, fpi_iterations=1):
    # eq. (19): the implicit equation (17) is solved by fixed point iteration
    # starting from x. The paper uses a single iteration.
    x_new = x
    for _ in range(fpi_iterations):
        x_new = x + dt * h.dxdt(x_new, t + dt)
    return x_new


def backward_euler_chord(h, x, t, dt, fpi_iterations=3):
    # The variant suggested after eq. (19): keep the LU of Gx(x, t + dt) and
    # only update Gt in the extra fixed point iterations.
    lu = factorize(h.Gx(x, t + dt), h.counter)
    x_new = x
    for _ in range(fpi_iterations):
        x_new = x - dt * lu(h.Gt(x_new))
    return x_new


def backward_euler_corrected(h, x, t, dt):
    # BE as a predictor, then one Newton step on G(x, t + dt) = 0 as a corrector.
    x_new = backward_euler(h, x, t, dt)
    return x_new - solve(h.Gx(x_new, t + dt), h.G(x_new, t + dt), h.counter)


def runge_kutta2(h, x, t, dt):
    # eqs. (20)-(21)
    k1 = h.dxdt(x, t)
    k2 = h.dxdt(x + dt * k1, t + dt)
    return x + dt / 2 * (k1 + k2)


def runge_kutta4(h, x, t, dt):
    k1 = h.dxdt(x, t)
    k2 = h.dxdt(x + dt / 2 * k1, t + dt / 2)
    k3 = h.dxdt(x + dt / 2 * k2, t + dt / 2)
    k4 = h.dxdt(x + dt * k3, t + dt)
    return x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def linear_first_step(h, x, t, dt):
    # eqs. (26)-(27): linearise G around x at t + dt. For the first step this
    # gives the same point as backward_euler.
    return x - solve(h.Gx(x, t + dt), h.G(x, t + dt), h.counter)


SCHEMES = {
    "FE": forward_euler,
    "BE": backward_euler,
    "BE-chord": backward_euler_chord,
    "BE-PC": backward_euler_corrected,
    "RK2": runge_kutta2,
    "RK4": runge_kutta4,
    "LIN": linear_first_step,
}
