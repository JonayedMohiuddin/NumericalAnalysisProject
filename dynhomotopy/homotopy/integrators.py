"""One step of each integration scheme, from (x, t) to t + dt."""

from ..linalg import solve


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


def runge_kutta2(h, x, t, dt):
    # eqs. (20)-(21)
    k1 = h.dxdt(x, t)
    k2 = h.dxdt(x + dt * k1, t + dt)
    return x + dt / 2 * (k1 + k2)


def linear_first_step(h, x, t, dt):
    # eqs. (26)-(27): linearise G around x at t + dt. For the first step this
    # gives the same point as backward_euler.
    return x - solve(h.Gx(x, t + dt), h.G(x, t + dt), h.counter)


SCHEMES = {
    "FE": forward_euler,
    "BE": backward_euler,
    "RK2": runge_kutta2,
    "LIN": linear_first_step,
}
