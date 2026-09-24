"""One-step integration schemes for the dynamic-homotopy ODE (Section 3).

Every scheme advances x_k at t_k to x_{k+1} at t_{k+1} = t_k + dt.

    FE  (16):  x_{k+1} = x_k + dt Phi(x_k, t_k)
    BE  (19):  x_{k+1} = x_k + dt Phi(x^_{k+1}, t_{k+1}),  x^_{k+1} = x_k,
               i.e. one fixed-point iteration of the implicit rule (17)
    RK2 (20-21): K1 = Phi(x_k, t_k),  K2 = Phi(x_k + dt K1, t_{k+1}),
               x_{k+1} = x_k + dt/2 (K1 + K2)

For the first step from t = 0 the paper uses either BE or the linear
approximation of the static homotopy (26)-(27); both give

    dx(0) = -[J(x0) + (1 - dt0)/dt0 K I]^{-1} g(x0)   (~ -[J + delta I]^{-1} g, eq. 25)
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from ..linalg import solve
from .fpv import FixedPointHomotopy

Step = Callable[[FixedPointHomotopy, np.ndarray, float, float], np.ndarray]


def forward_euler(h: FixedPointHomotopy, x: np.ndarray, t: float, dt: float) -> np.ndarray:
    return x + dt * h.phi(x, t)


def backward_euler(h: FixedPointHomotopy, x: np.ndarray, t: float, dt: float,
                   fpi_iterations: int = 1) -> np.ndarray:
    """Implicit Euler solved approximately by fixed-point iteration (eq. 19).

    The paper uses a single iteration seeded with x^_{k+1} = x_k, costing one
    LU factorisation per pathway point.
    """
    x_hat = x
    for _ in range(fpi_iterations):
        x_hat = x + dt * h.phi(x_hat, t + dt)
    return x_hat


def runge_kutta2(h: FixedPointHomotopy, x: np.ndarray, t: float, dt: float) -> np.ndarray:
    k1 = h.phi(x, t)
    k2 = h.phi(x + dt * k1, t + dt)
    return x + 0.5 * dt * (k1 + k2)


def linear_first_step(h: FixedPointHomotopy, x: np.ndarray, t: float, dt: float) -> np.ndarray:
    """Linearised static homotopy at t1 = dt (eqs. 26-27): G(x0) + Gx(x0) dx = 0."""
    t1 = t + dt
    return x - solve(h.Gx(x, t1), h.G(x, t1), h.counter)


SCHEMES: dict[str, Step] = {
    "FE": forward_euler,
    "BE": backward_euler,
    "RK2": runge_kutta2,
    "LIN": linear_first_step,
}
