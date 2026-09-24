import numpy as np
import scipy.sparse as sp

from ..linalg import LUCounter, solve


class FixedPointHomotopy:
    """Homotopy of Section 3.1 with the fixed point function g0(x) = K (x - x0).

    G(x, t)  = t g(x) + (1 - t) K (x - x0)      eq. (12), (22)
    Gx(x, t) = t J(x) + (1 - t) K I             eq. (23)
    Gt(x)    = g(x) - K (x - x0)                eq. (24)

    Along the path dx/dt = -Gx^-1 Gt (eq. 15).

    `scale` replaces K I by K diag(scale), so each equation can get its own shift.
    """

    def __init__(self, problem, x0, K, counter=None, scale=None):
        self.problem = problem
        self.x0 = x0.astype(float)
        self.K = K
        self.counter = counter or LUCounter()
        self.scale = np.ones(problem.n) if scale is None else scale
        self.shift = sp.diags(self.scale, format="csc")

    def G(self, x, t):
        return t * self.problem.g(x) + (1 - t) * self.K * self.scale * (x - self.x0)

    def Gx(self, x, t):
        return (t * self.problem.jacobian(x) + (1 - t) * self.K * self.shift).tocsc()

    def Gt(self, x):
        return self.problem.g(x) - self.K * self.scale * (x - self.x0)

    def dxdt(self, x, t):
        return -solve(self.Gx(x, t), self.Gt(x), self.counter)


class NewtonHomotopy(FixedPointHomotopy):
    """Newton homotopy G(x, t) = g(x) - (1 - t) g(x0), a standard alternative to the FPV function.

    Here Gx = J(x) has no regularising shift, so it is included only for comparison.
    """

    def __init__(self, problem, x0, K=None, counter=None):
        super().__init__(problem, x0, 0.0, counter)
        self.g_x0 = problem.g(self.x0)

    def G(self, x, t):
        return self.problem.g(x) - (1 - t) * self.g_x0

    def Gx(self, x, t):
        return self.problem.jacobian(x)

    def Gt(self, x):
        return self.g_x0


def jacobian_diagonal_scale(problem, x0):
    """|diag J(x0)| normalised to mean 1, so K keeps the same average size."""
    d = np.abs(problem.jacobian(x0).diagonal())
    d[d == 0] = d[d > 0].min() if np.any(d > 0) else 1.0
    return d / d.mean()
