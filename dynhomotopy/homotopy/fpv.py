"""Fixed-point-vector (FPV) homotopy function of Section 3.1.

    G(x, t)  = t g(x) + (1 - t) g0(x),      g0(x) = K (x - x0)          (12), (22)
    Gx(x, t) = t J(x) + (1 - t) K I                                    (23)
    Gt(x, t) = g(x) - g0(x)                                            (24)

Differentiating G(x(t), t) = 0 along the pathway gives the ODE (15)

    Gx dx/dt = -Gt,     x(0) = x0,     0 <= t <= 1,

whose right-hand side  Phi(x, t) = -Gx^{-1} Gt  costs one LU factorisation.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from ..linalg import LUCounter, solve
from ..problem import NonlinearProblem


class FixedPointHomotopy:
    """G(x, t) built from a problem g(x) = 0, its initial estimate x0 and factor K."""

    def __init__(self, problem: NonlinearProblem, x0: np.ndarray, K: float,
                 counter: LUCounter | None = None):
        self.problem = problem
        self.x0 = np.array(x0, dtype=float)
        self.K = float(K)
        self.counter = counter if counter is not None else LUCounter()
        self._eye = sp.identity(problem.n, format="csc")

    def g0(self, x: np.ndarray) -> np.ndarray:
        return self.K * (x - self.x0)

    def G(self, x: np.ndarray, t: float) -> np.ndarray:
        return t * self.problem.g(x) + (1.0 - t) * self.g0(x)

    def Gx(self, x: np.ndarray, t: float) -> sp.csc_matrix:
        return (t * self.problem.jacobian(x) + (1.0 - t) * self.K * self._eye).tocsc()

    def Gt(self, x: np.ndarray) -> np.ndarray:
        return self.problem.g(x) - self.g0(x)

    def phi(self, x: np.ndarray, t: float) -> np.ndarray:
        """ODE right-hand side dx/dt = Phi(x, t) = -Gx(x, t)^{-1} Gt(x)."""
        return -solve(self.Gx(x, t), self.Gt(x), self.counter)
