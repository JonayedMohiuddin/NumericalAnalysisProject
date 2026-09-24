"""Abstract nonlinear system g(x) = 0 shared by every solver in the package.

Both the tutorial example (Section 3.3) and the power-flow problem (Section 2)
implement this interface, so Newton-Raphson and the dynamic-homotopy
integrators are written once and applied to either.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
import scipy.sparse as sp


class NonlinearProblem(ABC):
    """A square nonlinear system g: R^n -> R^n with a sparse Jacobian."""

    name: str = "problem"

    @property
    @abstractmethod
    def n(self) -> int:
        """Number of unknowns / equations."""

    @abstractmethod
    def g(self, x: np.ndarray) -> np.ndarray:
        """Mismatch vector g(x)."""

    @abstractmethod
    def jacobian(self, x: np.ndarray) -> sp.csc_matrix:
        """Jacobian J(x) = dg/dx (eq. 3)."""

    def norm(self, x: np.ndarray) -> float:
        """Infinity norm of the mismatch, ||g(x)||_inf, used for every table/figure."""
        return float(np.linalg.norm(self.g(x), np.inf))


class TutorialProblem(NonlinearProblem):
    """Generic 2x2 tutorial system of Section 3.3 (taken from ref. [29]).

        g1(x1, x2) = x1^2 + x2^2 - 2 x1 x2 - 1 = 0
        g2(x1, x2) = x1 + x2 - 2               = 0

    The root reached from x0 = [1, 1] is x* = [1.5, 0.5]; the Jacobian is
    singular at x0 itself.
    """

    name = "tutorial"
    x0 = np.array([1.0, 1.0])

    @property
    def n(self) -> int:
        return 2

    def g(self, x: np.ndarray) -> np.ndarray:
        x1, x2 = x
        return np.array([x1**2 + x2**2 - 2 * x1 * x2 - 1.0, x1 + x2 - 2.0])

    def jacobian(self, x: np.ndarray) -> sp.csc_matrix:
        x1, x2 = x
        return sp.csc_matrix(np.array([[2 * x1 - 2 * x2, 2 * x2 - 2 * x1], [1.0, 1.0]]))
