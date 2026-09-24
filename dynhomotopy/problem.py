from abc import ABC, abstractmethod

import numpy as np
import scipy.sparse as sp


class NonlinearProblem(ABC):
    """A square system g(x) = 0. All solvers and integrators work on this interface."""

    name = "problem"

    @property
    @abstractmethod
    def n(self):
        pass

    @abstractmethod
    def g(self, x):
        pass

    @abstractmethod
    def jacobian(self, x):
        pass

    def norm(self, x):
        return float(np.linalg.norm(self.g(x), np.inf))


class TutorialProblem(NonlinearProblem):
    """The 2x2 example of Section 3.3.

    g1 = x1^2 + x2^2 - 2 x1 x2 - 1
    g2 = x1 + x2 - 2

    The Jacobian is singular at the initial guess x0 = [1, 1].
    """

    name = "tutorial"
    x0 = np.array([1.0, 1.0])

    @property
    def n(self):
        return 2

    def g(self, x):
        x1, x2 = x
        return np.array([x1**2 + x2**2 - 2 * x1 * x2 - 1, x1 + x2 - 2])

    def jacobian(self, x):
        x1, x2 = x
        return sp.csc_matrix([[2 * x1 - 2 * x2, 2 * x2 - 2 * x1], [1.0, 1.0]])
