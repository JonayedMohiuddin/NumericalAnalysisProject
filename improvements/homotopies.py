"""Other choices for the easy function g0 (the paper's future work, Section 5).

Both have the same methods as dynhomotopy's FixedPointHomotopy, so the step
rules work with any of them. Neither helped in our tests.
"""

import numpy as np
import scipy.sparse as sp

from dynhomotopy.homotopy.fpv import FixedPointHomotopy


class ScaledHomotopy(FixedPointHomotopy):
    """g0(x) = K D (x - x0), with D = |diag J(x0)| scaled to mean 1.

    Each equation gets a shift in proportion to its own diagonal, in the
    spirit of the conditioning step of ref. [29].
    """

    def __init__(self, problem, x0, K, counter=None):
        super().__init__(problem, x0, K, counter)
        d = np.abs(problem.jacobian(self.x0).diagonal())
        d[d == 0] = d[d > 0].min()
        self.scale = d / d.mean()
        self.eye = sp.diags(self.scale, format="csc")  # Gx uses this in place of I

    def G(self, x, t):
        return t * self.problem.g(x) + (1 - t) * self.K * self.scale * (x - self.x0)

    def Gt(self, x):
        return self.problem.g(x) - self.K * self.scale * (x - self.x0)


class NewtonHomotopy(FixedPointHomotopy):
    """G(x, t) = g(x) - (1 - t) g(x0), the standard Newton homotopy (ref. [36]).

    Gx = J(x) has no K I shift, so the ill-conditioned Jacobian is not helped.
    K is ignored.
    """

    def __init__(self, problem, x0, K=0.0, counter=None):
        super().__init__(problem, x0, K, counter)
        self.g_x0 = problem.g(self.x0)

    def G(self, x, t):
        return self.problem.g(x) - (1 - t) * self.g_x0

    def Gx(self, x, t):
        return self.problem.jacobian(x)

    def Gt(self, x):
        return self.g_x0


HOMOTOPIES = {
    "fpv": FixedPointHomotopy,
    "scaled": ScaledHomotopy,
    "newton": NewtonHomotopy,
}
