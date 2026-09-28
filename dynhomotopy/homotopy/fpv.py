import scipy.sparse as sp

from ..linalg import LUCounter, solve


class FixedPointHomotopy:
    """
    Homotopy with the fixed point function g0(x) = K (x - x0).
    G(x, t)  = t g(x) + (1 - t) K (x - x0)     
    Gx(x, t) = t J(x) + (1 - t) K I         
    Gt(x)    = g(x) - K (x - x0)              
    Along the path dx/dt = -Gx^-1 Gt
    """

    def __init__(self, problem, x0, K, counter=None):
        self.problem = problem
        self.x0 = x0.astype(float)
        self.K = K
        self.counter = counter or LUCounter()
        self.eye = sp.identity(problem.n, format="csc")

    def G(self, x, t):
        return t * self.problem.g(x) + (1 - t) * self.K * (x - self.x0)

    def Gx(self, x, t):
        return (t * self.problem.jacobian(x) + (1 - t) * self.K * self.eye).tocsc()

    def Gt(self, x):
        return self.problem.g(x) - self.K * (x - self.x0)

    def dxdt(self, x, t):
        return -solve(self.Gx(x, t), self.Gt(x), self.counter)
