import time

import numpy as np

from ..linalg import LUCounter, SingularMatrixError, solve
from ..results import SolveResult


def optimal_multiplier(a, c):
    """Step length from Iwamoto and Tamura [10].

    With a = g(x) and c = g(x + dx) for the Newton step dx, the mismatch along
    the step is modelled as g(x + mu dx) = (1 - mu) a + mu^2 c. The returned mu
    in (0, 2] minimises the squared norm of this model.
    """
    aa, ac, cc = a @ a, a @ c, c @ c
    roots = np.roots([2 * cc, -3 * ac, aa + 2 * ac, -aa])
    real = roots[np.abs(roots.imag) < 1e-9].real
    candidates = [mu for mu in real if 0 < mu <= 2] or [1.0]
    return min(candidates, key=lambda mu: np.sum(((1 - mu) * a + mu**2 * c) ** 2))


def newton_raphson(problem, x0, tol=1e-8, max_it=10, record_states=False, multiplier=False):
    """Newton-Raphson, eq. (4). max_it = 10 is MATPOWER's default.

    With multiplier=True every step is scaled by Iwamoto's optimal multiplier.
    """
    start = time.perf_counter()
    counter = LUCounter()
    x = np.array(x0, dtype=float)
    f = problem.g(x)
    norms = [np.linalg.norm(f, np.inf)]
    states = [x] if record_states else None
    converged = norms[0] < tol
    it = 0

    while not converged and it < max_it:
        it += 1
        try:
            dx = -solve(problem.jacobian(x), f, counter)
        except SingularMatrixError:
            break
        f_new = problem.g(x + dx)
        if multiplier and np.all(np.isfinite(f_new)):
            mu = optimal_multiplier(f, f_new)
            if mu != 1.0:
                dx *= mu
                f_new = problem.g(x + dx)
        x = x + dx
        f = f_new
        norms.append(np.linalg.norm(f, np.inf))
        if record_states:
            states.append(x)
        if not np.isfinite(norms[-1]):
            break
        converged = norms[-1] < tol

    return SolveResult("NR", converged, it, x, norms, counter.factorizations,
                       time.perf_counter() - start, states)
