"""Newton-Raphson with the optimal multiplier of Iwamoto and Tamura (ref. [10] of the paper)."""

import time

import numpy as np

from dynhomotopy.linalg import LUCounter, SingularMatrixError, solve
from dynhomotopy.results import SolveResult


def optimal_multiplier(a, c):
    """Step length mu for the Newton step dx.

    With a = g(x) and c = g(x + dx), the mismatch along the step is modelled
    as g(x + mu dx) = (1 - mu) a + mu^2 c. We return the mu in (0, 2] that
    minimises the squared norm of this model. Setting its derivative to zero
    gives a cubic in mu.
    """
    aa, ac, cc = a @ a, a @ c, c @ c
    roots = np.roots([2 * cc, -3 * ac, aa + 2 * ac, -aa])
    real = roots[np.abs(roots.imag) < 1e-9].real
    candidates = [mu for mu in real if 0 < mu <= 2] or [1.0]
    return min(candidates, key=lambda mu: np.sum(((1 - mu) * a + mu**2 * c) ** 2))


def newton_raphson_om(problem, x0, tol=1e-8, max_it=10, record_states=False):
    """Same as dynhomotopy.solvers.newton_raphson, but every step is scaled by mu.

    Close to the solution mu is 1, so the fast convergence of NR is kept.
    The multiplier needs g(x + dx), which NR computes anyway, so a step only
    costs one extra mismatch evaluation when mu is not 1, and no extra LU.
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
        f_full = problem.g(x + dx)
        mu = optimal_multiplier(f, f_full) if np.all(np.isfinite(f_full)) else 1.0
        x = x + mu * dx
        f = f_full if mu == 1.0 else problem.g(x)
        norms.append(np.linalg.norm(f, np.inf))
        if record_states:
            states.append(x)
        if not np.isfinite(norms[-1]):
            break
        converged = norms[-1] < tol

    return SolveResult("NR-OM", converged, it, x, norms, counter.factorizations,
                       time.perf_counter() - start, states)
