import time

import numpy as np

from ..linalg import LUCounter, SingularMatrixError, solve
from ..results import SolveResult


def newton_raphson(problem, x0, tol=1e-8, max_it=10, record_states=False):
    """Newton-Raphson, eq. (4). max_it = 10 is MATPOWER's default."""
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
            x = x - solve(problem.jacobian(x), f, counter)
        except SingularMatrixError:
            break
        f = problem.g(x)
        norms.append(np.linalg.norm(f, np.inf))
        if record_states:
            states.append(x)
        if not np.isfinite(norms[-1]):
            break
        converged = norms[-1] < tol

    return SolveResult("NR", converged, it, x, norms, counter.factorizations,
                       time.perf_counter() - start, states)
