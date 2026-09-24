"""Classical Newton-Raphson (Section 2.1.1, eq. 4), a port of MATPOWER ``newtonpf``."""

from __future__ import annotations

import time

import numpy as np

from ..linalg import LUCounter, SingularMatrixError, solve
from ..problem import NonlinearProblem
from ..results import SolveResult


def newton_raphson(problem: NonlinearProblem, x0: np.ndarray, tol: float = 1e-8,
                   max_it: int = 10, record_states: bool = False) -> SolveResult:
    """x_{i+1} = x_i - J(x_i)^{-1} g(x_i) until ||g||_inf < tol.

    ``max_it`` = 10 is MATPOWER's default; a run that does not meet ``tol``
    within it (or produces a singular Jacobian / non-finite state) is a fail.
    """
    start = time.perf_counter()
    counter = LUCounter()
    x = np.array(x0, dtype=float)
    f = problem.g(x)
    norms = [float(np.linalg.norm(f, np.inf))]
    states = [x.copy()] if record_states else None
    converged = norms[0] < tol
    it = 0
    while not converged and it < max_it:
        it += 1
        try:
            dx = -solve(problem.jacobian(x), f, counter)
        except SingularMatrixError:
            break
        x = x + dx
        f = problem.g(x)
        norms.append(float(np.linalg.norm(f, np.inf)))
        if states is not None:
            states.append(x.copy())
        if not np.isfinite(norms[-1]):
            break
        converged = norms[-1] < tol
    return SolveResult("NR", converged, it, x, norms, counter.factorizations,
                       time.perf_counter() - start, states)
