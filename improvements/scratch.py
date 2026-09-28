"""The paper's method (BE path, then NR or FDXB) with every linear system solved by our own LU.

These mirror dynhomotopy's solvers step by step, so their results can be
compared one to one with the SciPy (SuperLU) versions.
"""

import time

import numpy as np

from dynhomotopy.powerflow.network import make_b_xb
from dynhomotopy.results import SolveResult

from .gauss import factorize


def backward_euler_path(problem, x0, K, times):
    """The paper's BE path (eq. 19, one fixed point iteration per step). Returns (x(1), LUs)."""
    x0 = np.asarray(x0, dtype=float)
    eye = np.eye(problem.n)
    x = x0
    for t, t_next in zip(times, times[1:]):
        gx = t_next * problem.jacobian(x).toarray() + (1 - t_next) * K * eye
        gt = problem.g(x) - K * (x - x0)
        x = x - (t_next - t) * factorize(gx)(gt)
    return x, len(times) - 1


def newton_raphson(problem, x0, tol=1e-8, max_it=10):
    start = time.perf_counter()
    x = np.array(x0, dtype=float)
    f = problem.g(x)
    norms = [np.linalg.norm(f, np.inf)]
    it = 0
    while norms[-1] >= tol and it < max_it and np.isfinite(norms[-1]):
        it += 1
        x = x - factorize(problem.jacobian(x))(f)
        f = problem.g(x)
        norms.append(np.linalg.norm(f, np.inf))
    return SolveResult("NR (own LU)", norms[-1] < tol, it, x, norms, it, time.perf_counter() - start)


def fast_decoupled_xb(pf, x0, tol=1e-8, max_it=100):
    """FDXB with B' and B'' factorised once by our LU. Counts P plus Q updates, like dynhomotopy."""
    start = time.perf_counter()
    b_p, b_pp = make_b_xb(pf.case)
    solve_p = factorize(b_p[pf.pvpq][:, pf.pvpq])
    solve_q = factorize(b_pp[pf.pq][:, pf.pq])
    n_angles = pf.npv + pf.npq

    def mismatch(x):
        v = pf.voltage(x)
        mis = pf.complex_mismatch(v) / np.abs(v)
        return mis[pf.pvpq].real, mis[pf.pq].imag

    def small(p, q):
        return np.linalg.norm(p, np.inf) < tol and np.linalg.norm(q, np.inf) < tol

    x = np.array(x0, dtype=float)
    norms = [pf.norm(x)]
    p, q = mismatch(x)
    converged = small(p, q)
    cycles = updates = 0
    while not converged and cycles < max_it and np.isfinite(norms[-1]):
        cycles += 1
        for part, rhs_solve in ((slice(0, n_angles), solve_p), (slice(n_angles, None), solve_q)):
            x[part] -= rhs_solve(p if part.start == 0 else q)
            updates += 1
            p, q = mismatch(x)
            norms.append(pf.norm(x))
            converged = small(p, q)
            if converged or not np.isfinite(norms[-1]):
                break
    return SolveResult("FDXB (own LU)", converged, updates, x, norms, 2, time.perf_counter() - start)
