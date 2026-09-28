import time

import numpy as np

from ..linalg import LUCounter, factorize
from ..powerflow.network import make_b_xb
from ..results import SolveResult


def fast_decoupled_xb(pf, x0, tol=1e-8, max_it=100, record_states=False):
    start = time.perf_counter()
    counter = LUCounter()
    b_p, b_pp = make_b_xb(pf.case)
    solve_p = factorize(b_p[pf.pvpq][:, pf.pvpq], counter)
    solve_q = factorize(b_pp[pf.pq][:, pf.pq], counter)
    n_angles = pf.npv + pf.npq

    def mismatch(x):
        v = pf.voltage(x)
        mis = pf.complex_mismatch(v) / np.abs(v)
        return mis[pf.pvpq].real, mis[pf.pq].imag

    def small(p, q):
        return np.linalg.norm(p, np.inf) < tol and np.linalg.norm(q, np.inf) < tol

    x = np.array(x0, dtype=float)
    norms = [pf.norm(x)]
    states = [x.copy()] if record_states else None
    p, q = mismatch(x)
    converged = small(p, q)
    p_its = half_steps = 0

    while not converged and p_its < max_it and np.isfinite(norms[-1]):
        p_its += 1
        for part in ("P", "Q"):
            if part == "P":
                x[:n_angles] -= solve_p(p)
            else:
                x[n_angles:] -= solve_q(q)
            half_steps += 1
            p, q = mismatch(x)
            norms.append(pf.norm(x))
            if record_states:
                states.append(x.copy())
            converged = small(p, q)
            if converged or not np.isfinite(norms[-1]):
                break

    return SolveResult("FDXB", converged, half_steps, x, norms, counter.factorizations,
                       time.perf_counter() - start, states)
