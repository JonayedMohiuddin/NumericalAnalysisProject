import time

import numpy as np

from ..linalg import LUCounter, factorize
from ..powerflow.network import make_b_xb
from ..results import SolveResult


def fast_decoupled_xb(pf, x0, tol=1e-8, max_it=100, record_states=False):
    """Fast decoupled power flow, XB version (eqs. 5-8), following MATPOWER's fdpf.

    B' and B'' are factorised once. Each iteration is a P-theta update
    followed by a Q-V update, and max_it limits the number of these cycles.

    The returned `iterations` is the number of P updates plus Q updates. This
    is not MATPOWER's return value: fdpf returns one count per P/Q cycle and
    only prints the P and Q counts separately. We use the sum because it is
    what the paper's Table 7 reports (e.g. 15 P + 14 Q updates = 29).
    """
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
