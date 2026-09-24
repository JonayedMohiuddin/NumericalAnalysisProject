"""Fast decoupled power flow, XB version (Section 2.1.1, eqs. 5-8).

Port of MATPOWER ``fdpf`` with ``makeB(alg=2)``. B' and B'' are constant, so
each is factorised exactly once, and the P-theta and Q-V half iterations
alternate:

    dtheta = -B'^{-1}  dP/V,        dV = -B''^{-1} dQ/V
"""

from __future__ import annotations

import time

import numpy as np

from ..linalg import LUCounter, factorize
from ..powerflow.model import PowerFlowProblem
from ..powerflow.network import make_b_xb
from ..results import SolveResult


def fast_decoupled_xb(pf: PowerFlowProblem, x0: np.ndarray, tol: float = 1e-8,
                      max_it: int = 100, record_states: bool = False) -> SolveResult:
    """FDXB from ``x0``. Convergence uses MATPOWER's criterion on dP/V and dQ/V.

    ``max_it`` bounds the number of P iterations (MATPOWER's ``i``). The
    reported ``iterations`` is the number of P- plus Q-half-iterations, which
    is what MATPOWER prints and what Table 7 of the paper tabulates.
    ``norms`` holds the true ||g(x)||_inf after every half iteration.
    """
    start = time.perf_counter()
    counter = LUCounter()
    bp, bpp = make_b_xb(pf.case)
    solve_p = factorize(bp[pf.pvpq][:, pf.pvpq], counter)
    solve_q = factorize(bpp[pf.pq][:, pf.pq], counter)

    x = np.array(x0, dtype=float)
    na = pf.npv + pf.npq

    def scaled_mismatch(x):
        v = pf.voltage(x)
        mis = pf.complex_mismatch(v) / np.abs(v)
        return mis[pf.pvpq].real, mis[pf.pq].imag

    def done(p, q):
        return np.linalg.norm(p, np.inf) < tol and np.linalg.norm(q, np.inf) < tol

    norms = [pf.norm(x)]
    states = [x.copy()] if record_states else None
    p, q = scaled_mismatch(x)
    converged = done(p, q)
    p_its = half_its = 0
    while not converged and p_its < max_it:
        p_its += 1
        for part, rhs_solve in ((slice(0, na), solve_p), (slice(na, None), solve_q)):
            half_its += 1
            x[part] -= rhs_solve(p if part.start == 0 else q)   # P-theta, then Q-V
            p, q = scaled_mismatch(x)
            norms.append(pf.norm(x))
            if states is not None:
                states.append(x.copy())
            converged = done(p, q)
            if converged or not np.isfinite(norms[-1]):
                break
        if not np.isfinite(norms[-1]):
            break
    return SolveResult("FDXB", converged, half_its, x, norms, counter.factorizations,
                       time.perf_counter() - start, states)
