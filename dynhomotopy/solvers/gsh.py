"""GSH-NR: static "guided solution homotopy" solved point by point with NR.

Comparison method of Tables 7 and 8, from ref. [12] (Freitas & Silva, 2022).
This paper only summarises it (Section 2.3), so it is reconstructed here
from that summary:

* h1 inserts fictitious elements on the Ybus diagonal. At h = 0 they absorb
  exactly the power the flat start fails to balance, so the flat start is the
  solution of the "easy" network and no mismatch drives power between buses.
  Here they are the shunts y_f = conj(S_sp - S_calc(x0)) / |V0|^2, scaled by
  (1 - h), with only the real part kept at PV buses and none at the slack.
* h2 reduces the impedances of the branches connected to the slack bus by the
  factor delta at h = 0 and restores them linearly: z(h) = z (delta + (1 - delta) h).

Both parameters follow the same schedule h = dh, 2 dh, ..., 1. Each point is
solved by NR from the previous one, and the iterations are summed (iter_NR in
Table 7). With dh = 1 and delta = 1 the method reduces to NR from a flat start,
as in the limit-case rows of Table 7.
"""

from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp

from ..powerflow import idx
from ..powerflow.model import PowerFlowProblem
from ..powerflow.network import make_ybus
from ..results import SolveResult
from .newton import newton_raphson


def gsh_nr(pf: PowerFlowProblem, x0: np.ndarray, dh: float, delta: float,
           tol: float = 1e-8, max_it: int = 10) -> SolveResult:
    start = time.perf_counter()
    case = pf.case
    f = case.branch[:, idx.F_BUS].astype(int)
    t = case.branch[:, idx.T_BUS].astype(int)
    near_slack = np.isin(f, pf.ref) | np.isin(t, pf.ref)

    def network(h: float) -> sp.csr_matrix:
        br = case.branch.copy()
        scale = delta + (1.0 - delta) * h
        br[near_slack, idx.BR_R] *= scale
        br[near_slack, idx.BR_X] *= scale
        return make_ybus(case.base_mva, case.bus, br)

    v0 = pf.voltage(x0)
    ds = pf.sbus - v0 * np.conj(network(0.0) @ v0)
    ds[pf.pv] = ds[pf.pv].real
    ds[pf.ref] = 0.0
    y_fict = np.conj(ds) / np.abs(v0) ** 2

    n_points = int(round(1.0 / dh))
    x = np.array(x0, dtype=float)
    norms, iterations, lus, converged = [pf.norm(x)], 0, 0, True
    for k in range(1, n_points + 1):
        h = min(k * dh, 1.0)
        ybus = network(h) + sp.diags((1.0 - h) * y_fict)
        sub = PowerFlowProblem(case, ybus=ybus, sbus=pf.sbus)
        with np.errstate(all="ignore"):
            res = newton_raphson(sub, x, tol, max_it)
        iterations += res.iterations
        lus += res.factorizations
        x = res.x
        norms.append(pf.norm(x))
        if not res.converged:
            converged = False
            break
    return SolveResult("GSH-NR", converged, iterations, x, norms, lus, time.perf_counter() - start)
