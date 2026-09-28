import time

import numpy as np
import scipy.sparse as sp
from pypower.idx_brch import BR_R, BR_X, F_BUS, T_BUS

from ..powerflow.model import PowerFlowProblem
from ..powerflow.network import make_ybus
from ..results import SolveResult
from .newton import newton_raphson


def gsh_nr(pf, x0, dh, delta, tol=1e-8, max_it=10):
    """Static homotopy GSH-NR from ref. [12], used for comparison in Tables 7 and 8.

    The paper only describes this method in words, so this is our version of it:

    - Fictitious shunts are added to the diagonal of Ybus so that the flat
      start solves the network exactly at h = 0. They are scaled by (1 - h).
    - The impedance of branches connected to the slack bus is multiplied by
      delta at h = 0 and brought back to its real value at h = 1.

    h goes from dh to 1 in steps of dh. Each point is solved with NR starting
    from the previous one, and the NR iterations are added up.
    With dh = 1 and delta = 1 this is just NR from a flat start.
    """
    start = time.perf_counter()
    case = pf.case
    f = case.branch[:, F_BUS].astype(int)
    t = case.branch[:, T_BUS].astype(int)
    near_slack = np.isin(f, pf.ref) | np.isin(t, pf.ref)

    def network(h):
        branch = case.branch.copy()
        branch[near_slack, BR_R] *= delta + (1 - delta) * h
        branch[near_slack, BR_X] *= delta + (1 - delta) * h
        return make_ybus(case.base_mva, case.bus, branch)

    # shunts that absorb the flat start mismatch (only P at PV buses, nothing at the slack)
    v0 = pf.voltage(x0)
    ds = pf.sbus - v0 * np.conj(network(0) @ v0)
    ds[pf.pv] = ds[pf.pv].real
    ds[pf.ref] = 0
    y_shunt = np.conj(ds) / np.abs(v0) ** 2

    x = np.array(x0, dtype=float)
    norms = [pf.norm(x)]
    iterations = lus = 0
    converged = True
    for k in range(1, round(1 / dh) + 1):
        h = min(k * dh, 1.0)
        ybus = network(h) + sp.diags((1 - h) * y_shunt)
        with np.errstate(all="ignore"):
            res = newton_raphson(PowerFlowProblem(case, ybus, pf.sbus), x, tol, max_it)
        iterations += res.iterations
        lus += res.factorizations
        x = res.x
        norms.append(pf.norm(x))
        if not res.converged:
            converged = False
            break

    return SolveResult("GSH-NR", converged, iterations, x, norms, lus, time.perf_counter() - start)
