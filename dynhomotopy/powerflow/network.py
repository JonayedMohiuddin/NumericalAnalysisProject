import numpy as np
import scipy.sparse as sp

from . import idx


def make_ybus(base_mva, bus, branch):
    """Bus admittance matrix using the MATPOWER branch model (taps and phase shifters)."""
    nb, nl = len(bus), len(branch)
    status = branch[:, idx.BR_STATUS]
    ys = status / (branch[:, idx.BR_R] + 1j * branch[:, idx.BR_X])
    bc = status * branch[:, idx.BR_B]

    tap = np.ones(nl, dtype=complex)
    has_tap = branch[:, idx.TAP] != 0
    tap[has_tap] = branch[has_tap, idx.TAP]
    tap *= np.exp(1j * np.deg2rad(branch[:, idx.SHIFT]))

    ytt = ys + 1j * bc / 2
    yff = ytt / (tap * np.conj(tap))
    yft = -ys / np.conj(tap)
    ytf = -ys / tap
    ysh = (bus[:, idx.GS] + 1j * bus[:, idx.BS]) / base_mva

    f = branch[:, idx.F_BUS].astype(int)
    t = branch[:, idx.T_BUS].astype(int)
    diag = np.arange(nb)
    rows = np.concatenate([f, f, t, t, diag])
    cols = np.concatenate([f, t, f, t, diag])
    vals = np.concatenate([yff, yft, ytf, ytt, ysh])
    return sp.csr_matrix((vals, (rows, cols)), shape=(nb, nb))


def make_sbus(case):
    """Specified net power injection at each bus in p.u. (generation minus load)."""
    gen, bus = case.gen, case.bus
    sg = np.zeros(case.nb, dtype=complex)
    np.add.at(sg, gen[:, idx.GEN_BUS].astype(int), gen[:, idx.PG] + 1j * gen[:, idx.QG])
    sd = bus[:, idx.PD] + 1j * bus[:, idx.QD]
    return (sg - sd) / case.base_mva


def make_b_xb(case):
    """B' and B'' for the XB fast decoupled method (MATPOWER makeB with alg = 2).

    B' ignores shunts, line charging, taps and resistance.
    B'' is the full network without phase shifters.
    """
    bus_p = case.bus.copy()
    bus_p[:, idx.BS] = 0
    br_p = case.branch.copy()
    br_p[:, idx.BR_B] = 0
    br_p[:, idx.TAP] = 1
    br_p[:, idx.BR_R] = 0
    b_p = -make_ybus(case.base_mva, bus_p, br_p).imag

    br_pp = case.branch.copy()
    br_pp[:, idx.SHIFT] = 0
    b_pp = -make_ybus(case.base_mva, case.bus, br_pp).imag
    return sp.csr_matrix(b_p), sp.csr_matrix(b_pp)
