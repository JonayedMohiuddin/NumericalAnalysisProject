import scipy.sparse as sp
from pypower.makeB import makeB
from pypower.makeSbus import makeSbus
from pypower.makeYbus import makeYbus


def make_ybus(base_mva, bus, branch):
    """Bus admittance matrix (MATPOWER branch model, with taps and phase shifters)."""
    return sp.csr_matrix(makeYbus(base_mva, bus, branch)[0])


def make_sbus(case):
    """Specified net power injection at each bus in p.u. (generation minus load)."""
    return makeSbus(case.base_mva, case.bus, case.gen)


def make_b_xb(case):
    """B' and B'' for the XB fast decoupled method (makeB with alg = 2)."""
    b_p, b_pp = makeB(case.base_mva, case.bus, case.branch, 2)
    return sp.csr_matrix(b_p), sp.csr_matrix(b_pp)
