"""The methods compared in Tables 7 and 8."""

from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.solvers import gsh_nr, newton_raphson

from .common import K_PAPER

BE_PATH = pathway.explicit(0, 0.005, 1)
RK2_PATH = pathway.explicit(0, 0.005, 0.01, 0.02, 0.12, 0.42, 0.72, 1)

# (dh1, delta) for GSH-NR, taken from Table 7
GSH_PARAMS = {
    "case18482": (0.5, 1.0),
    "case27318": (0.25, 0.125),
    "case36964": (0.25, 1.0),
    "case54636": (0.01, 0.01),
    "case109272": (0.01, 0.01),
    "case69limit": (1.0, 1.0),
    "case141limit": (1.0, 1.0),
    "case_ACTIVSg500limit": (1.0, 1.0),
    "case_ACTIVSg2000limit": (1.0, 1.0),
}


def methods(name):
    """Map each method name to a function pf -> (converged, iterations)."""
    dh, delta = GSH_PARAMS[name]

    def outcome(res):
        if res is None:
            return False, 0
        return res.converged, res.iterations

    def hybrid(method, times, refiner):
        return lambda pf: outcome(solve_hybrid(pf, pf.flat_start(), K_PAPER, times, method, refiner).refine)

    return {
        "NR-MAT": lambda pf: outcome(newton_raphson(pf, pf.case_start())),
        "NR-flat": lambda pf: outcome(newton_raphson(pf, pf.flat_start())),
        "GSH-NR": lambda pf: outcome(gsh_nr(pf, pf.flat_start(), dh, delta)),
        "BE(NR)": hybrid("BE", BE_PATH, "NR"),
        "RK2(NR)": hybrid("RK2", RK2_PATH, "NR"),
        "BE(FDXB)": hybrid("BE", BE_PATH, "FDXB"),
        "RK2(FDXB)": hybrid("RK2", RK2_PATH, "FDXB"),
    }
