"""Solver configurations of Section 4.7 shared by Table 7 (iterations) and Table 8 (CPU time)."""

from __future__ import annotations

from typing import Callable

from common import K_PAPER
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.powerflow import PowerFlowProblem
from dynhomotopy.solvers import gsh_nr, newton_raphson

# BE: only two pathway points; RK2 needs small early steps (seven points)
BE_PATH = pathway.explicit(0, 0.005, 1.0)
RK2_PATH = pathway.explicit(0, 0.005, 0.01, 0.02, 0.12, 0.42, 0.72, 1.0)

# GSH-NR parameters (dh1, delta) per case, as listed in Table 7
GSH_PARAMS = {
    "case18482": (0.5, 1.0), "case27318": (0.25, 0.125), "case36964": (0.25, 1.0),
    "case54636": (0.01, 0.01), "case109272": (0.01, 0.01),
    "case69limit": (1.0, 1.0), "case141limit": (1.0, 1.0),
    "case_ACTIVSg500limit": (1.0, 1.0), "case_ACTIVSg2000limit": (1.0, 1.0),
}


def methods(name: str) -> dict[str, Callable[[PowerFlowProblem], tuple[bool, int]]]:
    """method label -> fn(pf) returning (converged, iterations of the final solver)."""
    dh, delta = GSH_PARAMS[name]

    def nr(start):
        return lambda pf: _it(newton_raphson(pf, getattr(pf, start)()))

    def hybrid(method, times, refiner):
        return lambda pf: _it(solve_hybrid(pf, pf.flat_start(), K_PAPER, times, method, refiner).refine)

    return {
        "NR-MAT": nr("case_start"),
        "NR-flat": nr("flat_start"),
        "GSH-NR": lambda pf: _it(gsh_nr(pf, pf.flat_start(), dh, delta)),
        "BE(NR)": hybrid("BE", BE_PATH, "NR"),
        "RK2(NR)": hybrid("RK2", RK2_PATH, "NR"),
        "BE(FDXB)": hybrid("BE", BE_PATH, "FDXB"),
        "RK2(FDXB)": hybrid("RK2", RK2_PATH, "FDXB"),
    }


def _it(res) -> tuple[bool, int]:
    return (res is not None and res.converged), (res.iterations if res is not None else 0)
