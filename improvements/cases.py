"""A wider test bed (the paper's 9 cases plus 25 more) and the reference solution of each case."""

import numpy as np

from dynhomotopy import datasets
from dynhomotopy.powerflow import PowerFlowProblem, load_case
from dynhomotopy.solvers import newton_raphson

from .multiplier import newton_raphson_om

EXTRA_ILL = ["case6024", "case6243", "case6748", "case7092", "case9961", "case10595", "case12110"]
EXTRA_LIMIT = [
    "case14limit", "case_ieee30limit", "case57limit", "case89pegaselimit", "case118limit",
    "case300limit", "case1354pegaselimit", "case2383wplimit", "case2736splimit",
    "case2737soplimit", "case2746woplimit", "case2746wplimit", "case2869pegaselimit",
    "case3012wplimit", "case3120splimit", "case3375wplimit", "case9241pegaselimit",
    "case13659pegaselimit",
]
EXTRA = EXTRA_ILL + EXTRA_LIMIT
ALL = datasets.ALL + EXTRA

_loaded = {}
_references = {}


def load(name):
    """A test system by name, or any MATPOWER case given as a path to a .mat file."""
    if name.endswith(".mat"):
        if name not in _loaded:
            _loaded[name] = PowerFlowProblem(load_case(name))
        return _loaded[name]
    if name in datasets.SYSTEMS:
        return datasets.load_problem(name)
    if name not in _loaded:
        record = datasets.ILL_CONDITIONED_RECORD if name in EXTRA_ILL else datasets.LIMIT_CASES_RECORD
        _loaded[name] = PowerFlowProblem(load_case(datasets.download(name, record)))
    return _loaded[name]


def reference_voltage(pf):
    """Bus voltages of the operating point reached from the case file's own initial guess.

    Power flow equations have more than one solution, and a small mismatch
    alone does not tell them apart. We take this solution as the intended
    operating point. It is found by NR (up to 50 iterations), with the
    optimal multiplier as a fallback for the stressed cases.
    """
    key = (pf.name, pf.case.base_mva)
    if key not in _references:
        res = newton_raphson(pf, pf.case_start(), max_it=50)
        if not res.converged:
            res = newton_raphson_om(pf, pf.case_start(), max_it=50)
        if not res.converged:
            raise RuntimeError(f"no reference solution for {pf.name}")
        _references[key] = pf.voltage(res.x)
    return _references[key]


def on_reference(pf, x, tol=1e-4):
    """True if state x is the reference operating point (largest voltage difference below tol pu)."""
    return float(np.max(np.abs(pf.voltage(x) - reference_voltage(pf)))) < tol
