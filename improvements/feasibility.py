"""Where in the (K, dt0) plane the paper's method reaches the reference solution."""

import numpy as np

from .solve import solve
from .cases import on_reference

FAILED, OTHER_ROOT, SOLVED = 0, 1, 2
DT0_VALUES = [0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2]
K_VALUES = list(np.logspace(-6, -1, 11))


def outcome(pf, K, dt0):
    res = solve(pf, pf.flat_start(), K, [0.0, dt0, 1.0])
    if not res.converged:
        return FAILED
    return SOLVED if on_reference(pf, res.refine.x) else OTHER_ROOT


def feasibility_grid(pf, K_values=K_VALUES, dt0_values=DT0_VALUES):
    """Outcome of the paper's method for every (dt0, K); rows are dt0, columns are K."""
    return np.array([[outcome(pf, K, dt0) for K in K_values] for dt0 in dt0_values])
