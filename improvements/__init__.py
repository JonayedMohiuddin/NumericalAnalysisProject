"""Our improvements to the paper's method, built on top of dynhomotopy without changing it.

    from improvements import Options, solve
    result = solve(pf, pf.flat_start(), K=1e-4, times=[0, 0.005, 1], options=Options(multiplier=True))
"""

from .solve import Options, solve
