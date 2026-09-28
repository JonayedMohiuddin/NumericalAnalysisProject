"""Our improvements to the paper's method, built on top of dynhomotopy without changing it.

Changes to the solver, each behind a switch of Options:
solve.py        solve() and Options
multiplier.py   NR with Iwamoto's optimal multiplier
steps.py        step rules (BE-chord, FE, RK2, RK4) and the Newton corrector
homotopies.py   the scaled and the Newton homotopy
path.py         fixed, adaptive and Richardson-controlled time points
cases.py        the 34-case test bed and the reference solutions

Studies of the method:
spectrum.py     power method and shifted inverse iteration for the Jacobian
feasibility.py  outcome of the paper's method on a (K, dt0) grid
tuning.py       golden-section search for delta and dt0
gauss.py        Gauss elimination and LU with partial pivoting, written from scratch
scratch.py      BE path, NR and FDXB using that LU

    from improvements import Options, solve
    result = solve(pf, pf.flat_start(), K=1e-4, times=[0, 0.005, 1], options=Options(multiplier=True))
"""

from .solve import Options, solve
