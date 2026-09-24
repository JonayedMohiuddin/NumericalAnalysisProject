"""Section 4.2.1: the first step to t1 = 0.005 with K = 1e-4.

An explicit first step moves x by -(dt0 / K) g(x0) = -50 g(x0), which blows up
when ||g(x0)|| is large. BE (or the linear approximation (27)) instead adds
delta = K / dt0 = 0.02 to the diagonal of the Jacobian.
"""

import numpy as np

from dynhomotopy.datasets import ALL, load_problem
from dynhomotopy.homotopy import integrate

from .common import DT0_PAPER, K_PAPER, label, write_table

RUNS = [(K_PAPER, "FE"), (K_PAPER, "RK2"), (K_PAPER, "BE"), (K_PAPER, "LIN"), (1.0, "FE")]


def first_step(pf, K, scheme):
    traj = integrate(pf, pf.flat_start(), K, [0.0, DT0_PAPER], scheme, first_step=None)
    if len(traj.norms) < 2:
        return np.inf, np.inf
    return traj.norms[-1], np.max(np.abs(traj.states[-1] - traj.states[0]))


def main(cases=ALL):
    rows = []
    for name in cases:
        pf = load_problem(name)
        row = [label(name), pf.norm(pf.flat_start())]
        for K, scheme in RUNS:
            row += first_step(pf, K, scheme)
        rows.append(row)

    header = ["case", "||g(x0)||"]
    for K, scheme in RUNS:
        header += [f"{scheme} K={K:g}: ||g(x(t1))||", f"{scheme} K={K:g}: max|dx|"]
    write_table("sec421_first_step", header, rows,
                f"Sec. 4.2.1: first step to t1 = {DT0_PAPER} (delta = K/dt0 = {K_PAPER / DT0_PAPER:g})")


if __name__ == "__main__":
    main()
