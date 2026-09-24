"""Section 4.2.1 - the first time step t1 = dt0 = 0.005 with K = 1e-4.

Explicit first steps move by dx(0) = -(dt0/K) g(x0) = -50 g(x0), which blows up
when ||g(x0)|| >> 1. The implicit step (BE, or the linear approximation (27))
adds the shift delta = K/dt0 = 0.02 to the Jacobian diagonal instead (eq. 25).
The table lists ||g(x(t1))|| for every first-step scheme, plus FE with K = 1.
"""

from __future__ import annotations

import numpy as np

from common import DT0_PAPER, K_PAPER, label, write_table
from dynhomotopy.datasets import ALL, load_problem
from dynhomotopy.homotopy import integrate


def first_step_norm(pf, K, scheme):
    tr = integrate(pf, pf.flat_start(), K, [0.0, DT0_PAPER], scheme, first_step=None)
    return tr.norms[-1] if len(tr.norms) > 1 else np.inf, np.max(np.abs(tr.states[-1] - tr.states[0]))


def main(cases=ALL):
    rows = []
    for name in cases:
        pf = load_problem(name)
        row = [label(name), pf.norm(pf.flat_start())]
        for K, scheme in [(K_PAPER, "FE"), (K_PAPER, "RK2"), (K_PAPER, "BE"), (K_PAPER, "LIN"), (1.0, "FE")]:
            nrm, step = first_step_norm(pf, K, scheme)
            row += [nrm, step]
        rows.append(row)
        print(label(name), "done")
    hdr = ["case", "||g(x0)||"]
    for tag in ["FE K=1e-4", "RK2 K=1e-4", "BE K=1e-4", "LIN (27) K=1e-4", "FE K=1"]:
        hdr += [f"{tag}: ||g(x(t1))||", f"{tag}: max|dx(0)|"]
    write_table("sec421_first_step", hdr, rows,
                f"Sec. 4.2.1: states at t1 = {DT0_PAPER} (delta = K/dt0 = {K_PAPER / DT0_PAPER:g})")


if __name__ == "__main__":
    main()
