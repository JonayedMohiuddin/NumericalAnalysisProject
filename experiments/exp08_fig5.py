"""Section 4.6 - Figure 5: sensitivity to the time step dt (> dt0).

BE with dt0 = 0.005, K = 1e-4 and constant dt in {0.995 (shown as 1), 0.5,
0.25, 0.125, 0.1}, then NR. Traces bus #6 of case109272.
"""

from __future__ import annotations

import numpy as np

from common import DT0_PAPER, K_PAPER, savefig, write_table
from dynhomotopy.datasets import load_problem
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid
from exp07_fig4 import CASE, draw_split, state_history, trace
from plots import SERIES, new_figure

BUS = 6
STEPS = [1.0, 0.5, 0.25, 0.125, 0.1]


def main():
    pf = load_problem(CASE)
    fig, (ax_v, ax_a) = new_figure(1, 2, width=11, height=4.2)
    rows = []
    for color, dt in zip(SERIES, STEPS):
        res = solve_hybrid(pf, pf.flat_start(), K_PAPER, pathway.constant_step(DT0_PAPER, dt), "BE", "NR",
                           record_states=True)
        states, norms, n_hom = state_history(res)
        vm, va = trace(pf, states, BUS)
        draw_split(ax_v, vm, n_hom, color, f"dt = {dt:g}")
        draw_split(ax_a, va, n_hom, color, f"dt = {dt:g}")
        rows.append([dt, n_hom, res.trajectory.factorizations, f"{vm[n_hom]:.4f}", f"{va[n_hom]:.2f}",
                     res.refine.iterations, res.factorizations, f"{vm[-1]:.4f}", f"{va[-1]:.2f}"])
    ax_v.set(xlabel="Iterations", ylabel="V (pu)", title=f"Fig. 5(a) voltage magnitude, bus {BUS}")
    ax_a.set(xlabel="Iterations", ylabel="Angle (degree)", title=f"Fig. 5(b) voltage angle, bus {BUS}")
    ax_v.legend()
    savefig(fig, "fig5_timestep_sensitivity")
    write_table("fig5_timestep", ["dt", "pathway steps", "homotopy LUs", "V at t=1", "angle at t=1 (deg)",
                                  "NR iterations", "total LUs", "V final", "angle final (deg)"], rows,
                f"Fig. 5 data: bus {BUS} of {CASE}, BE + NR")


if __name__ == "__main__":
    main()
