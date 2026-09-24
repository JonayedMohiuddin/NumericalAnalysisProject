"""Section 4.5 - Figure 4: states and mismatch for case109272 with a high time step.

BE with dt0 = 0.005, K = 1e-4 and dt = 0.995 (pathway {0, 0.005, 1}), then NR.
The magnitude and angle of buses 2-6 (all PQ) are traced over homotopy points
k = 0..2 (dotted) and NR iterations 3.. (solid).
The paper draws the mismatch on a second y-axis. Here it gets its own panel
below the states, on a shared x-axis.
"""

from __future__ import annotations

import numpy as np

from common import DT0_PAPER, K_PAPER, savefig, write_table
from dynhomotopy.datasets import load_problem
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid
from plots import SERIES, new_figure

CASE = "case109272"
BUSES = [2, 3, 4, 5, 6]


def state_history(res):
    """x at every homotopy point, then after every NR iteration, plus the joint norm list."""
    states = res.trajectory.states + res.refine.states[1:]
    norms = res.trajectory.norms + res.refine.norms[1:]
    return np.array(states), np.array(norms), len(res.trajectory.states) - 1


def trace(pf, states, bus):
    b = pf.case.internal_bus(bus)
    return states[:, pf.magnitude_index(b)], np.rad2deg(states[:, pf.angle_index(b)])


def draw_split(ax, values, n_hom, color, label):
    """Dotted over the homotopy points, solid over the NR iterations."""
    k = np.arange(len(values))
    ax.plot(k[: n_hom + 1], values[: n_hom + 1], ":", color=color)
    ax.plot(k[n_hom:], values[n_hom:], "-", color=color, label=label)


def main():
    pf = load_problem(CASE)
    times = pathway.constant_step(DT0_PAPER, 1.0)          # [0, 0.005, 1.0]
    res = solve_hybrid(pf, pf.flat_start(), K_PAPER, times, "BE", "NR", record_states=True)
    states, norms, n_hom = state_history(res)

    fig, (ax_v, ax_a, ax_g) = new_figure(3, 1, width=6.4, height=8.4, sharex=True)
    rows = []
    for color, bus in zip(SERIES, BUSES):
        vm, va = trace(pf, states, bus)
        draw_split(ax_v, vm, n_hom, color, f"Bus-{bus}")
        draw_split(ax_a, va, n_hom, color, f"Bus-{bus}")
        rows.append([f"Bus-{bus}", *[f"{v:.4f} pu / {a:.2f} deg" for v, a in zip(vm, va)]])
    ax_g.semilogy(np.arange(len(norms)), norms, "o-", color=SERIES[7], markeredgecolor="#fcfcfb")
    ax_g.axhline(1e-8, color="#52514e", lw=0.8, ls=(0, (2, 3)))
    ax_g.text(0.1, 2e-8, "tolerance 1e-8", color="#52514e", fontsize=8)
    ax_v.set(ylabel="V (pu)", title=f"Fig. 4 {CASE}: BE (dotted, k = 0..{n_hom}) then NR (solid)")
    ax_a.set(ylabel="Angle (degree)")
    ax_g.set(ylabel=r"$\|g(x)\|_\infty$", xlabel="Iterations")
    ax_v.legend(ncol=5, fontsize=8, loc="upper right")
    savefig(fig, "fig4_case109272_states")
    write_table("fig4_states", ["bus"] + [f"k={i}" for i in range(len(norms))],
                rows + [["||g||_inf", *[f"{n:.2e}" for n in norms]]],
                f"Fig. 4 data: V / angle of buses 2-6 of {CASE} (NR iterations: {res.refine.iterations})")


if __name__ == "__main__":
    main()
