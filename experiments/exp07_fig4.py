"""Figure 4: voltages of buses 2-6 of case109272 with only two homotopy points.

BE with dt0 = 0.005 and one large step to t = 1, then NR. Points 0-2 are the
homotopy (dotted) and the rest are NR iterations (solid). The paper puts the
mismatch on a second y axis; here it gets its own panel.
"""

import numpy as np

from dynhomotopy.datasets import load_problem
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid

from .common import DT0_PAPER, K_PAPER, savefig, write_table
from .plots import BACKGROUND, COLORS, new_figure

CASE = "case109272"
BUSES = [2, 3, 4, 5, 6]


def history(res):
    """States and norms over the homotopy points followed by the NR iterations."""
    states = res.trajectory.states + res.refine.states[1:]
    norms = res.trajectory.norms + res.refine.norms[1:]
    return np.array(states), np.array(norms), len(res.trajectory.states) - 1


def bus_voltage(pf, states, bus):
    b = pf.case.internal_bus(bus)
    return states[:, pf.magnitude_index(b)], np.rad2deg(states[:, pf.angle_index(b)])


def plot_split(ax, values, last_homotopy_point, color, label):
    k = np.arange(len(values))
    split = last_homotopy_point + 1
    ax.plot(k[:split], values[:split], ":", color=color)
    ax.plot(k[split - 1:], values[split - 1:], "-", color=color, label=label)


def main():
    pf = load_problem(CASE)
    times = pathway.constant_step(DT0_PAPER, 1.0)
    res = solve_hybrid(pf, pf.flat_start(), K_PAPER, times, "BE", "NR", record_states=True)
    states, norms, n_hom = history(res)

    fig, (ax_v, ax_a, ax_g) = new_figure(3, 1, width=6.4, height=8.4, sharex=True)
    rows = []
    for color, bus in zip(COLORS, BUSES):
        vm, va = bus_voltage(pf, states, bus)
        plot_split(ax_v, vm, n_hom, color, f"Bus-{bus}")
        plot_split(ax_a, va, n_hom, color, f"Bus-{bus}")
        rows.append([f"Bus-{bus}"] + [f"{v:.4f} pu / {a:.2f} deg" for v, a in zip(vm, va)])
    ax_g.semilogy(norms, "o-", color=COLORS[7], markeredgecolor=BACKGROUND)
    ax_g.axhline(1e-8, color="#52514e", lw=0.8, ls=":")
    ax_g.text(0.1, 2e-8, "tolerance 1e-8", color="#52514e", fontsize=8)

    ax_v.set(ylabel="V (pu)", title=f"Fig. 4 {CASE}: BE (dotted) then NR (solid)")
    ax_a.set(ylabel="Angle (degree)")
    ax_g.set(ylabel=r"$\|g(x)\|_\infty$", xlabel="Iterations")
    ax_v.legend(ncol=5, fontsize=8, loc="upper right")
    savefig(fig, "fig4_case109272_states")

    rows.append(["||g||_inf"] + [f"{n:.2e}" for n in norms])
    write_table("fig4_states", ["bus"] + [f"k={i}" for i in range(len(norms))], rows,
                f"Fig. 4 data: buses 2-6 of {CASE} (NR iterations: {res.refine.iterations})")


if __name__ == "__main__":
    main()
