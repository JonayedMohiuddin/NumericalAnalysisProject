"""Figure 5: effect of the time step dt on bus 6 of case109272 (BE then NR, K = 1e-4)."""

from dynhomotopy.datasets import load_problem
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid

from .common import DT0_PAPER, K_PAPER, savefig, write_table
from .exp07_fig4 import CASE, bus_voltage, history, plot_split
from .plots import COLORS, new_figure

BUS = 6
STEPS = [1.0, 0.5, 0.25, 0.125, 0.1]


def main():
    pf = load_problem(CASE)
    fig, (ax_v, ax_a) = new_figure(1, 2, width=11, height=4.2)
    rows = []
    for color, dt in zip(COLORS, STEPS):
        times = pathway.constant_step(DT0_PAPER, dt)
        res = solve_hybrid(pf, pf.flat_start(), K_PAPER, times, "BE", "NR", record_states=True)
        states, norms, n_hom = history(res)
        vm, va = bus_voltage(pf, states, BUS)
        plot_split(ax_v, vm, n_hom, color, f"dt = {dt:g}")
        plot_split(ax_a, va, n_hom, color, f"dt = {dt:g}")
        rows.append([dt, n_hom, res.trajectory.factorizations, f"{vm[n_hom]:.4f}", f"{va[n_hom]:.2f}",
                     res.refine.iterations, res.factorizations, f"{vm[-1]:.4f}", f"{va[-1]:.2f}"])

    ax_v.set(xlabel="Iterations", ylabel="V (pu)", title=f"Fig. 5(a) voltage magnitude, bus {BUS}")
    ax_a.set(xlabel="Iterations", ylabel="Angle (degree)", title=f"Fig. 5(b) voltage angle, bus {BUS}")
    ax_v.legend()
    savefig(fig, "fig5_timestep_sensitivity")
    write_table("fig5_timestep",
                ["dt", "homotopy steps", "homotopy LUs", "V at t=1", "angle at t=1 (deg)",
                 "NR iterations", "total LUs", "final V", "final angle (deg)"],
                rows, f"Fig. 5 data: bus {BUS} of {CASE}")


if __name__ == "__main__":
    main()
