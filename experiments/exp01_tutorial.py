"""Figures 1 and 2: the 2x2 tutorial example of Section 3.3.

Fig. 1(a) runs NR alone from x0 = [1, 1 - eps]. Figs. 1(b)-(d) and 2 run the
homotopy from x0 = [1, 1] (dt0 = 0.01) and then NR. The tolerance is 1e-5.
Homotopy points are drawn at their time t in [0, 1] and NR iteration i at 1 + i.
"""

import numpy as np

from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.problem import TutorialProblem
from dynhomotopy.solvers import newton_raphson

from .common import pad, savefig, write_table
from .paper_values import FIG1_NR_ITERS, FIG1A_ITERS
from .plots import BACKGROUND, COLORS, METHOD_COLOR, METHOD_MARKER, new_figure

TOL = 1e-5
METHODS = ["FE", "BE", "RK2"]
HYBRID_RUNS = [("fig1b", 0.5, 0.05), ("fig1c", 0.5, 0.01), ("fig1d", 0.5, 0.005), ("fig2", 0.125, 0.05)]


def plot_hybrid(runs, title, filename):
    fig, ax = new_figure()
    for m, r in runs.items():
        traj = r.trajectory
        ax.scatter(traj.times[:-1], np.log10(traj.norms[:-1]), color=METHOD_COLOR[m],
                   marker=METHOD_MARKER[m], s=36, zorder=3, edgecolor=BACKGROUND, label=m)
        steps = 1 + np.arange(len(r.refine.norms))
        ax.plot(steps, np.log10(r.refine.norms), color=METHOD_COLOR[m], marker=METHOD_MARKER[m],
                markeredgecolor=BACKGROUND, label=f"N-{m}")
    ax.set(xlabel="Iteration", ylabel=r"$\log_{10}\|g(x)\|_\infty$", title=title)
    ax.legend(ncol=2, fontsize=8)
    savefig(fig, filename)


def main():
    problem = TutorialProblem()

    nr_runs = {eps: newton_raphson(problem, np.array([1, 1 - eps]), TOL, max_it=50)
               for eps in [0.005, 0.01, 0.05]}
    fig, ax = new_figure()
    for color, (eps, r) in zip(COLORS, nr_runs.items()):
        ax.plot(np.log10(r.norms), marker="o", color=color, markeredgecolor=BACKGROUND,
                label=rf"$\epsilon$={eps}")
    ax.set(xlabel="Iteration", ylabel=r"$\log_{10}\|g(x)\|_\infty$", title="Fig. 1(a) NR solver")
    ax.legend()
    savefig(fig, "fig1a_tutorial_nr")

    hybrid_runs = {}
    for filename, dt, K in HYBRID_RUNS:
        times = pathway.constant_step(0.01, dt)
        runs = {m: solve_hybrid(problem, problem.x0, K, times, m, "NR", TOL, max_it=50) for m in METHODS}
        hybrid_runs[(dt, K)] = runs
        figure = filename.replace("fig", "Fig. ")
        plot_hybrid(runs, f"{figure}: dt = {dt}, K = {K}", f"{filename}_tutorial_hybrid")

    def root(x):
        return f"({x[0]:.4f}, {x[1]:.4f})"

    rows = [[f"NR only, eps = {eps}", "-", "-", r.iterations, FIG1A_ITERS[eps], root(r.x)]
            for eps, r in nr_runs.items()]
    for (dt, K), runs in hybrid_runs.items():
        for m, r in runs.items():
            rows.append([f"N-{m}", dt, K, r.refine.iterations, FIG1_NR_ITERS[(dt, K)][m], root(r.refine.x)])
    write_table("fig1_fig2_tutorial_iterations",
                ["run", "dt", "K", "NR iterations (ours)", "NR iterations (paper plot)", "root"], rows,
                "Tutorial: NR iterations to reach ||g|| < 1e-5")

    rows = []
    for (dt, K), runs in hybrid_runs.items():
        for m, r in runs.items():
            points = [f"{t:g}: {n:.3g}" for t, n in zip(r.trajectory.times, r.trajectory.norms)]
            rows.append([m, dt, K] + pad(points, 10, ""))
    write_table("fig1_fig2_tutorial_pathways", ["method", "dt", "K"] + [f"point {i}" for i in range(10)],
                rows, "Tutorial: ||g(x(t_k))|| along each path (t: norm)")


if __name__ == "__main__":
    main()
