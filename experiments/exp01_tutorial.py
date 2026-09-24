"""Section 3.3 - generic tutorial example: Figures 1(a)-(d) and 2.

(a) NR alone from x0 = [1, 1 - eps] (the Jacobian is singular at [1, 1]).
(b)-(d) dynamic homotopy from x0 = [1, 1] with dt0 = 0.01, dt = 0.5 and
        K = 0.05, 0.01, 0.005, followed by NR (curves N-FE, N-BE, N-RK2).
Fig. 2  same with dt = 0.125 and K = 0.05.
The tutorial uses tolerance 1e-5. Homotopy points are drawn at their time
t_k in [0, 1]; NR iteration i is drawn at abscissa 1 + i.
"""

from __future__ import annotations

import numpy as np

from common import pad, savefig, write_table
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.problem import TutorialProblem
from dynhomotopy.solvers import newton_raphson
from paper_values import FIG1_NR_ITERS, FIG1A_ITERS
from plots import METHOD_COLOR, METHOD_MARKER, SERIES, new_figure

TOL = 1e-5
METHODS = ["FE", "BE", "RK2"]


def nr_alone(problem, epsilons):
    return {eps: newton_raphson(problem, np.array([1.0, 1.0 - eps]), TOL, max_it=50) for eps in epsilons}


def hybrid_runs(problem, dt, K):
    times = pathway.constant_step(0.01, dt)
    return {m: solve_hybrid(problem, problem.x0, K, times, m, "NR", TOL, max_it=50) for m in METHODS}


def draw_hybrid(ax, runs, title):
    for m, r in runs.items():
        tr, nr = r.trajectory, r.refine
        ax.scatter(tr.times[:-1], np.log10(tr.norms[:-1]), color=METHOD_COLOR[m], marker=METHOD_MARKER[m],
                   s=36, zorder=3, edgecolor="#fcfcfb", linewidth=1, label=m)
        xs = 1 + np.arange(len(nr.norms))
        ax.plot(xs, np.log10(nr.norms), color=METHOD_COLOR[m], marker=METHOD_MARKER[m],
                markeredgecolor="#fcfcfb", label=f"N-{m}")
    ax.set(xlabel="Iteration", ylabel=r"$\log_{10}\|g(x)\|_\infty$", title=title)
    ax.legend(ncol=2, fontsize=8)


def main():
    problem = TutorialProblem()

    # ---- Figure 1(a) -------------------------------------------------------------
    eps_runs = nr_alone(problem, [0.005, 0.01, 0.05])
    fig, ax = new_figure()
    for color, (eps, r) in zip(SERIES, eps_runs.items()):
        ax.plot(np.arange(len(r.norms)), np.log10(r.norms), marker="o", color=color,
                markeredgecolor="#fcfcfb", label=rf"$\epsilon$={eps}")
    ax.set(xlabel="Iteration", ylabel=r"$\log_{10}\|g(x)\|_\infty$", title="Fig. 1(a) NR solver")
    ax.legend()
    savefig(fig, "fig1a_tutorial_nr")

    # ---- Figures 1(b)-(d) and 2 ----------------------------------------------------
    configs = [("fig1b", 0.5, 0.05), ("fig1c", 0.5, 0.01), ("fig1d", 0.5, 0.005), ("fig2", 0.125, 0.05)]
    all_runs = {}
    for name, dt, K in configs:
        runs = hybrid_runs(problem, dt, K)
        all_runs[(dt, K)] = runs
        fig, ax = new_figure()
        draw_hybrid(ax, runs, f"{name.replace('fig', 'Fig. ')}: dt = {dt}, K = {K}")
        savefig(fig, f"{name}_tutorial_hybrid")

    # ---- summary table: NR iterations, ours vs paper ----------------------------------
    rows = [[f"NR only, eps = {eps}", "-", "-", r.iterations, FIG1A_ITERS[eps],
             f"({r.x[0]:.4f}, {r.x[1]:.4f})"] for eps, r in eps_runs.items()]
    for (dt, K), runs in all_runs.items():
        for m, r in runs.items():
            tr = r.trajectory
            rows.append([f"N-{m}", dt, K, r.refine.iterations, FIG1_NR_ITERS[(dt, K)][m],
                         f"({r.refine.x[0]:.4f}, {r.refine.x[1]:.4f})"])
    write_table("fig1_fig2_tutorial_iterations",
                ["run", "dt", "K", "NR iterations (ours)", "NR iterations (paper, read from plot)", "root"],
                rows, "Tutorial (Sec. 3.3): NR iterations to reach ||g|| < 1e-5")

    norm_rows = []
    for (dt, K), runs in all_runs.items():
        for m, r in runs.items():
            norm_rows.append([m, dt, K] + pad([f"{t:g}: {n:.3g}" for t, n in zip(r.trajectory.times, r.trajectory.norms)], 10, ""))
    write_table("fig1_fig2_tutorial_pathways", ["method", "dt", "K"] + [f"point {i}" for i in range(10)], norm_rows,
                "Tutorial: ||g(x(t_k))|| along each pathway (t: norm)")


if __name__ == "__main__":
    main()
