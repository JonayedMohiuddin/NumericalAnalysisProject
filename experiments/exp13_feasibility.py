"""Proposal: where in the (K, dt0) plane does the paper's method work?

The paper's method (BE path {0, dt0, 1}, then NR) is run on a grid of K and
dt0 for four cases. Each cell is solved (reference operating point), other
root, or failed. Lines of constant delta = K / dt0 are diagonals of the map.
"""

import numpy as np
from matplotlib.colors import ListedColormap

from improvements import solve
from improvements.cases import load, on_reference

from .common import savefig, write_table
from .plots import COLORS, new_figure

CASES = ["case18482", "case36964", "case6024", "case_ACTIVSg2000limit"]
DT0 = [0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2]
K = list(np.logspace(-6, -1, 11))
PAPER_POINTS = [(0.005, 1e-4), (0.05, 1e-3), (0.1, 2e-3)]
FAILED, OTHER_ROOT, SOLVED = 0, 1, 2


def outcome(pf, k, dt0):
    res = solve(pf, pf.flat_start(), k, [0.0, dt0, 1.0])
    if not res.converged:
        return FAILED
    return SOLVED if on_reference(pf, res.refine.x) else OTHER_ROOT


def plot(ax, grid, name):
    cmap = ListedColormap(["#e4e3df", COLORS[3], COLORS[2]])
    ax.imshow(grid, origin="lower", cmap=cmap, vmin=0, vmax=2, aspect="auto",
              extent=(np.log10(K[0]) - 0.25, np.log10(K[-1]) + 0.25, -0.5, len(DT0) - 0.5))
    for dt0, k in PAPER_POINTS:
        ax.plot(np.log10(k), DT0.index(dt0), "o", color=COLORS[6], markersize=8)
    log_k = np.linspace(np.log10(K[0]), np.log10(K[-1]), 50)
    rows = np.interp(log_k - np.log10(0.02), np.log10(DT0), np.arange(len(DT0)), left=np.nan, right=np.nan)
    ax.plot(log_k, rows, "--", color=COLORS[6], lw=1)
    ax.set_yticks(range(len(DT0)), [f"{d:g}" for d in DT0])
    ax.set(xlabel="log10 K", ylabel="dt0", title=name)
    ax.grid(False)


def main(cases=CASES):
    fig, axes = new_figure(2, 2, width=11, height=9)
    fig.subplots_adjust(hspace=0.35)
    rows = []
    for ax, name in zip(axes.ravel(), cases):
        pf = load(name)
        grid = np.array([[outcome(pf, k, dt0) for k in K] for dt0 in DT0])
        plot(ax, grid, name)
        solved = grid == SOLVED
        delta = np.array([[k / dt0 for k in K] for dt0 in DT0])
        rows.append([name, f"{solved.sum()} of {grid.size}", int((grid == OTHER_ROOT).sum()),
                     f"{delta[solved].min():.2g} to {delta[solved].max():.2g}" if solved.any() else "-",
                     "yes" if grid[DT0.index(0.005), K.index(min(K, key=lambda k: abs(k - 1e-4)))] == SOLVED else "no"])
    fig.suptitle("Grey: fails, yellow: other root, green: reference solution. "
                 "Dots: the paper's settings, dashed: delta = 0.02", fontsize=10)
    savefig(fig, "feasibility_maps")
    write_table("feasibility", ["case", "cells solved", "cells on another root", "delta range that solves",
                                "paper's (dt0, K) solves"], rows,
                "Feasible region of the paper's method in the (K, dt0) plane")


if __name__ == "__main__":
    main()
