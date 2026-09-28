"""Proposal: where in the (K, dt0) plane does the paper's method work?

The paper's method (BE path {0, dt0, 1}, then NR) is run on a grid of K and
dt0 for four cases (improvements/feasibility.py). Each cell is solved
(reference operating point), other root, or failed. Lines of constant
delta = K / dt0 are diagonals of the map.
"""

import numpy as np
from matplotlib.colors import ListedColormap

from improvements.feasibility import DT0_VALUES, K_VALUES, OTHER_ROOT, SOLVED, feasibility_grid
from improvements.cases import load

from .common import savefig, write_table
from .plots import COLORS, new_figure

CASES = ["case18482", "case36964", "case6024", "case_ACTIVSg2000limit"]
PAPER_POINTS = [(0.005, 1e-4), (0.05, 1e-3), (0.1, 2e-3)]


def plot(ax, grid, name):
    cmap = ListedColormap(["#e4e3df", COLORS[3], COLORS[2]])
    ax.imshow(grid, origin="lower", cmap=cmap, vmin=0, vmax=2, aspect="auto",
              extent=(np.log10(K_VALUES[0]) - 0.25, np.log10(K_VALUES[-1]) + 0.25,
                      -0.5, len(DT0_VALUES) - 0.5))
    for dt0, k in PAPER_POINTS:
        ax.plot(np.log10(k), DT0_VALUES.index(dt0), "o", color=COLORS[6], markersize=8)
    log_k = np.linspace(np.log10(K_VALUES[0]), np.log10(K_VALUES[-1]), 50)
    rows = np.interp(log_k - np.log10(0.02), np.log10(DT0_VALUES), np.arange(len(DT0_VALUES)),
                     left=np.nan, right=np.nan)
    ax.plot(log_k, rows, "--", color=COLORS[6], lw=1)
    ax.set_yticks(range(len(DT0_VALUES)), [f"{d:g}" for d in DT0_VALUES])
    ax.set(xlabel="log10 K", ylabel="dt0", title=name)
    ax.grid(False)


def main(cases=CASES):
    fig, axes = new_figure(2, 2, width=11, height=9)
    fig.subplots_adjust(hspace=0.35)
    rows = []
    for ax, name in zip(axes.ravel(), cases):
        grid = feasibility_grid(load(name))
        plot(ax, grid, name)
        solved = grid == SOLVED
        delta = np.array([[k / dt0 for k in K_VALUES] for dt0 in DT0_VALUES])
        paper_cell = grid[DT0_VALUES.index(0.005), int(np.argmin(np.abs(np.log10(K_VALUES) + 4)))]
        rows.append([name, f"{solved.sum()} of {grid.size}", int((grid == OTHER_ROOT).sum()),
                     f"{delta[solved].min():.2g} to {delta[solved].max():.2g}" if solved.any() else "-",
                     "yes" if paper_cell == SOLVED else "no"])
    fig.suptitle("Grey: fails, yellow: other root, green: reference solution. "
                 "Dots: the paper's settings, dashed: delta = 0.02", fontsize=10)
    savefig(fig, "feasibility_maps")
    write_table("feasibility", ["case", "cells solved", "cells on another root", "delta range that solves",
                                "paper's (dt0, K) solves"], rows,
                "Feasible region of the paper's method in the (K, dt0) plane")


if __name__ == "__main__":
    main()
