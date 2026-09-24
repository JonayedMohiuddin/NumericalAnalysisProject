"""Plot style shared by all figure scripts.

Categorical colours come from a fixed-order, colour-vision-deficiency-checked
palette, and every method also gets its own marker so identity never relies
on colour alone.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
TEXT, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

METHOD_COLOR = {"FE": SERIES[0], "BE": SERIES[1], "RK2": SERIES[2], "NR": SERIES[6]}
METHOD_MARKER = {"FE": "o", "BE": "s", "RK2": "^", "NR": "D"}
CASE_STYLE = {"109k": "-", "36k": ":", "lim": "--"}

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": MUTED, "axes.labelcolor": TEXT, "text.color": TEXT,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2, "lines.markersize": 7,
    "legend.frameon": False, "font.size": 10, "axes.titlesize": 11,
})


def new_figure(nrows: int = 1, ncols: int = 1, width: float = 6.4, height: float = 4.2, **kw):
    return plt.subplots(nrows, ncols, figsize=(width, height), **kw)
