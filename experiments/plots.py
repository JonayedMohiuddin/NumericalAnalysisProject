import matplotlib.pyplot as plt

plt.switch_backend("Agg")

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
BACKGROUND = "#fcfcfb"

METHOD_COLOR = {"FE": COLORS[0], "BE": COLORS[1], "RK2": COLORS[2]}
METHOD_MARKER = {"FE": "o", "BE": "s", "RK2": "^"}
CASE_LINESTYLE = {"109k": "-", "36k": ":", "lim": "--"}

plt.rcParams.update({
    "figure.facecolor": BACKGROUND,
    "axes.facecolor": BACKGROUND,
    "axes.edgecolor": "#52514e",
    "xtick.color": "#52514e",
    "ytick.color": "#52514e",
    "axes.grid": True,
    "grid.color": "#e4e3df",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 2,
    "lines.markersize": 7,
    "legend.frameon": False,
    "font.size": 10,
})


def new_figure(nrows=1, ncols=1, width=6.4, height=4.2, **kwargs):
    return plt.subplots(nrows, ncols, figsize=(width, height), **kwargs)
