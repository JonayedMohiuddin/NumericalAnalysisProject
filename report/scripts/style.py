"""Shared matplotlib style for the report figures."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
REPORT = REPO / "report"
FIG = REPORT / "figures"
GEN = REPORT / "generated"
TABLES = REPO / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
GEN.mkdir(parents=True, exist_ok=True)

MAROON = "#1D3557"  # navy, used for "the paper" reference marks
BLUE = "#3A7CA5"
CORAL = "#E07A5F"
GREEN = "#3D9970"
GOLD = "#5DADE2"  # sky blue (replaces the old gold accent)
SKY = "#2E90D1"
PALE = "#A9D6F5"
NAVY2 = "#5B7DB1"
PURPLE = "#6C5B7B"
GREY = "#8D8A85"
LIGHT = "#E9E6E1"
INK = "#2E2B28"
BG = "#FFFFFF"

METHOD = {"FE": BLUE, "RK2": GREEN, "BE": CORAL, "LIN": GOLD}

plt.rcParams.update({
    "figure.facecolor": BG,
    "axes.facecolor": BG,
    "axes.edgecolor": "#5b5752",
    "axes.labelcolor": INK,
    "xtick.color": "#5b5752",
    "ytick.color": "#5b5752",
    "axes.grid": True,
    "grid.color": "#ECEAE6",
    "grid.linewidth": 0.8,
    "axes.axisbelow": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.titlecolor": INK,
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "lines.linewidth": 1.8,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "pdf.fonttype": 42,
})


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)
    print("figure", name)


def write_tex(name, text):
    (GEN / f"{name}.tex").write_text(text, encoding="utf-8")
    print("table ", name)
