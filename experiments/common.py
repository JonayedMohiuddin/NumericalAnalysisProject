import csv
import time
from pathlib import Path

import numpy as np

from dynhomotopy.datasets import SYSTEMS

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

# values used by the paper from Section 4.3 on
K_PAPER = 1e-4
DT0_PAPER = 0.005


def label(name):
    return SYSTEMS[name].label


def fmt(v, digits=2):
    """Format a number roughly the way the paper prints it (532, 7.2, 2e-04)."""
    if v is None:
        return "-"
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return str(v)
    if not np.isfinite(v):
        return "inf"
    a = abs(v)
    if 1e-3 <= a < 100 and len(repr(float(v))) <= 5:
        return repr(float(v))
    if a != 0 and (a < 1e-2 or a >= 1e4):
        return f"{v:.{digits - 1}e}"
    if a >= 100:
        return f"{v:.0f}"
    return f"{v:.{digits}g}"


def pad(values, length, filler="-"):
    return list(values[:length]) + [filler] * max(0, length - len(values))


def markdown_cell(text):
    # a "|" inside a cell (as in ||g||) would start a new column
    return text.replace("|", "\\|")


def write_table(name, header, rows, caption=""):
    """Save a table as results/tables/<name>.csv and .md and print it."""
    rows = [[fmt(c) for c in row] for row in rows]
    with open(TABLES / f"{name}.csv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)

    lines = [f"**{caption}**", ""] if caption else []
    lines.append("| " + " | ".join(markdown_cell(h) for h in header) + " |")
    lines.append("|" + "|".join("---" for _ in header) + "|")
    lines += ["| " + " | ".join(markdown_cell(c) for c in row) + " |" for row in rows]
    text = "\n".join(lines) + "\n"
    (TABLES / f"{name}.md").write_text(text, encoding="utf-8")
    print(text)


def benchmark(functions, reps):
    """Median run time of each function.

    The functions are run in turn, one of each per round, so that changes in
    CPU speed during the run affect all of them equally.
    """
    samples = {key: [] for key in functions}
    for _ in range(reps):
        for key, fn in functions.items():
            start = time.perf_counter()
            fn()
            samples[key].append(time.perf_counter() - start)
    return {key: float(np.median(s)) for key, s in samples.items()}


def savefig(fig, name):
    path = FIGURES / f"{name}.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    print(f"saved {path.relative_to(ROOT)}")
