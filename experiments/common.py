"""Shared helpers for the experiment scripts: paths, formatting, table output, timing."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dynhomotopy.datasets import SYSTEMS  # noqa: E402

TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

K_PAPER = 1e-4      # factor K used from Section 4.3 on
DT0_PAPER = 0.005   # initial time step dt0 used from Section 4.3 on


def label(name: str) -> str:
    return SYSTEMS[name].label if name in SYSTEMS else name


def fmt(v, digits: int = 2) -> str:
    """Compact number formatting close to the paper's tables (e.g. 2e-04, 7.2, 532)."""
    if v is None:
        return "-"
    if isinstance(v, str):
        return v
    if not np.isfinite(v):
        return "inf"
    if isinstance(v, (int, np.integer)):
        return str(v)
    a = abs(v)
    if 1e-3 <= a < 100 and len(repr(float(v))) <= 5:
        return repr(float(v))  # short parameters such as 0.125 or 0.005 print exactly
    if a != 0 and (a < 1e-2 or a >= 1e4):
        return f"{v:.{digits - 1}e}"
    if a >= 100:
        return f"{v:.0f}"
    return f"{v:.{digits}g}"


def write_table(name: str, header: list[str], rows: list[list], caption: str = "") -> str:
    """Write ``results/tables/<name>.csv`` and ``.md`` and return the markdown."""
    cells = [[fmt(c) if not isinstance(c, str) else c for c in r] for r in rows]
    with open(TABLES / f"{name}.csv", "w", encoding="utf-8") as fh:
        fh.write(",".join(header) + "\n")
        for r in cells:
            fh.write(",".join(r) + "\n")
    md = []
    if caption:
        md += [f"**{caption}**", ""]
    md.append("| " + " | ".join(header) + " |")
    md.append("|" + "|".join("---" for _ in header) + "|")
    md += ["| " + " | ".join(r) + " |" for r in cells]
    text = "\n".join(md) + "\n"
    (TABLES / f"{name}.md").write_text(text, encoding="utf-8")
    print(text)
    return text


def pad(values: list, length: int, filler="-") -> list:
    """Right-pad a list to ``length`` (missing NR iterations print as '-')."""
    return list(values[:length]) + [filler] * max(0, length - len(values))


def benchmark(fns: dict, reps: int) -> dict:
    """Median wall time of every callable in ``fns`` over ``reps`` rounds.

    The paper averages 100 runs. Here the callables run round-robin (one of
    each per round) so that clock-frequency drift hits all methods alike, and
    the median discards outliers.
    """
    samples = {k: [] for k in fns}
    for _ in range(reps):
        for k, fn in fns.items():
            t = time.perf_counter()
            fn()
            samples[k].append(time.perf_counter() - t)
    return {k: float(np.median(v)) for k, v in samples.items()}


def savefig(fig, name: str) -> None:
    path = FIGURES / f"{name}.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    print(f"saved {path.relative_to(ROOT)}")
