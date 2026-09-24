"""Reproduce every figure and table of the paper, in the paper's order.

    python run_all.py              # everything (about 20-30 min; timing runs dominate)
    python run_all.py --only 3 7   # selected experiments
    python run_all.py --reps 20    # more repetitions for the CPU-time tables

Outputs land in results/tables (CSV + Markdown) and results/figures (PNG).
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "experiments"))

EXPERIMENTS = [
    ("exp01_tutorial", "Sec. 3.3  - Figures 1 and 2 (tutorial example)"),
    ("exp02_first_step", "Sec. 4.2.1 - the first time step t1"),
    ("exp03_tables2to4", "Sec. 4.3  - Tables 2, 3, 4 (FE / RK2 / BE, large dt)"),
    ("exp04_table5", "Sec. 4.3  - Table 5 (progressive time steps)"),
    ("exp05_sec44_pathways", "Sec. 4.4  - pathways with dt0 = 0.05, K = 0.001"),
    ("exp06_fig3_table6", "Sec. 4.4  - Figure 3 and Table 6"),
    ("exp07_fig4", "Sec. 4.5  - Figure 4 (high time step, case109272)"),
    ("exp08_fig5", "Sec. 4.6  - Figure 5 (time-step sensitivity)"),
    ("exp09_table7", "Sec. 4.7  - Table 7 (iterations of all methods)"),
    ("exp10_table8", "Sec. 4.8  - Table 8 (CPU time)"),
]
TIMED = {"exp06_fig3_table6", "exp10_table8"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", type=int, help="experiment numbers (1-10) to run")
    ap.add_argument("--reps", type=int, default=5, help="repetitions for the CPU-time tables")
    args = ap.parse_args()
    for i, (module, title) in enumerate(EXPERIMENTS, start=1):
        if args.only and i not in args.only:
            continue
        print(f"\n{'=' * 78}\n[{i}] {title}\n{'=' * 78}")
        start = time.perf_counter()
        mod = importlib.import_module(module)
        if module in TIMED:
            mod.main(reps=args.reps)
        else:
            mod.main()
        print(f"[{i}] finished in {time.perf_counter() - start:.1f} s")


if __name__ == "__main__":
    main()
