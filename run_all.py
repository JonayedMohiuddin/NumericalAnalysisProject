"""Run the experiments that reproduce the paper's figures and tables.

    python run_all.py              # all of them
    python run_all.py --only 3 7   # only experiments 3 and 7
    python run_all.py --reps 20    # more repetitions for the timing tables

A single experiment can also be run on its own, e.g.
    python -m experiments.exp03_tables2to4
"""

import argparse
import importlib
import time

EXPERIMENTS = [
    ("exp01_tutorial", "Figures 1 and 2 (tutorial example)"),
    ("exp02_first_step", "Section 4.2.1 (first time step)"),
    ("exp03_tables2to4", "Tables 2, 3 and 4"),
    ("exp04_table5", "Table 5"),
    ("exp05_sec44_pathways", "Section 4.4 (paths with dt0 = 0.05, K = 0.001)"),
    ("exp06_fig3_table6", "Figure 3 and Table 6"),
    ("exp07_fig4", "Figure 4"),
    ("exp08_fig5", "Figure 5"),
    ("exp09_table7", "Table 7"),
    ("exp10_table8", "Table 8"),
    ("exp11_improvements", "Our improvements, every on/off combination"),
    ("exp12_spectrum", "Proposal: spectrum of the Jacobian"),
    ("exp13_feasibility", "Proposal: feasible region in the (K, dt0) plane"),
    ("exp14_tuning", "Proposal: golden-section tuning and a spectral rule for delta"),
    ("exp15_scratch_lu", "Proposal: NR and FDXB with our own LU"),
]
TIMED = {"exp06_fig3_table6", "exp10_table8"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", nargs="*", type=int, help="experiment numbers to run (1-15)")
    parser.add_argument("--reps", type=int, default=5, help="repetitions for the timing tables")
    args = parser.parse_args()

    for i, (module, title) in enumerate(EXPERIMENTS, start=1):
        if args.only and i not in args.only:
            continue
        print(f"\n[{i}] {title}")
        start = time.perf_counter()
        experiment = importlib.import_module(f"experiments.{module}")
        if module in TIMED:
            experiment.main(reps=args.reps)
        else:
            experiment.main()
        print(f"[{i}] done in {time.perf_counter() - start:.1f} s")


if __name__ == "__main__":
    main()
