"""Run one proposal feature on one case.

    python -m extensions spectrum case18482
    python -m extensions feasibility case6024
    python -m extensions tune case6024
    python -m extensions scratch case_ACTIVSg500limit --refiner FDXB

Richardson step control, the fifth feature, is a switch of the improvements
package:  python -m improvements case6024 --richardson
"""

import argparse

import numpy as np

from improvements import solve
from improvements.cases import ALL, load, on_reference

from . import scratch
from .feasibility import DT0_VALUES, K_VALUES, feasibility_grid
from .spectrum import spectrum_summary
from .tuning import tune


def run_spectrum(pf, args):
    s = spectrum_summary(pf.jacobian(pf.flat_start()).tocsc(), args.delta)
    print(f"largest eigenvalue modulus of J (power method):     {s['lambda_max']:.4g}")
    print(f"smallest eigenvalue modulus of J (inverse iteration): {s['lambda_min']:.4g}")
    print(f"smallest of J + {args.delta:g} I (shifted inverse iteration): {s['shifted_min']:.4g}")
    print(f"ratio lambda_max / lambda_min: J {s['ratio_J']:.3g}, J + {args.delta:g} I {s['ratio_shifted']:.3g}")


def run_feasibility(pf, args):
    grid = feasibility_grid(pf)
    symbols = {0: ".", 1: "o", 2: "#"}
    print("# = reference solution, o = other root, . = fails")
    print("dt0 \\ log10 K " + " ".join(f"{np.log10(k):4.1f}" for k in K_VALUES))
    for dt0, row in reversed(list(zip(DT0_VALUES, grid))):
        print(f"{dt0:<14g} " + " ".join(f"{symbols[c]:>4s}" for c in row))


def run_tune(pf, args):
    delta, dt0, evaluations = tune(pf)
    res = solve(pf, pf.flat_start(), delta * dt0, [0.0, dt0, 1.0])
    solved = res.converged and on_reference(pf, res.refine.x)
    paper = solve(pf, pf.flat_start(), 1e-4, [0.0, 0.005, 1.0])
    paper_solved = paper.converged and on_reference(pf, paper.refine.x)
    print(f"tuned delta = {delta:.3g}, dt0 = {dt0:.3g}  ({evaluations} solves in the search)")
    print(f"tuned: solved = {solved}, LUs = {res.factorizations}")
    print(f"paper (delta = 0.02, dt0 = 0.005): solved = {paper_solved}, LUs = {paper.factorizations}")


def run_scratch(pf, args):
    path = [0.0, 0.005, 1.0]
    x1, _ = scratch.backward_euler_path(pf, pf.flat_start(), 1e-4, path)
    if args.refiner == "NR":
        ours = scratch.newton_raphson(pf, x1)
    else:
        ours = scratch.fast_decoupled_xb(pf, x1)
    theirs = solve(pf, pf.flat_start(), 1e-4, path, refiner=args.refiner)
    print(f"{args.refiner} iterations: our LU {ours.iterations}, SuperLU {theirs.refine.iterations}")
    print(f"largest difference in the solution: {np.max(np.abs(ours.x - theirs.refine.x)):.1e}")
    print(f"time: our LU {ours.time:.3f} s, SuperLU {theirs.refine.time:.3f} s")


FEATURES = {"spectrum": run_spectrum, "feasibility": run_feasibility, "tune": run_tune,
            "scratch": run_scratch}


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m extensions", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("feature", choices=list(FEATURES))
    parser.add_argument("case", choices=ALL, metavar="case", help="test system, e.g. case18482")
    parser.add_argument("--delta", type=float, default=0.02, help="shift for 'spectrum' (default 0.02)")
    parser.add_argument("--refiner", choices=["NR", "FDXB"], default="NR", help="for 'scratch'")
    args = parser.parse_args(argv)
    FEATURES[args.feature](load(args.case), args)


if __name__ == "__main__":
    main()
