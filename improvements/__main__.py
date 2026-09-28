"""Run the paper's method and a chosen set of modifications on one case,
or one of the studies (spectrum, feasibility, tune, scratch).

    python -m improvements case18482 --corrector
    python -m improvements case6024 --multiplier --richardson
    python -m improvements case36964 --dt0 0.1 --K 2e-3 --corrector --multiplier
    python -m improvements path/to/mycase.mat --corrector

    python -m improvements spectrum case18482
    python -m improvements feasibility case6024
    python -m improvements tune case6024
    python -m improvements scratch case_ACTIVSg500limit --refiner FDXB
"""

import argparse
import sys

import numpy as np

from . import scratch
from .cases import ALL, load, on_reference
from .feasibility import DT0_VALUES, K_VALUES, feasibility_grid
from .homotopies import HOMOTOPIES
from .solve import Options, solve
from .spectrum import spectrum_summary
from .steps import STEPS
from .tuning import tune


def describe(pf, res):
    if res.refine is None:
        return "homotopy failed", "-", "-", "-"
    if not res.converged:
        outcome = "NR failed"
    else:
        outcome = "solved" if on_reference(pf, res.refine.x) else "other root"
    return outcome, len(res.trajectory.times) - 1, res.refine.iterations, res.factorizations


def run_solve(pf, args):
    options = Options(multiplier=args.multiplier, corrector=args.corrector, adaptive=args.adaptive,
                      richardson=args.richardson, richardson_tol=args.richardson_tol,
                      step=args.step, homotopy=args.homotopy)
    path = [0.0, args.dt0, 1.0]
    print(f"{args.case}: dt0 = {args.dt0}, K = {args.K}, refiner = {args.refiner}\n")
    print(f"{'method':32s} {'outcome':12s} {'path steps':>10s} {'NR its':>7s} {'LUs':>5s}")
    for label, opts in (("paper", Options()), (options.name(), options)):
        res = solve(pf, pf.flat_start(), args.K, path, opts, refiner=args.refiner)
        outcome, steps, its, lus = describe(pf, res)
        print(f"{label:32s} {outcome:12s} {steps!s:>10s} {its!s:>7s} {lus!s:>5s}")


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


STUDIES = {"spectrum": run_spectrum, "feasibility": run_feasibility, "tune": run_tune,
           "scratch": run_scratch}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    study = argv.pop(0) if argv and argv[0] in STUDIES else None

    parser = argparse.ArgumentParser(prog="python -m improvements", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", help="test system (e.g. case18482) or a path to a .mat file")
    parser.add_argument("--multiplier", action="store_true", help="optimal multiplier in the final NR")
    parser.add_argument("--corrector", action="store_true", help="Newton corrector after every step")
    control = parser.add_mutually_exclusive_group()
    control.add_argument("--adaptive", action="store_true", help="halve steps while the mismatch grows")
    control.add_argument("--richardson", action="store_true", help="Richardson error control of the steps")
    parser.add_argument("--richardson-tol", type=float, default=0.5)
    parser.add_argument("--step", choices=list(STEPS), default="BE")
    parser.add_argument("--homotopy", choices=list(HOMOTOPIES), default="fpv")
    parser.add_argument("--dt0", type=float, default=0.005, help="first time step (default 0.005)")
    parser.add_argument("--K", type=float, default=1e-4, help="factor K (default 1e-4)")
    parser.add_argument("--refiner", choices=["NR", "FDXB"], default="NR")
    parser.add_argument("--delta", type=float, default=0.02, help="shift for 'spectrum' (default 0.02)")
    args = parser.parse_args(argv)
    if args.case not in ALL and not args.case.endswith(".mat"):
        parser.error(f"unknown case {args.case}; use one of {', '.join(ALL)} or a .mat file")

    run = STUDIES[study] if study else run_solve
    run(load(args.case), args)


if __name__ == "__main__":
    main()
