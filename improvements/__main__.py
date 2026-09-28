"""Run the paper's method and a chosen set of modifications on one case.

    python -m improvements case18482 --corrector
    python -m improvements case6024 --multiplier --richardson
    python -m improvements case36964 --dt0 0.1 --K 2e-3 --corrector --multiplier
    python -m improvements path/to/mycase.mat --corrector
"""

import argparse

from .cases import ALL, load, on_reference
from .homotopies import HOMOTOPIES
from .solve import Options, solve
from .steps import STEPS


def describe(pf, res):
    if res.refine is None:
        return "homotopy failed", "-", "-", "-"
    if not res.converged:
        outcome = "NR failed"
    else:
        outcome = "solved" if on_reference(pf, res.refine.x) else "other root"
    return outcome, len(res.trajectory.times) - 1, res.refine.iterations, res.factorizations


def main(argv=None):
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
    args = parser.parse_args(argv)
    if args.case not in ALL and not args.case.endswith(".mat"):
        parser.error(f"unknown case {args.case}; use one of {', '.join(ALL)} or a .mat file")

    options = Options(multiplier=args.multiplier, corrector=args.corrector, adaptive=args.adaptive,
                      richardson=args.richardson, richardson_tol=args.richardson_tol,
                      step=args.step, homotopy=args.homotopy)
    pf = load(args.case)
    path = [0.0, args.dt0, 1.0]
    print(f"{args.case}: dt0 = {args.dt0}, K = {args.K}, refiner = {args.refiner}\n")
    print(f"{'method':32s} {'outcome':12s} {'path steps':>10s} {'NR its':>7s} {'LUs':>5s}")
    for label, opts in (("paper", Options()), (options.name(), options)):
        res = solve(pf, pf.flat_start(), args.K, path, opts, refiner=args.refiner)
        outcome, steps, its, lus = describe(pf, res)
        print(f"{label:32s} {outcome:12s} {steps!s:>10s} {its!s:>7s} {lus!s:>5s}")


if __name__ == "__main__":
    main()
