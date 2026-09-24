"""Shared code for Tables 2-5, which list ||g|| along the path and then for 3 NR iterations."""

from dynhomotopy.datasets import load_problem
from dynhomotopy.hybrid import solve_hybrid

from .common import K_PAPER, fmt, pad

NR_SHOWN = 3


def norm_row(name, times, method, K=K_PAPER):
    pf = load_problem(name)
    res = solve_hybrid(pf, pf.flat_start(), K, times, method, "NR")
    path = pad(res.trajectory.norms, len(times), "diverged")
    nr = res.refine.norms[1:] if res.refine else []
    return path + pad(nr, NR_SHOWN), res


def status(res):
    if res.refine is None:
        return "homotopy failed"
    return f"yes ({res.refine.iterations})" if res.converged else "no"


def ours_vs_paper(ours, paper):
    if paper is None:
        return [fmt(v) for v in ours]
    return [f"{fmt(o)} / {fmt(p)}" for o, p in zip(ours, paper)]


def header(times):
    return [f"t={t:g}" for t in times] + [f"NR {i}" for i in range(1, NR_SHOWN + 1)]
