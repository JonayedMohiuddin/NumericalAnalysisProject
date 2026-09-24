"""Norm-evolution tables (Tables 2-5): ||g|| along the pathway, then 3 NR iterations."""

from __future__ import annotations

from common import K_PAPER, fmt, label, pad
from dynhomotopy.datasets import load_problem
from dynhomotopy.hybrid import solve_hybrid

NR_SHOWN = 3  # the paper tabulates the first three NR iterations


def norm_row(name: str, times: list[float], method: str, K: float = K_PAPER):
    """[||g(x(t_k))|| for every t_k] + [NR iteration 1..3] and the hybrid result."""
    pf = load_problem(name)
    res = solve_hybrid(pf, pf.flat_start(), K, times, method, "NR")
    traj = pad(res.trajectory.norms, len(times), "diverged")
    nr = res.refine.norms[1:] if res.refine else []
    return traj + pad(nr, NR_SHOWN), res


def status(res) -> str:
    if res.refine is None:
        return "homotopy failed"
    return f"yes ({res.refine.iterations})" if res.converged else "no"


def comparison_rows(ours: list, paper: list | None) -> list[str]:
    """'ours / paper' cells."""
    if paper is None:
        return [fmt(v) for v in ours]
    return [f"{fmt(o)} / {fmt(p)}" for o, p in zip(ours, paper)]


def header(times: list[float]) -> list[str]:
    return [f"t={t:g}" for t in times] + [f"NR {i}" for i in range(1, NR_SHOWN + 1)]


def case_label(name: str) -> str:
    return label(name)
