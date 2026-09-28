"""Golden-section search for delta = K / dt0 and for dt0, minimising LU factorizations."""

import math

import numpy as np

from .solve import Options, solve

INV_PHI = (math.sqrt(5) - 1) / 2
FAIL_COST = 50.0


def golden_section(f, a, b, tol=0.05, max_evals=25):
    """Minimise f on [a, b]. Assumes f is unimodal there, but always returns the best point seen.

    Returns (best x, best f, number of evaluations).
    """
    seen = {}

    def value(x):
        if x not in seen:
            seen[x] = f(x)
        return seen[x]

    c, d = b - INV_PHI * (b - a), a + INV_PHI * (b - a)
    while b - a > tol and len(seen) < max_evals:
        if value(c) <= value(d):
            b, d = d, c
            c = b - INV_PHI * (b - a)
        else:
            a, c = c, d
            d = a + INV_PHI * (b - a)
    best = min(seen, key=seen.get)
    return best, seen[best], len(seen)


def lu_cost(pf, delta, dt0, options=Options()):
    """Objective for the search: total LUs of homotopy + NR for the path {0, dt0, 1}.

    A small tie-break, 0.1 log10 ||g(x(1))||, lets the search move inside a
    plateau of equal LU counts. Failed runs cost FAIL_COST plus the tie-break,
    so the search still prefers smaller mismatches among failures.
    """
    res = solve(pf, pf.flat_start(), delta * dt0, [0.0, dt0, 1.0], options)
    g1 = res.trajectory.norms[-1] if not res.trajectory.failed else 1e6
    tie_break = 0.1 * np.log10(max(g1, 1e-16))
    return (res.factorizations if res.converged else FAIL_COST) + tie_break


def tune(pf, dt0=0.005, log_delta=(-4.0, 1.0), log_dt0=(-3.0, -0.5)):
    """Tune log10(delta) with dt0 fixed, then log10(dt0) with that delta fixed.

    Returns (delta, dt0, evaluations). Each evaluation is one full solve, so
    the search itself costs roughly 6 LUs per evaluation.
    """
    u, _, evals_delta = golden_section(lambda u: lu_cost(pf, 10**u, dt0), *log_delta)
    delta = 10**u
    w, _, evals_dt0 = golden_section(lambda w: lu_cost(pf, delta, 10**w), *log_dt0)
    return delta, 10**w, evals_delta + evals_dt0
