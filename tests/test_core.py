"""Sanity and regression tests. Run with:  python -m pytest tests"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dynhomotopy.datasets import load_problem  # noqa: E402
from dynhomotopy.homotopy import integrate, pathway  # noqa: E402
from dynhomotopy.hybrid import solve_hybrid  # noqa: E402
from dynhomotopy.problem import TutorialProblem  # noqa: E402
from dynhomotopy.solvers import fast_decoupled_xb, newton_raphson  # noqa: E402


@pytest.mark.parametrize("name", ["case69limit", "case_ACTIVSg500limit"])
def test_jacobian_matches_finite_differences(name):
    pf = load_problem(name)
    rng = np.random.default_rng(0)
    x = pf.flat_start() + 0.05 * rng.standard_normal(pf.n)
    J = pf.jacobian(x).toarray()
    h = 1e-7
    for j in rng.choice(pf.n, 8, replace=False):
        e = np.zeros(pf.n)
        e[j] = h
        fd = (pf.g(x + e) - pf.g(x - e)) / (2 * h)
        assert np.allclose(J[:, j], fd, atol=1e-5 * max(1, np.abs(fd).max()))


def test_pathway_construction():
    assert pathway.constant_step(0.005, 0.5, t2=0.01) == [0.0, 0.005, 0.01, 0.51, 1.0]
    assert pathway.constant_step(0.005, 1.0) == [0.0, 0.005, 1.0]
    assert len(pathway.constant_step(0.005, 0.1)) == 12


def test_be_first_step_equals_linear_approximation():
    """Section 3.2: eq. (25) (BE, one FPI) and eq. (27) (linearised static homotopy) agree."""
    pf = load_problem("case_ACTIVSg2000limit")
    be = integrate(pf, pf.flat_start(), 1e-4, [0, 0.005], "BE", None)
    lin = integrate(pf, pf.flat_start(), 1e-4, [0, 0.005], "LIN", None)
    assert np.allclose(be.states[-1], lin.states[-1], atol=1e-10)


def test_tutorial_nr_iterations():
    p = TutorialProblem()
    r = newton_raphson(p, np.array([1.0, 0.995]), tol=1e-5, max_it=50)
    assert r.converged and r.iterations == 11
    assert np.allclose(r.x, [1.5, 0.5], atol=1e-6)


def test_reference_nr_and_flat_start():
    pf = load_problem("case_ACTIVSg2000limit")
    assert newton_raphson(pf, pf.case_start()).iterations == 3     # Table 7, NR-MAT
    assert newton_raphson(pf, pf.flat_start()).iterations == 5     # Table 7, NR-flat


def test_table4_be_row_2000limit():
    """Table 4, case2000limit: 54, 3.9, 2, 2.8, 1.5 | 0.01, 8e-7, 9e-13."""
    pf = load_problem("case_ACTIVSg2000limit")
    res = solve_hybrid(pf, pf.flat_start(), 1e-4, [0, 0.005, 0.01, 0.51, 1.0], "BE", "NR")
    got = res.trajectory.norms + res.refine.norms[1:]
    paper = [54, 3.9, 2, 2.8, 1.5, 0.01, 8e-7, 9e-13]
    assert res.refine.iterations == 3
    for g, p in zip(got[:6], paper[:6]):
        assert abs(np.log10(g) - np.log10(p)) < 0.1


def test_fdxb_counts_half_iterations():
    pf = load_problem("case_ACTIVSg2000limit")
    res = solve_hybrid(pf, pf.flat_start(), 1e-4, [0, 0.005, 1.0], "BE", "FDXB")
    assert res.converged and res.refine.iterations == 22           # Table 7, BE(FDXB)
    ref = fast_decoupled_xb(pf, pf.case_start())
    assert ref.converged and ref.norms[-1] < 1e-7
