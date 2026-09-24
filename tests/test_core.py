import numpy as np
import pytest

from dynhomotopy.datasets import load_problem
from dynhomotopy.homotopy import integrate, pathway
from dynhomotopy.hybrid import solve_hybrid
from dynhomotopy.problem import TutorialProblem
from dynhomotopy.solvers import fast_decoupled_xb, newton_raphson


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


def test_pathway():
    assert pathway.constant_step(0.005, 0.5, t2=0.01) == [0.0, 0.005, 0.01, 0.51, 1.0]
    assert pathway.constant_step(0.005, 1.0) == [0.0, 0.005, 1.0]
    assert len(pathway.constant_step(0.005, 0.1)) == 12


def test_first_step_be_equals_linear_approximation():
    # eqs. (25) and (27) give the same point
    pf = load_problem("case_ACTIVSg2000limit")
    be = integrate(pf, pf.flat_start(), 1e-4, [0, 0.005], "BE", None)
    lin = integrate(pf, pf.flat_start(), 1e-4, [0, 0.005], "LIN", None)
    assert np.allclose(be.states[-1], lin.states[-1], atol=1e-10)


def test_tutorial_nr():
    r = newton_raphson(TutorialProblem(), np.array([1.0, 0.995]), tol=1e-5, max_it=50)
    assert r.converged and r.iterations == 11
    assert np.allclose(r.x, [1.5, 0.5], atol=1e-6)


def test_nr_iterations_match_table7():
    pf = load_problem("case_ACTIVSg2000limit")
    assert newton_raphson(pf, pf.case_start()).iterations == 3
    assert newton_raphson(pf, pf.flat_start()).iterations == 5


def test_table4_row_for_2000_bus_case():
    pf = load_problem("case_ACTIVSg2000limit")
    res = solve_hybrid(pf, pf.flat_start(), 1e-4, [0, 0.005, 0.01, 0.51, 1.0], "BE", "NR")
    ours = res.trajectory.norms + res.refine.norms[1:]
    paper = [54, 3.9, 2, 2.8, 1.5, 0.01]
    assert res.refine.iterations == 3
    for a, b in zip(ours, paper):
        assert abs(np.log10(a) - np.log10(b)) < 0.1


def test_fdxb():
    pf = load_problem("case_ACTIVSg2000limit")
    res = solve_hybrid(pf, pf.flat_start(), 1e-4, [0, 0.005, 1.0], "BE", "FDXB")
    assert res.converged and res.refine.iterations == 22
    assert fast_decoupled_xb(pf, pf.case_start()).converged
