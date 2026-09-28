import numpy as np
import pytest
import scipy.sparse as sp

from dynhomotopy.datasets import load_problem
from dynhomotopy.homotopy.fpv import FixedPointHomotopy
from dynhomotopy.homotopy.integrators import backward_euler
from dynhomotopy.hybrid import solve_hybrid
from improvements import Options, scratch, solve
from improvements.cases import load
from improvements.gauss import gauss_solve, lu_factor, lu_solve
from improvements.homotopies import NewtonHomotopy, ScaledHomotopy
from improvements.multiplier import optimal_multiplier
from improvements.spectrum import inverse_iteration, power_method
from improvements.steps import with_corrector
from improvements.tuning import golden_section

PATH = [0.0, 0.005, 1.0]


def test_all_options_off_is_the_paper_method():
    pf = load_problem("case_ACTIVSg2000limit")
    paper = solve_hybrid(pf, pf.flat_start(), 1e-4, PATH)
    ours = solve(pf, pf.flat_start(), 1e-4, PATH, Options())
    assert np.allclose(paper.trajectory.x_final, ours.trajectory.x_final)
    assert paper.refine.norms == ours.refine.norms


def test_optimal_multiplier():
    a = np.array([1.0, -2.0])
    # linear problem: the full Newton step solves it, so mu = 1
    assert optimal_multiplier(a, np.zeros(2)) == pytest.approx(1.0)
    # a step that overshoots (the mismatch flips sign and grows) gets shortened
    assert optimal_multiplier(a, -3 * a) < 1


def test_corrector_moves_the_point_closer_to_the_path():
    pf = load_problem("case_ACTIVSg2000limit")
    h = FixedPointHomotopy(pf, pf.flat_start(), 1e-4)
    plain = backward_euler(h, h.x0, 0.0, 0.005)
    corrected = with_corrector(backward_euler)(h, h.x0, 0.0, 0.005)
    residual = lambda x: np.linalg.norm(h.G(x, 0.005), np.inf)
    assert residual(corrected) < residual(plain)


def test_adaptive_keeps_the_direct_jump_when_it_works():
    pf = load_problem("case_ACTIVSg500limit")
    fixed = solve(pf, pf.flat_start(), 1e-4, PATH, Options())
    adaptive = solve(pf, pf.flat_start(), 1e-4, PATH, Options(adaptive=True))
    assert adaptive.trajectory.times == PATH
    assert np.allclose(adaptive.trajectory.x_final, fixed.trajectory.x_final)


def test_multiplier_and_corrector_solve_case14limit():
    # the paper's method needs 11 NR iterations here and fails the limit of 10
    pf = load("case14limit")
    assert not solve(pf, pf.flat_start(), 1e-4, PATH).converged
    assert solve(pf, pf.flat_start(), 1e-4, PATH, Options(multiplier=True)).converged
    assert solve(pf, pf.flat_start(), 1e-4, PATH, Options(corrector=True)).converged


def test_explicit_step_starts_with_backward_euler():
    pf = load_problem("case_ACTIVSg2000limit")
    paper = solve(pf, pf.flat_start(), 1e-4, PATH)
    rk4 = solve(pf, pf.flat_start(), 1e-4, PATH, Options(step="RK4"))
    assert np.allclose(rk4.trajectory.states[1], paper.trajectory.states[1])


@pytest.mark.parametrize("homotopy", [ScaledHomotopy, NewtonHomotopy])
def test_other_homotopies_reach_the_original_problem_at_t1(homotopy):
    pf = load_problem("case69limit")
    h = homotopy(pf, pf.flat_start(), 1e-4)
    x = pf.flat_start() + 0.01
    assert np.allclose(h.G(x, 1.0), pf.g(x))
    assert np.allclose(h.Gx(x, 1.0).toarray(), pf.jacobian(x).toarray())


def test_gauss_and_lu_match_numpy():
    rng = np.random.default_rng(3)
    a = rng.standard_normal((40, 40))
    b = rng.standard_normal(40)
    x = np.linalg.solve(a, b)
    lu, perm = lu_factor(a)
    assert np.allclose(gauss_solve(a, b), x)
    assert np.allclose(lu_solve(lu, perm, b), x)


def test_pivoting_handles_a_zero_on_the_diagonal():
    a = np.array([[0.0, 1.0], [1.0, 1.0]])
    assert np.allclose(gauss_solve(a, [1.0, 2.0]), [1.0, 1.0])


def test_power_method_and_inverse_iteration():
    a = sp.diags([0.5, 2.0, -7.0, 3.0], format="csc")
    largest, _ = power_method(lambda v: a @ v, 4)
    smallest, _ = inverse_iteration(a)
    closest_to_2_1, _ = inverse_iteration(a, shift=2.1)
    assert np.isclose(largest, 7.0, rtol=1e-4)
    assert np.isclose(smallest, 0.5, rtol=1e-4)
    assert np.isclose(closest_to_2_1, 0.1, rtol=1e-3)


def test_golden_section_finds_the_minimum_of_a_parabola():
    x, fx, _ = golden_section(lambda x: (x - 1.3) ** 2, -2, 4, tol=1e-4, max_evals=60)
    assert abs(x - 1.3) < 1e-3


def test_own_lu_gives_the_same_newton_iterates():
    pf = load_problem("case_ACTIVSg500limit")
    x1, _ = scratch.backward_euler_path(pf, pf.flat_start(), 1e-4, [0, 0.005, 1.0])
    ours = scratch.newton_raphson(pf, x1)
    theirs = solve(pf, pf.flat_start(), 1e-4, [0, 0.005, 1.0])
    assert ours.iterations == theirs.refine.iterations
    assert np.allclose(ours.x, theirs.refine.x, atol=1e-10)


def test_richardson_reaches_t1():
    pf = load_problem("case_ACTIVSg2000limit")
    res = solve(pf, pf.flat_start(), 1e-4, [0, 0.005, 1.0], Options(richardson=True))
    assert res.trajectory.times[-1] == 1.0 and res.converged
