import numpy as np
import scipy.sparse as sp

from dynhomotopy.datasets import load_problem
from extensions import scratch
from extensions.gauss import gauss_solve, lu_factor, lu_solve
from extensions.spectrum import inverse_iteration, power_method
from extensions.tuning import golden_section
from improvements import Options, solve


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
