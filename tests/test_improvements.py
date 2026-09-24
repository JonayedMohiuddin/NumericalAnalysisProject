import numpy as np
import pytest

from dynhomotopy.datasets import load_problem
from dynhomotopy.homotopy.fpv import FixedPointHomotopy
from dynhomotopy.homotopy.integrators import backward_euler
from dynhomotopy.hybrid import solve_hybrid
from improvements import Options, solve
from improvements.cases import load
from improvements.corrector import backward_euler_corrected
from improvements.multiplier import optimal_multiplier

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
    corrected = backward_euler_corrected(h, h.x0, 0.0, 0.005)
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
