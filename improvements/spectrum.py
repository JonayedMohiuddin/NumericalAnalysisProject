"""Eigenvalue estimates for the Jacobian: power method and (shifted) inverse iteration.

The first homotopy step solves with J + delta I (eq. 25), whose eigenvalues
are those of J moved by delta. These routines estimate the largest and the
smallest eigenvalue moduli of J, which say how ill-conditioned J is and how
large delta has to be to move the small eigenvalues away from zero.
"""

import numpy as np
import scipy.sparse as sp

from dynhomotopy.linalg import factorize


def power_method(apply, n, max_it=500, tol=1e-6, seed=0):
    """Largest eigenvalue modulus of the linear map `apply`.

    Iterating with the map applied twice, sqrt(||A A x|| / ||x||) converges to
    the largest modulus even when it belongs to a complex pair or a +-lambda
    pair, where the plain power method oscillates. Returns (modulus, iterations).
    """
    x = np.random.default_rng(seed).standard_normal(n)
    x /= np.linalg.norm(x)
    estimate = 0.0
    for it in range(1, max_it + 1):
        y = apply(apply(x))
        norm = np.linalg.norm(y)
        new_estimate = np.sqrt(norm)
        x = y / norm
        if abs(new_estimate - estimate) <= tol * new_estimate:
            return new_estimate, it
        estimate = new_estimate
    return estimate, max_it


def inverse_iteration(a, shift=0.0, max_it=200, tol=1e-6):
    """Distance from `shift` to the nearest eigenvalue of a (shifted inverse iteration).

    Runs the power method on (a - shift I)^-1, using one LU of a - shift I.
    With shift = 0 this is the smallest eigenvalue modulus of a.
    Returns (distance, iterations).
    """
    n = a.shape[0]
    solve = factorize(a - shift * sp.identity(n, format="csc"))
    modulus, it = power_method(solve, n, max_it, tol)
    return 1 / modulus, it


def spectrum_summary(jacobian, delta):
    """Eigenvalue moduli of J and J + delta I, and the resulting eigenvalue ratios."""
    n = jacobian.shape[0]
    lam_max, it_max = power_method(lambda v: jacobian @ v, n)
    lam_min, it_min = inverse_iteration(jacobian)
    shifted_min, it_shift = inverse_iteration(jacobian, shift=-delta)
    return {
        "lambda_max": lam_max,
        "lambda_min": lam_min,
        "shifted_min": shifted_min,
        "ratio_J": lam_max / lam_min,
        "ratio_shifted": (lam_max + delta) / shifted_min,
        "iterations": it_max + it_min + it_shift,
    }
