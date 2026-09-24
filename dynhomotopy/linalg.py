"""Sparse LU factorisation with a factorisation counter.

The paper's cost argument is "one LU factorisation per pathway point"
(Section 3.2), so every factorisation in the package goes through here and is
counted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class SingularMatrixError(RuntimeError):
    """Raised when a matrix cannot be factorised."""


@dataclass
class LUCounter:
    factorizations: int = 0


def factorize(a: sp.spmatrix, counter: LUCounter | None = None):
    """Return a solve(b) callable for the sparse LU factors of ``a``."""
    if counter is not None:
        counter.factorizations += 1
    try:
        lu = spla.splu(sp.csc_matrix(a))
    except RuntimeError as exc:  # "Factor is exactly singular"
        raise SingularMatrixError(str(exc)) from exc
    return lu.solve


def solve(a: sp.spmatrix, b: np.ndarray, counter: LUCounter | None = None) -> np.ndarray:
    """Factorise ``a`` once and solve a x = b."""
    return factorize(a, counter)(b)
