from dataclasses import dataclass

import scipy.sparse as sp
import scipy.sparse.linalg as spla


class SingularMatrixError(RuntimeError):
    pass


@dataclass
class LUCounter:
    factorizations: int = 0


def factorize(a, counter=None):
    # Sparse LU of a. Returns a function that solves a x = b.
    if counter is not None:
        counter.factorizations += 1
    try:
        lu = spla.splu(sp.csc_matrix(a))
    except RuntimeError as exc:
        raise SingularMatrixError(str(exc)) from exc
    return lu.solve


def solve(a, b, counter=None):
    return factorize(a, counter)(b)
