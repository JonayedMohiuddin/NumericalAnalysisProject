"""Dense linear solvers written from scratch: Gauss elimination and LU with partial pivoting.

Both work on dense numpy arrays and cost O(n^3), so they are only practical
for the smaller test systems. The inner loops are vectorised over rows.
"""

import numpy as np


class SingularError(ValueError):
    pass


def gauss_solve(a, b):
    """Solve a x = b by Gauss elimination with partial pivoting and back substitution."""
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    for k in range(n - 1):
        p = k + np.argmax(np.abs(a[k:, k]))
        if a[p, k] == 0:
            raise SingularError(f"zero pivot in column {k}")
        if p != k:
            a[[k, p]] = a[[p, k]]
            b[[k, p]] = b[[p, k]]
        m = a[k + 1:, k] / a[k, k]
        a[k + 1:, k:] -= np.outer(m, a[k, k:])
        b[k + 1:] -= m * b[k]
    return back_substitution(a, b)


def back_substitution(u, y):
    x = np.zeros_like(y)
    for i in range(len(y) - 1, -1, -1):
        if u[i, i] == 0:
            raise SingularError(f"zero pivot in row {i}")
        x[i] = (y[i] - u[i, i + 1:] @ x[i + 1:]) / u[i, i]
    return x


def forward_substitution(l_unit, b):
    """Solve L y = b where L has a unit diagonal (stored below the diagonal)."""
    y = np.array(b, dtype=float)
    for i in range(1, len(y)):
        y[i] -= l_unit[i, :i] @ y[:i]
    return y


def lu_factor(a):
    """Doolittle LU with partial pivoting: P a = L U.

    Returns (lu, perm): L (unit diagonal, not stored) and U share the array
    `lu`, and perm[i] is the original row that ended up in row i.
    """
    lu = np.array(a, dtype=float)
    n = lu.shape[0]
    perm = np.arange(n)
    for k in range(n - 1):
        p = k + np.argmax(np.abs(lu[k:, k]))
        if lu[p, k] == 0:
            raise SingularError(f"zero pivot in column {k}")
        if p != k:
            lu[[k, p]] = lu[[p, k]]
            perm[[k, p]] = perm[[p, k]]
        lu[k + 1:, k] /= lu[k, k]
        lu[k + 1:, k + 1:] -= np.outer(lu[k + 1:, k], lu[k, k + 1:])
    if lu[-1, -1] == 0:
        raise SingularError("zero pivot in the last column")
    return lu, perm


def lu_solve(lu, perm, b):
    y = forward_substitution(lu, np.asarray(b, dtype=float)[perm])
    return back_substitution(lu, y)


def factorize(a):
    """Same interface as dynhomotopy.linalg.factorize: returns a function that solves a x = b."""
    lu, perm = lu_factor(a.toarray() if hasattr(a, "toarray") else a)
    return lambda b: lu_solve(lu, perm, b)
