"""Proposal: the spectrum of the Jacobian at the flat start.

For every case, the largest eigenvalue modulus of J (power method), the
smallest (inverse iteration) and the smallest of J + delta I for the paper's
delta = 0.02 (shifted inverse iteration at -delta). ARPACK (scipy eigs) is
used as a check. The ratio lambda_max / lambda_min shows how ill-conditioned
J is, and how much the shift in the first step (eq. 25) improves it.
"""

import numpy as np
import scipy.sparse.linalg as spla

from extensions.spectrum import spectrum_summary
from improvements.cases import ALL, load

from .common import write_table

DELTA = 0.02


def arpack(jacobian):
    small = spla.eigs(jacobian, k=1, sigma=0, return_eigenvectors=False)
    large = spla.eigs(jacobian, k=1, which="LM", return_eigenvectors=False, maxiter=5000)
    return float(np.abs(large[0])), float(np.abs(small[0]))


def main(cases=ALL):
    rows = []
    for name in cases:
        pf = load(name)
        jac = pf.jacobian(pf.flat_start()).tocsc()
        s = spectrum_summary(jac, DELTA)
        big, small = arpack(jac)
        rows.append([name, pf.n, f"{s['lambda_max']:.4g}", f"{big:.4g}", f"{s['lambda_min']:.4g}",
                     f"{small:.4g}", f"{s['shifted_min']:.4g}", f"{s['ratio_J']:.2e}",
                     f"{s['ratio_shifted']:.2e}"])
    write_table("spectrum",
                ["case", "n", "lambda_max (power)", "lambda_max (ARPACK)", "lambda_min (inverse it.)",
                 "lambda_min (ARPACK)", "min of J + 0.02 I", "ratio for J", "ratio for J + 0.02 I"],
                rows, "Eigenvalue moduli of the Jacobian at the flat start")


if __name__ == "__main__":
    main()
