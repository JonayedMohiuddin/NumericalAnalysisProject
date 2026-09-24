"""Section 4.3 - Table 5: reduced initial steps with progressive increments.

Pathway t = {0, 0.005, 0.01, 0.02, 0.12, 0.62, 1.0}, K = 1e-4, for
case36964, case109272 and case_ACTIVSg2000limit, then 3 NR iterations.
"""

from __future__ import annotations

from common import write_table
from dynhomotopy.homotopy import pathway
from norm_tables import case_label, comparison_rows, header, norm_row, status
from paper_values import TABLE5

TIMES = pathway.explicit(0, 0.005, 0.01, 0.02, 0.12, 0.62, 1.0)
CASES = ["case36964", "case109272", "case_ACTIVSg2000limit"]


def main(cases=CASES):
    rows = []
    for name in cases:
        for method in ["FE", "RK2", "BE"]:
            ours, res = norm_row(name, TIMES, method)
            rows.append([case_label(name), method] + comparison_rows(ours, TABLE5.get((name, method)))
                        + [status(res)])
    write_table("table5", ["case", "solver"] + header(TIMES) + ["NR converged (iters)"], rows,
                "Table 5: ||g(x)||_inf with progressive time steps, K = 1e-4  (ours / paper)")


if __name__ == "__main__":
    main()
