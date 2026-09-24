"""Section 4.3 - Tables 2, 3, 4: FE, RK2 and BE with a large time step.

Pathway t = {0, 0.005, 0.01, 0.51, 1.0}: t1 = dt0 = 0.005 by BE, t2 = 2 dt0 by
the scheme under test, then dt = 0.5. K = 1e-4. NR refines x(1).
Each cell prints "ours / paper".
"""

from __future__ import annotations

from common import DT0_PAPER, write_table
from dynhomotopy.datasets import ALL
from dynhomotopy.homotopy import pathway
from norm_tables import case_label, comparison_rows, header, norm_row, status
from paper_values import TABLE2_FE, TABLE3_RK2, TABLE4_BE

TIMES = pathway.constant_step(DT0_PAPER, 0.5, t2=2 * DT0_PAPER)   # [0, .005, .01, .51, 1]
TABLES = [("table2_FE", "FE", TABLE2_FE), ("table3_RK2", "RK2", TABLE3_RK2), ("table4_BE", "BE", TABLE4_BE)]


def main(cases=ALL):
    for fname, method, paper in TABLES:
        rows = []
        for name in cases:
            ours, res = norm_row(name, TIMES, method)
            rows.append([case_label(name)] + comparison_rows(ours, paper.get(name)) + [status(res)])
        write_table(fname, ["case"] + header(TIMES) + ["NR converged (iters)"], rows,
                    f"{fname}: ||g(x)||_inf, {method} + NR, dt0 = {DT0_PAPER}, dt = 0.5  (ours / paper)")


if __name__ == "__main__":
    main()
