"""Tables 2, 3 and 4: FE, RK2 and BE with a large time step, followed by NR.

Path: t = 0, 0.005 (BE), 0.01 (method under test), 0.51, 1.0 with K = 1e-4.
"""

from dynhomotopy.datasets import ALL
from dynhomotopy.homotopy import pathway

from .common import DT0_PAPER, label, write_table
from .norm_tables import header, norm_row, ours_vs_paper, status
from .paper_values import TABLE2_FE, TABLE3_RK2, TABLE4_BE

TIMES = pathway.constant_step(DT0_PAPER, 0.5, t2=2 * DT0_PAPER)
TABLES = [("table2_FE", "FE", TABLE2_FE), ("table3_RK2", "RK2", TABLE3_RK2), ("table4_BE", "BE", TABLE4_BE)]


def main(cases=ALL):
    for filename, method, paper in TABLES:
        rows = []
        for name in cases:
            ours, res = norm_row(name, TIMES, method)
            rows.append([label(name)] + ours_vs_paper(ours, paper.get(name)) + [status(res)])
        write_table(filename, ["case"] + header(TIMES) + ["NR converged (iterations)"], rows,
                    f"{method} + NR, dt0 = {DT0_PAPER}, dt = 0.5 (ours / paper)")


if __name__ == "__main__":
    main()
