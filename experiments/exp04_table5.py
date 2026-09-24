"""Table 5: small first steps that grow along the path, followed by NR (K = 1e-4)."""

from dynhomotopy.homotopy import pathway

from .common import label, write_table
from .norm_tables import header, norm_row, ours_vs_paper, status
from .paper_values import TABLE5

TIMES = pathway.explicit(0, 0.005, 0.01, 0.02, 0.12, 0.62, 1.0)
CASES = ["case36964", "case109272", "case_ACTIVSg2000limit"]


def main(cases=CASES):
    rows = []
    for name in cases:
        for method in ["FE", "RK2", "BE"]:
            ours, res = norm_row(name, TIMES, method)
            paper = TABLE5.get((name, method))
            rows.append([label(name), method] + ours_vs_paper(ours, paper) + [status(res)])
    write_table("table5", ["case", "solver"] + header(TIMES) + ["NR converged (iterations)"], rows,
                "Table 5: growing time steps, K = 1e-4 (ours / paper)")


if __name__ == "__main__":
    main()
