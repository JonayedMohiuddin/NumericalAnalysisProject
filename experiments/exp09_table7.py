"""Table 7: iterations needed by each method on every test system.

NR-MAT starts from the initial guess stored in the case file, every other
method starts from a flat start. For BE and RK2 the table gives the
iterations of the refining solver (NR or FDXB) after the homotopy.
"""

from dynhomotopy.datasets import ALL, load_problem

from .common import label, write_table
from .paper_values import TABLE7
from .setups import GSH_PARAMS, methods

COLUMNS = ["NR-MAT", "NR-flat", "GSH-NR", "BE(NR)", "RK2(NR)", "BE(FDXB)", "RK2(FDXB)"]
PAPER_COLUMN = {"NR-MAT": 0, "NR-flat": 1, "GSH-NR": 4, "BE(NR)": 5, "RK2(NR)": 6,
                "BE(FDXB)": 7, "RK2(FDXB)": 8}


def cell(name, column, pf):
    converged, its = methods(name)[column](pf)
    paper = TABLE7[name][PAPER_COLUMN[column]]
    return f"{its if converged else 'fail'} / {paper}"


def main(cases=ALL):
    rows = []
    for name in cases:
        pf = load_problem(name)
        dh, delta = GSH_PARAMS[name]
        rows.append([label(name), dh, delta] + [cell(name, c, pf) for c in COLUMNS])
    write_table("table7", ["case", "GSH dh1", "GSH delta"] + COLUMNS, rows,
                "Table 7: iterations (ours / paper)")

    # The two small limit cases on the 10 MVA base stored in their files
    rows = []
    for name in ["case69limit", "case141limit"]:
        pf = load_problem(name, paper_base=False)
        rows.append([label(name)] + [cell(name, c, pf) for c in ["NR-MAT", "NR-flat", "BE(NR)"]])
    write_table("table7_native_base", ["case (10 MVA base)", "NR-MAT", "NR-flat", "BE(NR)"], rows,
                "Table 7 on the base MVA stored in the file (ours / paper)")


if __name__ == "__main__":
    main()
