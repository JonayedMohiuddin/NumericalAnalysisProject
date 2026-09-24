"""Section 4.7 - Table 7: iterations of every method on all test systems.

NR-MAT : NR from the case file's own initial estimate (reference)
NR-flat: NR from a flat start
GSH-NR : static homotopy of ref. [12] (reconstruction, see dynhomotopy/solvers/gsh.py)
BE / RK2 (NR | FDXB): dynamic homotopy with K = 1e-4 (BE: t = {0, 0.005, 1};
         RK2: seven points), then the refining solver; iter_ref is reported.
Each cell prints "ours / paper".
"""

from __future__ import annotations

from common import label, write_table
from dynhomotopy.datasets import ALL, LIMIT, load_problem
from paper_values import TABLE7
from setups import GSH_PARAMS, methods

COLUMNS = ["NR-MAT", "NR-flat", "GSH-NR", "BE(NR)", "RK2(NR)", "BE(FDXB)", "RK2(FDXB)"]
PAPER_COL = {"NR-MAT": 0, "NR-flat": 1, "GSH-NR": 4, "BE(NR)": 5, "RK2(NR)": 6, "BE(FDXB)": 7, "RK2(FDXB)": 8}


def main(cases=ALL):
    rows = []
    for name in cases:
        pf = load_problem(name)
        paper = TABLE7[name]
        dh, delta = GSH_PARAMS[name]
        row = [label(name), dh, delta]
        for col, fn in methods(name).items():
            ok, it = fn(pf)
            row.append(f"{it if ok else 'fail'} / {paper[PAPER_COL[col]]}")
        rows.append(row)
        print("done", name)
    write_table("table7", ["case", "GSH dh1", "GSH delta"] + COLUMNS, rows,
                "Table 7: iterations (iter_NR or iter_ref), ours / paper")

    # NR-MAT and NR-flat for the two cases whose file stores baseMVA = 10
    native = []
    for name in [c for c in LIMIT if c in cases and c in ("case69limit", "case141limit")]:
        pf = load_problem(name, paper_base=False)
        m = methods(name)
        native.append([label(name)] + [
            (lambda r: f"{r[1] if r[0] else 'fail'} / {TABLE7[name][PAPER_COL[c]]}")(m[c](pf))
            for c in ["NR-MAT", "NR-flat", "BE(NR)"]])
    if native:
        write_table("table7_native_base", ["case (file base 10 MVA)", "NR-MAT", "NR-flat", "BE(NR)"], native,
                    "Table 7 check on the native 10 MVA base (ours / paper)")


if __name__ == "__main__":
    main()
