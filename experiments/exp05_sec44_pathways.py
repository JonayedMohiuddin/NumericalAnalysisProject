"""Section 4.4 - larger dt0 and K with the same perturbation delta = K/dt0 = 0.02.

dt0 = 0.05, K = 0.001; x(t1) always by BE. The pathways tested in the paper:
  3 points : {0, 0.05, 1}
  4 points : {0, 0.05, 0.1, 1}, {0, 0.05, 0.15, 1}, {0, 0.05, 0.2, 1}
  5 points : {0, 0.05, 0.1, t3, 1} with t3 in {0.20, 0.25, 0.30}
Paper's finding: with 3-4 points only BE yields an x^(0) from which NR and
FDXB converge; with 5 points all schemes do (except FE for t3 = 0.25 on
case36964). Last row: dt0 = 0.1, K = 0.002 (delta still 0.02) fails.
Each cell: NR / FDXB iteration count, or "fail".
"""

from __future__ import annotations

from common import label, write_table
from dynhomotopy.datasets import ALL, load_problem
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid

K_44 = 1e-3
PATHWAYS = [
    pathway.explicit(0, 0.05, 1.0),
    pathway.explicit(0, 0.05, 0.1, 1.0),
    pathway.explicit(0, 0.05, 0.15, 1.0),
    pathway.explicit(0, 0.05, 0.2, 1.0),
    pathway.explicit(0, 0.05, 0.1, 0.2, 1.0),
    pathway.explicit(0, 0.05, 0.1, 0.25, 1.0),
    pathway.explicit(0, 0.05, 0.1, 0.3, 1.0),
]
EXTRA = [(2e-3, pathway.explicit(0, 0.1, 1.0)), (2e-3, pathway.explicit(0, 0.1, 0.2, 0.3, 1.0))]


def cell(pf, K, times, method):
    out = []
    for refiner in ["NR", "FDXB"]:
        r = solve_hybrid(pf, pf.flat_start(), K, times, method, refiner)
        out.append(str(r.refine.iterations) if r.converged else "fail")
    return " / ".join(out)


def main(cases=ALL):
    configs = [(K_44, t) for t in PATHWAYS] + EXTRA
    header = ["K", "pathway"] + [f"{label(c)} {m}" for c in cases for m in ["FE", "RK2", "BE"]]
    rows = []
    for K, times in configs:
        row = [K, "{" + ", ".join(f"{t:g}" for t in times) + "}"]
        for name in cases:
            pf = load_problem(name)
            row += [cell(pf, K, times, m) for m in ["FE", "RK2", "BE"]]
        rows.append(row)
        print("done", row[1])
    write_table("sec44_pathways", header, rows,
                "Sec. 4.4: refiner iterations NR / FDXB from x(1) of each scheme (delta = 0.02)")


if __name__ == "__main__":
    main()
