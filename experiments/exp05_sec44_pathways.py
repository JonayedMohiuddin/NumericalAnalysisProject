"""Section 4.4: larger dt0 and K, keeping delta = K / dt0 = 0.02.

The paper tries dt0 = 0.05, K = 0.001 with paths of 3, 4 and 5 points and
reports that with few points only BE gives NR and FDXB a usable start.
It also reports that dt0 = 0.1, K = 0.002 does not work.
Each cell shows the NR / FDXB iterations after the homotopy, or "fail".
"""

from dynhomotopy.datasets import ALL, load_problem
from dynhomotopy.homotopy import pathway
from dynhomotopy.hybrid import solve_hybrid

from .common import label, write_table

RUNS = [
    (1e-3, pathway.explicit(0, 0.05, 1)),
    (1e-3, pathway.explicit(0, 0.05, 0.1, 1)),
    (1e-3, pathway.explicit(0, 0.05, 0.15, 1)),
    (1e-3, pathway.explicit(0, 0.05, 0.2, 1)),
    (1e-3, pathway.explicit(0, 0.05, 0.1, 0.2, 1)),
    (1e-3, pathway.explicit(0, 0.05, 0.1, 0.25, 1)),
    (1e-3, pathway.explicit(0, 0.05, 0.1, 0.3, 1)),
    (2e-3, pathway.explicit(0, 0.1, 1)),
    (2e-3, pathway.explicit(0, 0.1, 0.2, 0.3, 1)),
]
METHODS = ["FE", "RK2", "BE"]


def iterations(pf, K, times, method):
    cells = []
    for refiner in ["NR", "FDXB"]:
        r = solve_hybrid(pf, pf.flat_start(), K, times, method, refiner)
        cells.append(str(r.refine.iterations) if r.converged else "fail")
    return " / ".join(cells)


def main(cases=ALL):
    rows = []
    for K, times in RUNS:
        row = [K, "{" + ", ".join(f"{t:g}" for t in times) + "}"]
        for name in cases:
            pf = load_problem(name)
            row += [iterations(pf, K, times, m) for m in METHODS]
        rows.append(row)
    header = ["K", "path"] + [f"{label(c)} {m}" for c in cases for m in METHODS]
    write_table("sec44_pathways", header, rows, "Sec. 4.4: NR / FDXB iterations after the homotopy")


if __name__ == "__main__":
    main()
