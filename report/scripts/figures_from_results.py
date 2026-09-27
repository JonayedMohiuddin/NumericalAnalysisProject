"""Figures and LaTeX tables built only from the CSV/MD files in results/tables.

Run from the repository root:  python report/scripts/figures_from_results.py
"""

import csv
import re
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from style import (BLUE, CORAL, GOLD, GREEN, GREY, INK, LIGHT, MAROON, METHOD, NAVY2, PALE, PURPLE, SKY,
                   TABLES, save, write_tex)


def read_csv(name):
    with open(TABLES / f"{name}.csv", encoding="utf-8") as fh:
        return list(csv.reader(fh))


def read_md(name):
    rows = []
    for line in (TABLES / f"{name}.md").read_text(encoding="utf-8").splitlines():
        if line.startswith("|") and not line.startswith("|---"):
            rows.append([c.strip() for c in line.strip("|").split("|")])
    return rows[0], rows[1:]


def num(s):
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        return None


SHORT = {"case_ACTIVSg500limit": "case500limit", "case_ACTIVSg2000limit": "case2000limit"}


def short(name):
    return SHORT.get(name, name)


# --------------------------------------------------------------------------- first step
def fig_first_step():
    rows = read_csv("sec421_first_step")
    head, body = rows[0], rows[1:]
    cases = [r[0] for r in body]
    g0 = np.array([float(r[1]) for r in body])
    col = {h: i for i, h in enumerate(head)}
    series = [("FE", "FE K=0.0001: ||g(x(t1))||"), ("RK2", "RK2 K=0.0001: ||g(x(t1))||"),
              ("BE", "BE K=0.0001: ||g(x(t1))||")]
    fig, ax = plt.subplots(figsize=(7.2, 3.1))
    x = np.arange(len(cases))
    w = 0.2
    ax.bar(x - 1.5 * w, g0, w, color=GREY, label=r"$\|g(x^{(0)})\|_\infty$ (flat start)")
    for k, (m, c) in enumerate(series):
        vals = np.array([float(r[col[c]]) for r in body])
        ax.bar(x + (k - 0.5) * w, vals, w, color=METHOD[m], label=f"{m} step to $t_1$")
    ax.set_yscale("log")
    ax.set_xticks(x, cases, rotation=25, ha="right")
    ax.set_ylabel(r"$\|g(x(t_1))\|_\infty$  (log scale)")
    ax.set_title(r"One step to $t_1=0.005$ with $K=10^{-4}$: explicit steps explode, BE contracts")
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.33))
    save(fig, "first_step")
    ratios = [float(r[col["FE K=0.0001: ||g(x(t1))||"]]) / float(r[1]) for r in body]
    print("  FE blow-up factor range", min(ratios), max(ratios))


# --------------------------------------------------------------------------- fidelity
def paired_cells(name, skip):
    rows = read_csv(name)
    out = []
    for r in rows[1:]:
        for cell in r[skip:-1]:
            parts = cell.split(" / ")
            if len(parts) != 2:
                continue
            a, b = num(parts[0]), num(parts[1])
            if a is not None and b is not None and a > 0 and b > 0:
                out.append((a, b))
    return out


def fig_fidelity():
    sets = [("Table 2 (FE)", "table2_FE", 1, METHOD["FE"], "o"),
            ("Table 3 (RK2)", "table3_RK2", 1, METHOD["RK2"], "^"),
            ("Table 4 (BE)", "table4_BE", 1, METHOD["BE"], "s"),
            ("Table 5 (all)", "table5", 2, PURPLE, "D")]
    fig, ax = plt.subplots(figsize=(4.9, 4.4))
    lo, hi = 1e-13, 1e16
    ax.fill_between([lo, hi], [lo / 2, hi / 2], [lo * 2, hi * 2], color=LIGHT, lw=0,
                    label="within a factor of 2")
    ax.plot([lo, hi], [lo, hi], color=GREY, lw=0.8)
    stats = {}
    allpairs = []
    for label, name, skip, color, marker in sets:
        p = np.array(paired_cells(name, skip))
        allpairs.append(p)
        err = np.abs(np.log10(p[:, 0]) - np.log10(p[:, 1]))
        stats[label] = (len(p), int(np.sum(err <= np.log10(2))))
        ax.scatter(p[:, 1], p[:, 0], s=16, color=color, marker=marker, alpha=0.85, lw=0,
                   label=f"{label}: {stats[label][1]}/{stats[label][0]} cells")
    ax.set(xscale="log", yscale="log", xlim=(lo, hi), ylim=(lo, hi),
           xlabel=r"value printed in the paper, $\|g\|_\infty$", ylabel=r"our value, $\|g\|_\infty$",
           title="Reproduction fidelity, Tables 2-5")
    ax.legend(loc="upper left", fontsize=7.4)
    save(fig, "fidelity")
    p = np.vstack(allpairs)
    err = np.abs(np.log10(p[:, 0]) - np.log10(p[:, 1]))
    print("  fidelity", stats, "total", len(p), "within2x", int(np.sum(err <= np.log10(2))),
          "within10%", int(np.sum(err <= np.log10(1.1))))
    return stats, len(p), int(np.sum(err <= np.log10(2)))


# --------------------------------------------------------------------------- tables 2-5 as LaTeX
def fmt_cell(s):
    s = s.strip()
    if s in ("-", ""):
        return "--"
    if s == "diverged":
        return r"\textit{div.}"
    v = num(s)
    if v is None:
        return s
    if v == 0:
        return "0"
    a = abs(v)
    if a < 1e-2 or a >= 1e4:
        m, e = f"{v:.1e}".split("e")
        return rf"${m}\mathrm{{e}}{{{int(e)}}}$".replace("{", "{").replace(".0\\mathrm", "\\mathrm")
    return f"{v:g}"


def sci(v):
    if v is None:
        return "--"
    a = abs(v)
    if a == 0:
        return "0"
    if a < 1e-2 or a >= 1e4:
        e = int(np.floor(np.log10(a)))
        m = v / 10**e
        ms = f"{m:.1f}".rstrip("0").rstrip(".")
        return rf"${ms}{{\times}}10^{{{e}}}$" if ms != "1" else rf"$10^{{{e}}}$"
    if a >= 100:
        return f"{v:.0f}"
    return f"{v:.2g}"


def cell_pair(s):
    parts = s.split(" / ")
    if len(parts) != 2:
        return sci(num(s)) if num(s) is not None else "--"
    a, b = parts
    ta = r"\textit{div.}" if a.strip() == "diverged" else ("--" if a.strip() == "-" else sci(num(a)))
    tb = "--" if b.strip() == "-" else sci(num(b))
    return ta, tb


def table_norms(name, out, caption_cols):
    rows = read_csv(name)
    head, body = rows[0], rows[1:]
    lines = []
    for r in body:
        case = short(r[0])
        ours, paper = [], []
        for c in r[1:-1]:
            a, b = cell_pair(c)
            ours.append(a)
            paper.append(b)
        status = r[-1].replace("yes", "conv.").replace("no", "no")
        lines.append(rf"\texttt{{{case}}} & ours & " + " & ".join(ours) + rf" & {status} \\")
        lines.append(r" & paper & " + " & ".join(paper) + r" & \\ \addlinespace[2pt]")
    write_tex(out, "\n".join(lines) + "\n")


def table5_tex():
    rows = read_csv("table5")
    body = rows[1:]
    lines = []
    last = None
    for r in body:
        case, solver = short(r[0]), r[1]
        ours, paper = [], []
        for c in r[2:-1]:
            a, b = cell_pair(c)
            ours.append(a)
            paper.append(b)
        first = rf"\texttt{{{case}}}" if case != last else ""
        if case != last and last is not None:
            lines.append(r"\midrule")
        last = case
        lines.append(rf"{first} & {solver} & ours & " + " & ".join(ours) + r" \\")
        lines.append(r" & & paper & " + " & ".join(paper) + r" \\")
    write_tex("table5_rows", "\n".join(lines) + "\n")


# --------------------------------------------------------------------------- section 4.4 grid
def fig_sec44():
    head, body = read_md("sec44_pathways")  # the CSV has unquoted commas in the path column
    cols = head[2:]
    labels = [f"K={r[0]}  " + r[1] for r in body]
    code = np.zeros((len(body), len(cols)))
    text = [["" for _ in cols] for _ in body]
    for i, r in enumerate(body):
        for j, c in enumerate(r[2:]):
            nr, fd = [s.strip() for s in c.split("/")]
            nr_ok, fd_ok = nr != "fail", fd != "fail"
            code[i, j] = 3 if nr_ok and fd_ok else 2 if nr_ok else 1 if fd_ok else 0
            text[i][j] = nr if nr_ok else ""
    cmap = ListedColormap(["#D9534F", PALE, SKY, "#5AA469"])
    fig, ax = plt.subplots(figsize=(7.4, 3.5))
    ax.imshow(code, cmap=cmap, vmin=-0.5, vmax=3.5, aspect="auto")
    for i in range(len(body)):
        for j in range(len(cols)):
            if text[i][j]:
                ax.text(j, i, text[i][j], ha="center", va="center", fontsize=6.3, color="white")
    ax.set_yticks(range(len(body)), labels, fontsize=7)
    ax.set_xticks(range(len(cols)), [c.split()[-1] for c in cols], fontsize=6.5)
    ax.grid(False)
    cases = [short(c.split()[0]).replace("case", "").replace("limit", "lim") for c in cols[::3]]
    for k, cname in enumerate(cases):
        ax.text(3 * k + 1, -0.9, cname, ha="center", va="bottom", fontsize=7, fontweight="bold",
                color=INK)
        if k:
            ax.axvline(3 * k - 0.5, color="white", lw=2.2)
    ax.axhline(6.5, color="white", lw=2.2)
    ax.legend(handles=[Patch(color="#5AA469", label="NR and FDXB converge (number = NR its)"),
                       Patch(color=SKY, label="only NR converges"),
                       Patch(color=PALE, label="only FDXB converges"),
                       Patch(color="#D9534F", label="both fail")],
              ncol=2, loc="upper center", bbox_to_anchor=(0.45, -0.1), fontsize=7.2)
    ax.set_title("Section 4.4 re-run: which (path, integrator) pairs give a usable start", pad=22)
    save(fig, "sec44_grid")
    by_method = defaultdict(lambda: [0, 0])
    for i in range(len(body)):
        for j, c in enumerate(cols):
            m = c.split()[-1]
            by_method[m][0] += code[i, j] >= 2
            by_method[m][1] += 1
    print("  sec44 NR-success by method", dict(by_method))
    # few-point paths (rows 0 and 7 are 3-point paths)
    for i in (0, 7):
        print("  3-point path", labels[i], {m: sum(code[i, j] >= 2 for j, c in enumerate(cols)
                                                  if c.endswith(m)) for m in ("FE", "RK2", "BE")})


# --------------------------------------------------------------------------- timing
def parse_time(cell):
    ours, paper = [s.strip() for s in cell.split(" / ")]
    m = re.search(r"\((\d+)%\)", ours)
    o = float(m.group(1)) if m else None
    p = float(paper.rstrip("%")) if paper.endswith("%") else None
    return o, p


def fig_timing():
    rows = read_csv("table8")
    head, body = rows[0], rows[1:]
    cases = head[1:]
    methods = ["BE(NR)", "RK2(NR)", "BE(FDXB)", "RK2(FDXB)"]
    colors = {"BE(NR)": SKY, "RK2(NR)": PALE, "BE(FDXB)": MAROON, "RK2(FDXB)": NAVY2}
    data = {r[0]: [parse_time(c) for c in r[1:]] for r in body}
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    x = np.arange(len(cases))
    w = 0.2
    for k, m in enumerate(methods):
        o = [v[0] for v in data[m]]
        p = [v[1] for v in data[m]]
        pos = x + (k - 1.5) * w
        ax.bar(pos, o, w * 0.92, color=colors[m], label=f"{m} ours")
        ax.scatter(pos, p, marker="_", s=110, color=INK, lw=1.8, zorder=4,
                   label="paper (MATLAB)" if k == 0 else None)
    ax.axhline(100, color=GREY, ls="--", lw=0.9)
    ax.text(len(cases) - 0.45, 104, "NR-MAT = 100%", fontsize=7, color=GREY, ha="right")
    ax.set_xticks(x, cases)
    ax.set_ylabel("run time, % of NR-MAT")
    ax.set_title("Table 8 re-run: cost relative to Newton from the case-file start")
    ax.legend(ncol=5, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    save(fig, "timing")


# --------------------------------------------------------------------------- improvements
# Order used in every table and figure: the paper, the corrector / multiplier / adaptive family,
# Richardson step control, then the ideas that did not work.
ORDER = ["paper", "PC", "OM + PC", "OM + PC + adaptive", "PC + adaptive", "OM", "OM + adaptive",
         "adaptive", "richardson", "OM + richardson", "richardson + FE", "richardson + RK2",
         "richardson + RK4", "scaled homotopy", "BE-chord", "newton homotopy", "RK4"]
FAMILY = {**{o: "helpful" for o in ORDER[1:8]}, **{o: "richardson" for o in ORDER[8:13]},
          **{o: "failed" for o in ORDER[13:]}, "paper": "paper"}
DISPLAY = {"richardson": "richardson (BE)", "RK4": "RK4 (fixed steps)"}
SETTINGS = ["S1 paper", "S2 Sec 4.4", "S3 Sec 4.4", "S4 weak K", "S5 strong K"]
SET_LABEL = {"S1 paper": "S1", "S2 Sec 4.4": "S2", "S3 Sec 4.4": "S3", "S4 weak K": "S4",
             "S5 strong K": "S5"}
SOLVED, OTHER = "solved", "other root"


def disp(o):
    return DISPLAY.get(o, o)


def load_raw():
    res, cases = {}, []
    with open(TABLES / "improvements_raw.csv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            res[(r["case"], r["setting"], r["options"])] = (r["outcome"], int(r["lu"]))
            if r["case"] not in cases:
                cases.append(r["case"])
    return res, cases


def improvement_stats(res, cases):
    def ok(c, s, o):
        return res[(c, s, o)][0] == SOLVED

    stats = {}
    for o in ORDER:
        per = [sum(ok(c, s, o) for c in cases) for s in SETTINGS]
        other = sum(res[(c, s, o)][0] == OTHER for c in cases for s in SETTINGS)
        gained = [(c, s) for c in cases for s in SETTINGS if ok(c, s, o) and not ok(c, s, "paper")]
        lost = [(c, s) for c in cases for s in SETTINGS if not ok(c, s, o) and ok(c, s, "paper")]
        common = [(c, s) for c in cases for s in SETTINGS if ok(c, s, o) and ok(c, s, "paper")]
        lu_o = sum(res[(c, s, o)][1] for c, s in common)
        lu_p = sum(res[(c, s, "paper")][1] for c, s in common)
        common1 = [c for c, s in common if s == "S1 paper"]
        stats[o] = dict(solved=sum(per), per=per, other=other, gained=gained, lost=lost,
                        lu_ratio=lu_o / lu_p if lu_p else np.nan, common=len(common),
                        s1=np.mean([res[(c, "S1 paper", o)][1] for c in common1]),
                        s1_paper=np.mean([res[(c, "S1 paper", "paper")][1] for c in common1]),
                        s1_n=len(common1))
    return stats


def gains_kinds():
    """(case, setting) -> 'rescue' (paper diverges), 'root' (paper reaches another root), 'speed' (slow)."""
    head, body = read_md("improvements_gains")
    kind = {}
    for r in body:
        kind[(r[1], r[0])] = ("rescue" if r[2].startswith("diverges") else
                              "root" if "another root" in r[2] else "speed")
    return kind


def split(gained, kind):
    return tuple(sum(kind.get(k) == t for k in gained) for t in ("rescue", "root", "speed"))


def fig_improvements():
    res, cases = load_raw()
    st = improvement_stats(res, cases)
    kind = gains_kinds()

    # ---- solved per configuration, stacked by setting, with other-root runs hatched
    fig, ax = plt.subplots(figsize=(7.3, 5.0))
    set_colors = [MAROON, "#2E5E8C", SKY, "#7DBBE6", "#BFE0F5"]
    y = np.arange(len(ORDER))[::-1]
    left = np.zeros(len(ORDER))
    for k, s in enumerate(SETTINGS):
        vals = np.array([st[o]["per"][k] for o in ORDER])
        ax.barh(y, vals, left=left, color=set_colors[k], height=0.66, label=SET_LABEL[s],
                edgecolor="white", lw=0.6)
        left += vals
    others = np.array([st[o]["other"] for o in ORDER])
    ax.barh(y, others, left=left, color="#FBE3B5", height=0.66, edgecolor=CORAL, lw=0.5, hatch="////",
            label="converged to another root")
    for yi, o in zip(y, ORDER):
        txt = f"{st[o]['solved']}"
        if o != "paper":
            txt += f"   +{len(st[o]['gained'])} / -{len(st[o]['lost'])}"
        ax.text(st[o]["solved"] + st[o]["other"] + 2, yi, txt, va="center", fontsize=7.4, color=INK)
    ax.axvline(st["paper"]["solved"], color=MAROON, ls="--", lw=1)
    ax.text(st["paper"]["solved"] - 1.5, len(ORDER) - 0.35, f"paper's method: {st['paper']['solved']}",
            color=MAROON, fontsize=7.4, ha="right")
    for b in (7.5, 12.5):
        ax.axhline(len(ORDER) - 1 - b + 0.0, color=GREY, lw=0.6, ls=":")
    ax.set_yticks(y, [disp(o) for o in ORDER])
    ax.set_xlim(0, 205)
    ax.set_xlabel("runs that reach the reference solution, out of 170 (34 cases x 5 settings)")
    ax.set_title("Every configuration on the 34-case test bed  (+gained / -lost vs. paper)", pad=10)
    ax.legend(ncol=6, loc="upper center", bbox_to_anchor=(0.42, -0.1), fontsize=7.2,
              handlelength=1.4, columnspacing=1.0)
    ax.grid(axis="y", visible=False)
    save(fig, "improve_solved")

    # ---- pareto: LU ratio vs solved (log x)
    fig, ax = plt.subplots(figsize=(6.2, 4.1))
    colors = {"paper": MAROON, "helpful": SKY, "richardson": CORAL, "failed": GREY}
    place = {
        "paper": (0.84, 118, "center"),
        "PC": (1.0, 178, "center"), "OM + PC": (0.8, 165, "center"),
        "OM + PC + adaptive": (1.5, 172, "left"), "PC + adaptive": (1.62, 150, "left"),
        "OM": (0.8, 136, "center"), "OM + adaptive": (0.79, 150, "center"),
        "adaptive": (1.22, 133, "left"),
        "richardson": (1.62, 141, "left"), "OM + richardson": (1.62, 161, "left"),
        "richardson + FE": (3.1, 168, "center"), "richardson + RK2": (6.2, 166, "center"),
        "richardson + RK4": (7.0, 136, "center"),
        "scaled homotopy": (1.22, 124, "left"), "BE-chord": (1.2, 92, "left"),
        "newton homotopy": (1.2, 78, "left"),
    }
    for o in ORDER:
        if o == "RK4":
            continue
        x, yv = st[o]["lu_ratio"], st[o]["solved"]
        ax.scatter(x, yv, s=46, color=colors[FAMILY[o]], zorder=3, edgecolor="white", lw=0.7)
        lx, ly, ha = place[o]
        ax.annotate(disp(o), (x, yv), (lx, ly), fontsize=7, ha=ha, va="center", color=INK,
                    arrowprops=dict(arrowstyle="-", color="#b7b3ad", lw=0.6, shrinkA=0, shrinkB=3))
    ax.text(11.5, 104, f"off-scale: RK4 (fixed steps) at ({st['RK4']['lu_ratio']:.2f}, {st['RK4']['solved']})",
            fontsize=6.8, color=GREY, ha="right", va="bottom")
    ax.set_xscale("log")
    ax.set_xticks([0.8, 1, 1.5, 2, 3, 5, 10], ["0.8", "1", "1.5", "2", "3", "5", "10"])
    ax.minorticks_off()
    ax.set_xlabel("LU factorizations relative to the paper's method (log scale)\n"
                  "(runs both solve, all five settings pooled)")
    ax.set_ylabel("runs solved (of 170)")
    ax.set_title("Robustness against cost: up and to the left is better")
    ax.axhline(st["paper"]["solved"], color=MAROON, lw=0.7, ls=":")
    ax.axvline(1.0, color=MAROON, lw=0.7, ls=":")
    ax.set_xlim(0.72, 12)
    ax.set_ylim(70, 182)
    ax.legend(handles=[Patch(color=SKY, label="corrector / multiplier / adaptive"),
                       Patch(color=CORAL, label="Richardson step control"),
                       Patch(color=GREY, label="did not work")],
              loc="lower right", fontsize=7)
    save(fig, "improve_pareto")

    # ---- per-case S1 heatmap
    fig, ax = plt.subplots(figsize=(7.3, 8.6))
    cmap = plt.get_cmap("Blues")
    norm = LogNorm(3, 150)
    rgba = np.zeros((len(cases), len(ORDER), 4))
    for i, c in enumerate(cases):
        for j, o in enumerate(ORDER):
            outcome, lu = res[(c, "S1 paper", o)]
            if outcome == SOLVED:
                rgba[i, j] = cmap(0.08 + 0.85 * norm(lu))
            else:
                rgba[i, j] = plt.matplotlib.colors.to_rgba("#FBE3B5" if outcome == OTHER else "#F4D6D2")
            if outcome == SOLVED:
                p_out, p_lu = res[(c, "S1 paper", "paper")]
                bold = p_out != SOLVED or lu < p_lu
                ax.text(j, i, f"{lu}", ha="center", va="center", fontsize=6.1,
                        color="white" if lu > 12 else INK, fontweight="bold" if bold and j else "normal")
            else:
                ax.text(j, i, "o" if outcome == OTHER else "x", ha="center", va="center", fontsize=6.6,
                        color="#A0522D" if outcome == OTHER else "#B03A2E")
    ax.imshow(rgba, aspect="auto")
    ax.set_xticks(range(len(ORDER)), [disp(o) for o in ORDER], fontsize=7)
    ax.set_yticks(range(len(cases)), [short(c) for c in cases], fontsize=7)
    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=45, ha="left")
    ax.grid(False)
    for yy in (8.5, 15.5):
        ax.axhline(yy, color="white", lw=2.5)
    for xx in (0.5, 7.5, 12.5):
        ax.axvline(xx, color="white", lw=2.5)
    for yy, lab in ((4, "paper's 9 cases"), (12, "7 extra ill-conditioned"), (25, "18 extra limit cases")):
        ax.text(len(ORDER) - 0.3, yy, lab, rotation=-90, va="center", ha="left", fontsize=7.4, color=GREY)
    save(fig, "improve_cases_s1")

    # ---- LaTeX table
    lines = []
    for o in ORDER:
        s = st[o]
        per = " & ".join(str(v) for v in s["per"])
        name = r"\textbf{paper}" if o == "paper" else disp(o)
        if o == "paper":
            gl = "-- & --"
            lu = f"{s['s1']:.2f}"
        else:
            r_, o_, sp_ = split(s["gained"], kind)
            gl = f"+{len(s['gained'])} ({r_}/{o_}/{sp_}) & $-${len(s['lost'])}"
            lu = f"{s['s1']:.2f}" + ("" if abs(s["s1_paper"] - st["paper"]["s1"]) < 1e-9
                                      else rf" \textcolor{{ink!55}}{{({s['s1_paper']:.2f})}}")
        lines.append(rf"{name} & {per} & \textbf{{{s['solved']}}} & {s['other']} & {gl} & {lu} & "
                     rf"{s['lu_ratio']:.2f} \\")
        if o in ("paper", "adaptive", "richardson + RK4"):
            lines.append(r"\midrule")
    write_tex("improve_table", "\n".join(lines) + "\n")

    for o in ORDER:
        s = st[o]
        print(f"  {o:20s} solved {s['solved']:3d} other {s['other']:2d} gained {len(s['gained']):2d} "
              f"{split(s['gained'], kind)} lost {len(s['lost']):3d} lu_ratio {s['lu_ratio']:.3f} "
              f"S1 {s['s1']:.2f} vs {s['s1_paper']:.2f} ({s['s1_n']})")
        if s["lost"] and len(s["lost"]) < 12:
            print("     lost:", s["lost"])
    never = [(c, s) for c in cases for s in SETTINGS
             if not any(res[(c, s, o)][0] == SOLVED for o in ORDER)]
    print("  never solved:", never)
    paper_other = [(c, s) for c in cases for s in SETTINGS if res[(c, s, "paper")][0] == OTHER]
    print("  paper other root:", paper_other)
    return st


# --------------------------------------------------------------------------- proposal: spectrum
ILL = {"case18482", "case27318", "case36964", "case54636", "case109272", "case6024", "case6243",
       "case6748", "case7092", "case9961", "case10595", "case12110"}


def fig_spectrum():
    rows = read_csv("spectrum")
    body = rows[1:]
    names = [short(r[0]) for r in body]
    lam_max_p, lam_max_a = [float(r[2]) for r in body], [float(r[3]) for r in body]
    lam_min, lam_min_a, shifted = [float(r[4]) for r in body], [float(r[5]) for r in body], [float(r[6]) for r in body]
    ratio, ratio_s = [float(r[7]) for r in body], [float(r[8]) for r in body]
    order = sorted(range(len(body)), key=lambda i: (body[i][0] not in ILL, -lam_min[i]))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.4, 5.4), sharey=True,
                                  gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.07))
    y = np.arange(len(order))[::-1]
    for yi, i in zip(y, order):
        c = MAROON if body[i][0] in ILL else SKY
        ax.plot([lam_min[i], shifted[i]], [yi, yi], color=c, lw=1.2, alpha=0.6)
        ax.scatter(lam_min[i], yi, color=c, s=18, zorder=3)
        ax.scatter(shifted[i], yi, facecolor="white", edgecolor=c, s=18, zorder=3, lw=1.1)
        ax2.plot([ratio_s[i], ratio[i]], [yi, yi], color=c, lw=1.2, alpha=0.6)
        ax2.scatter(ratio[i], yi, color=c, s=18, zorder=3)
        ax2.scatter(ratio_s[i], yi, facecolor="white", edgecolor=c, s=18, zorder=3, lw=1.1)
    ax.axvline(0.02, color=CORAL, ls="--", lw=0.9)
    ax.text(0.022, 13.5, r"$\delta=0.02$", color=CORAL, fontsize=7.5, va="bottom")
    ax.set_xscale("log")
    ax.set_yticks(y, [names[i] for i in order], fontsize=6.8)
    ax.set_xlabel(r"smallest eigenvalue modulus $|\lambda|_{\min}$")
    ax.set_title(r"$|\lambda|_{\min}$: $J(x^{(0)})$ filled, $J+0.02I$ open", fontsize=9)
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$|\lambda|_{\max}/|\lambda|_{\min}$")
    ax2.set_title("Eigenvalue ratio, before and after", fontsize=9)
    for a in (ax, ax2):
        a.grid(axis="y", visible=False)
        a.axhline(len(order) - 12.5, color=GREY, lw=0.6, ls=":")
    fig.legend(handles=[Line2D([], [], color=MAROON, marker="o", ls="", label="ill-conditioned (12)"),
                        Line2D([], [], color=SKY, marker="o", ls="", label="loading limit (22)")],
               ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.03), fontsize=7.4)
    save(fig, "spectrum")
    ill = [i for i in range(len(body)) if body[i][0] in ILL]
    lim = [i for i in range(len(body)) if body[i][0] not in ILL]
    err_max = [abs(lam_max_p[i] - lam_max_a[i]) / lam_max_a[i] for i in range(len(body))]
    err_min = [abs(lam_min[i] - lam_min_a[i]) / lam_min_a[i] for i in range(len(body))]
    print("  spectrum ill lam_min", min(lam_min[i] for i in ill), max(lam_min[i] for i in ill))
    print("  spectrum lim lam_min", sorted(round(lam_min[i], 4) for i in lim))
    print("  spectrum ill ratio_s", min(ratio_s[i] for i in ill), max(ratio_s[i] for i in ill))
    print("  power method max rel err", max(err_max), "cases >1%:",
          [(body[i][0], round(err_max[i], 3)) for i in range(len(body)) if err_max[i] > 0.01])
    print("  inverse iteration cases >1e-4:", [(body[i][0], round(err_min[i], 4))
                                                for i in range(len(body)) if err_min[i] > 1e-4])


# --------------------------------------------------------------------------- proposal: tuning
def fig_tuning():
    rows = read_csv("tuning_per_case")
    body = rows[1:]
    paper = {"case18482", "case27318", "case36964", "case54636", "case109272", "case69limit",
             "case141limit", "case_ACTIVSg500limit", "case_ACTIVSg2000limit"}
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    for r in body:
        lam, delta, ok = float(r[1]), float(r[2]), r[5] == "yes"
        c = MAROON if r[0] in ILL else SKY
        ax.scatter(lam, delta, s=34 if r[0] in paper else 22, zorder=3,
                   marker="o" if ok else "X", color=c if r[0] in paper else "white",
                   edgecolor=c, lw=1.1)
    xs = np.logspace(-3.6, 0.5, 20)
    ax.plot(xs, 1.86 * xs, color=CORAL, lw=1, ls="--")
    ax.text(0.5, 1.86 * 0.5 * 1.5, r"spectral rule $\delta=1.86\,|\lambda|_{\min}$", color=CORAL,
            fontsize=7, ha="right", va="bottom")
    ax.axhline(0.02, color=GREEN, lw=1, ls=":")
    ax.text(3.2, 0.024, r"paper $\delta=0.02$", color=GREEN, fontsize=7, ha="right", va="bottom")
    ax.set(xscale="log", yscale="log", xlabel=r"$|\lambda|_{\min}$ of $J(x^{(0)})$",
           ylabel=r"tuned $\delta = K/\Delta t_0$", title="Golden-section tuned shift against the spectrum")
    ax.set_ylim(5e-5, 30)
    ax.legend(handles=[Line2D([], [], color=MAROON, marker="o", ls="", label="ill-conditioned"),
                       Line2D([], [], color=SKY, marker="o", ls="", label="loading limit"),
                       Line2D([], [], color=INK, marker="o", ls="", label="filled: paper's 9 (training)"),
                       Line2D([], [], color=INK, marker="o", mfc="white", ls="", label="open: 25 others (test)"),
                       Line2D([], [], color=INK, marker="X", ls="", label="not solved at tuned point")],
              fontsize=6.8, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3)
    save(fig, "tuning")
    ill_ratio = [float(r[2]) / float(r[1]) for r in body if r[0] in ILL]
    print("  tuned delta / lam_min on ill cases", min(ill_ratio), max(ill_ratio))


# --------------------------------------------------------------------------- figure 5 trade-off
def fig_fig5():
    rows = read_csv("fig5_timestep")
    body = rows[1:]
    dt = [r[0] for r in body]
    hom = np.array([int(r[2]) for r in body])
    nr = np.array([int(r[5]) for r in body])
    ang = [float(r[4]) for r in body]
    fig, ax = plt.subplots(figsize=(4.8, 2.9))
    x = np.arange(len(dt))
    ax.bar(x, hom, 0.55, color=SKY, label="homotopy LUs")
    ax.bar(x, nr, 0.55, bottom=hom, color=MAROON, label="NR LUs")
    for xi, h, n, a in zip(x, hom, nr, ang):
        ax.text(xi, h + n + 0.3, f"{h + n}", ha="center", fontsize=8, color=INK)
    ax.set_xticks(x, [rf"$\Delta t$={d}" for d in dt])
    ax.set_ylabel("LU factorizations")
    ax.set_ylim(0, 16.5)
    ax.set_title("case109272: smaller steps buy nothing")
    ax.legend(loc="upper left")
    save(fig, "fig5_tradeoff")


if __name__ == "__main__":
    fig_first_step()
    fig_fidelity()
    table_norms("table2_FE", "table2_rows", None)
    table_norms("table3_RK2", "table3_rows", None)
    table_norms("table4_BE", "table4_rows", None)
    table5_tex()
    fig_sec44()
    fig_timing()
    fig_improvements()
    fig_spectrum()
    fig_tuning()
    fig_fig5()
