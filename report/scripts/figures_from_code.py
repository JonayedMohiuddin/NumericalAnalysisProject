"""Figures that call the project's own code (dynhomotopy / improvements).

Needs the .mat test systems in data/cases (downloaded on first use).
Run from anywhere:  python report/scripts/figures_from_code.py
"""

import json
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from style import (BLUE, CORAL, GEN, GOLD, GREEN, GREY, INK, LIGHT, MAROON, METHOD, NAVY2, PALE, PURPLE, SKY,
                   REPO, save, write_tex)

sys.path.insert(0, str(REPO))

from dynhomotopy import datasets  # noqa: E402
from dynhomotopy.homotopy import pathway  # noqa: E402
from dynhomotopy.homotopy.fpv import FixedPointHomotopy  # noqa: E402
from dynhomotopy.homotopy.integrators import backward_euler  # noqa: E402
from dynhomotopy.hybrid import solve_hybrid  # noqa: E402
from dynhomotopy.problem import TutorialProblem  # noqa: E402
from dynhomotopy.solvers import newton_raphson  # noqa: E402
from improvements import Options, solve  # noqa: E402
from improvements import path as imp_path  # noqa: E402
from improvements.cases import ALL, EXTRA_ILL, load, on_reference  # noqa: E402
from improvements.steps import with_corrector  # noqa: E402

SHORT = {"case_ACTIVSg500limit": "case500limit", "case_ACTIVSg2000limit": "case2000limit"}
BUILT = {
    "case18482": "9241pegase x2", "case27318": "13659pegase x2", "case36964": "9241pegase x4",
    "case54636": "13659pegase x4", "case109272": "13659pegase x8",
    "case6024": "3012wp x2", "case6243": "3375wp + 2869pegase", "case6748": "3375wp x2",
    "case7092": "2869pegase x2 + 1354pegase", "case9961": "2869pegase x3 + 1354pegase",
    "case10595": "9241pegase + 1354pegase", "case12110": "9241pegase + 2869pegase",
}


def group(name):
    if name in datasets.SYSTEMS:
        return "paper-ill" if datasets.SYSTEMS[name].kind == "ill" else "paper-limit"
    return "extra-ill" if name in EXTRA_ILL else "extra-limit"


# --------------------------------------------------------------------------- dataset statistics
def systems():
    cache = GEN / "systems.json"
    if cache.exists():
        return json.loads(cache.read_text())
    out = []
    for name in ALL:
        pf = load(name)
        x0 = pf.flat_start()
        J = pf.jacobian(x0)
        nr = newton_raphson(pf, x0)
        mat = newton_raphson(pf, pf.case_start())
        out.append(dict(name=name, group=group(name), nb=int(pf.case.nb), npv=int(pf.npv),
                        npq=int(pf.npq), n=int(pf.n), nnz=int(J.nnz),
                        nbranch=int(pf.case.branch.shape[0]), g0=float(pf.norm(x0)),
                        flat_ok=bool(nr.converged), flat_its=int(nr.iterations),
                        mat_ok=bool(mat.converged), mat_its=int(mat.iterations)))
        print("  stats", name, out[-1]["nb"], out[-1]["flat_ok"], flush=True)
    cache.write_text(json.dumps(out, indent=1))
    return out


def table_systems(stats):
    lines, last = [], None
    heads = {"paper-ill": "Paper, ill-conditioned (record 3514739)",
             "paper-limit": "Paper, loading limit (record 3491654)",
             "extra-ill": "Added by us, ill-conditioned (record 3514739)",
             "extra-limit": "Added by us, loading limit (record 3491654)"}
    for s in stats:
        if s["group"] != last:
            if last is not None:
                lines.append(r"\midrule")
            lines.append(rf"\multicolumn{{8}}{{l}}{{\textit{{{heads[s['group']]}}}}} \\")
            last = s["group"]
        name = SHORT.get(s["name"], s["name"]).replace("_", r"\_")
        built = BUILT.get(s["name"], "")
        flat = f"{s['flat_its']}" if s["flat_ok"] else r"\textcolor{failred}{fail}"
        g0 = f"{s['g0']:.3g}" if s["g0"] >= 0.01 else f"{s['g0']:.1e}"
        lines.append(rf"\texttt{{{name}}} & {built} & {s['nb']:,} & {s['npv']:,} & {s['n']:,} & "
                     rf"{s['nnz']:,} & {g0} & {flat} \\".replace(",", "{,}"))
    write_tex("systems_rows", "\n".join(lines) + "\n")


def fig_systems(stats):
    colors = {"paper-ill": MAROON, "paper-limit": SKY, "extra-ill": NAVY2, "extra-limit": PALE}
    stats = sorted(stats, key=lambda s: s["nb"])
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.4, 4.9), sharey=True,
                                  gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.06))
    y = np.arange(len(stats))
    for yi, s in zip(y, stats):
        ax.barh(yi, s["nb"], color=colors[s["group"]], height=0.72)
        ax2.scatter(s["g0"], yi, color=colors[s["group"]], s=22, zorder=3,
                    marker="o" if s["flat_ok"] else "X")
    ax.set_xscale("log")
    ax.set_yticks(y, [SHORT.get(s["name"], s["name"]) for s in stats], fontsize=6.8)
    ax.set_xlabel("buses (log scale)")
    ax.set_title("Size of the 34 test systems")
    ax.grid(axis="y", visible=False)
    ax2.set_xscale("log")
    ax2.set_xlabel(r"flat-start mismatch $\|g(x^{(0)})\|_\infty$")
    ax2.set_title("How far the flat start is")
    ax2.grid(axis="y", visible=False)
    handles = [Line2D([], [], color=c, lw=6, label=l) for l, c in
               [("paper, ill-conditioned", MAROON), ("paper, loading limit", SKY),
                ("extra, ill-conditioned", NAVY2), ("extra, loading limit", PALE)]]
    handles += [Line2D([], [], color=INK, marker="o", ls="", label="NR from flat start converges"),
                Line2D([], [], color=INK, marker="X", ls="", label="NR from flat start fails")]
    fig.legend(handles=handles, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.06), fontsize=7.2)
    save(fig, "systems")


# --------------------------------------------------------------------------- sparsity
def fig_spy():
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.75))
    for ax, name in zip(axes, ["case_ACTIVSg2000limit", "case18482", "case109272"]):
        pf = load(name)
        J = pf.jacobian(pf.flat_start()).tocoo()
        ax.scatter(J.col, J.row, s=0.02 if J.nnz > 2e5 else 0.08, color=MAROON, lw=0,
                   rasterized=True)
        ax.set_xlim(0, pf.n)
        ax.set_ylim(pf.n, 0)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_color(GREY)
        nv = pf.npv + pf.npq
        ax.axhline(nv, color=SKY, lw=0.7)
        ax.axvline(nv, color=SKY, lw=0.7)
        ax.set_title(SHORT.get(name, name), fontsize=9)
        ax.set_xlabel(f"n = {pf.n:,}   nnz = {J.nnz:,}", fontsize=7.5)
    save(fig, "spy")


# --------------------------------------------------------------------------- tutorial geometry
def trace_true_path(problem, K, steps=800):
    h = FixedPointHomotopy(problem, problem.x0, K)
    step = with_corrector(backward_euler)
    ts = np.linspace(0, 1, steps + 1)
    ts = np.concatenate([[0], 0.01 * np.linspace(0, 1, 60)[1:], np.linspace(0.01, 1, steps)[1:]])
    x = problem.x0.astype(float)
    xs = [x]
    for t, t2 in zip(ts, ts[1:]):
        x = step(h, x, t, t2 - t)
        for _ in range(3):  # extra corrector sweeps keep it on the curve
            from dynhomotopy.linalg import solve as lsolve
            x = x - lsolve(h.Gx(x, t2), h.G(x, t2))
        xs.append(x)
    return np.array(xs)


def fig_tutorial():
    problem = TutorialProblem()
    g = np.linspace(-0.2, 2.2, 400)
    X1, X2 = np.meshgrid(g, g)
    G1 = (X1 - X2) ** 2 - 1
    G2 = X1 + X2 - 2
    N = np.log10(np.maximum(np.abs(G1), np.abs(G2)) + 1e-12)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.75), sharey=True)
    for ax, K in zip(axes, [0.05, 0.005]):
        cs = ax.contourf(X1, X2, N, levels=np.linspace(-2.5, 1, 15), cmap="Blues_r", alpha=0.55,
                         extend="both")
        ax.plot(g, g, color=GREY, lw=1, ls=(0, (4, 3)))
        ax.text(1.95, 2.08, r"$\det J = 0$", color=GREY, fontsize=7.5, ha="right")
        path = trace_true_path(problem, K)
        ax.plot(path[:, 0], path[:, 1], color=INK, lw=2.2, label="exact homotopy path")
        times = pathway.constant_step(0.01, 0.5)
        for m in ["FE", "RK2", "BE"]:
            r = solve_hybrid(problem, problem.x0, K, times, m, "NR", 1e-5, max_it=50,
                             record_states=True)
            S = np.array(r.trajectory.states)
            R = np.array(r.refine.states)
            ax.plot(S[:, 0], S[:, 1], "--", color=METHOD[m], lw=1.3, marker="o", ms=4,
                    label=f"{m}, 3 steps")
            ax.plot(R[:, 0], R[:, 1], ":", color=METHOD[m], lw=1.1, marker=".", ms=4)
            ax.scatter(*R[-1], marker="*", s=90, color=METHOD[m], zorder=5, edgecolor="white", lw=0.5)
        for root in [(1.5, 0.5), (0.5, 1.5)]:
            ax.scatter(*root, s=60, facecolor="none", edgecolor=INK, lw=1.2, zorder=6)
        ax.scatter(1, 1, s=40, color=INK, zorder=6)
        ax.annotate(r"$x^{(0)}$", (1, 1), (1.05, 1.12), fontsize=8)
        ax.annotate("root A", (1.5, 0.5), (1.57, 0.55), fontsize=7.5)
        ax.annotate("root B", (0.5, 1.5), (0.57, 1.55), fontsize=7.5)
        ax.set(xlim=(-0.2, 2.2), ylim=(-0.2, 2.2), xlabel="$x_1$", title=f"K = {K}")
        ax.set_aspect("equal")
        ax.grid(False)
    axes[0].set_ylabel("$x_2$")
    axes[0].legend(loc="upper left", fontsize=6.8, framealpha=0.85, frameon=True)
    cb = fig.colorbar(cs, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label(r"$\log_{10}\|g(x)\|_\infty$", fontsize=8)
    save(fig, "tutorial_geometry")


# --------------------------------------------------------------------------- rescue histories
def history(r):
    t = r.trajectory
    hom = list(t.norms)
    nr = list(r.refine.norms[1:]) if r.refine else []
    return hom, nr


def outcome_label(pf, r):
    if not r.converged:
        return "fail"
    return "solved" if on_reference(pf, r.refine.x) else "other root"


def fig_rescues():
    runs = [("case36964", 0.05, 1e-3, "rescue: case36964, S2\n($\\Delta t_0=0.05$, $K=10^{-3}$)"),
            ("case3012wplimit", 0.005, 1e-4, "root fixed: case3012wplimit, S1\n($\\Delta t_0=0.005$, $K=10^{-4}$)"),
            ("case14limit", 0.005, 1e-4, "speed-up: case14limit, S1\n($\\Delta t_0=0.005$, $K=10^{-4}$)")]
    cfgs = [("paper", Options(), MAROON, "o"), ("OM", Options(multiplier=True), SKY, "s"),
            ("PC", Options(corrector=True), GREEN, "^"),
            ("OM + PC", Options(multiplier=True, corrector=True), CORAL, "D")]
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.1), sharey=True)
    summary = {}
    for ax, (name, dt0, K, title) in zip(axes, runs):
        pf = load(name)
        for label, opt, color, marker in cfgs:
            r = solve(pf, pf.flat_start(), K, [0.0, dt0, 1.0], opt)
            hom, nr = history(r)
            k_h = np.arange(len(hom))
            ax.plot(k_h, np.log10(hom), ":", color=color, marker=marker, ms=3.5, lw=1.1)
            k_n = np.arange(len(hom) - 1, len(hom) + len(nr))
            out = outcome_label(pf, r)
            ax.plot(k_n, np.log10(np.maximum([hom[-1]] + nr, 1e-16)), "-", color=color, marker=marker,
                    ms=3.5, lw=1.4, label=f"{label} ({out}, {r.factorizations} LU)")
            summary[(name, label)] = (out, r.factorizations, len(hom) - 1,
                                      r.refine.iterations if r.refine else None)
        ax.axhline(-8, color=GREY, lw=0.8, ls="--")
        ax.set_title(title, fontsize=8)
        ax.legend(fontsize=5.9, loc="lower left")
        ax.set_ylim(-14, 9)
    axes[0].set_ylabel(r"$\log_{10}\|g(x)\|_\infty$")
    fig.supxlabel("step index: homotopy points (dotted), then NR iterations (solid)", fontsize=8, y=-0.02)
    save(fig, "rescues")
    for k, v in summary.items():
        print("  rescue", k, v)


def fig_richardson():
    """Richardson-controlled paths on case6024 (S1), the one case only Richardson control solves in S1."""
    name, dt0, K = "case6024", 0.005, 1e-4
    pf = load(name)
    cfgs = [("paper (fixed path)", Options(), MAROON, "s"),
            ("richardson (BE)", Options(richardson=True), CORAL, "o"),
            ("richardson + FE", Options(richardson=True, step="FE"), BLUE, "^"),
            ("richardson + RK2", Options(richardson=True, step="RK2"), GREEN, "v")]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw=dict(wspace=0.3))
    for label, opt, color, marker in cfgs:
        r = solve(pf, pf.flat_start(), K, [0.0, dt0, 1.0], opt)
        t = np.array(r.trajectory.times)
        n = np.array(r.trajectory.norms)
        out = outcome_label(pf, r)
        its = r.refine.iterations if r.refine else None
        ax.plot(t, np.log10(n), "-", color=color, marker=marker, ms=3, lw=1.1,
                label=f"{label}: {len(t) - 1} pts, {r.factorizations} LU, {out}")
        if len(t) > 2:
            ax2.plot(t[1:-1], np.diff(t)[1:], "-", color=color, marker=marker, ms=3, lw=1.1)
        print("  richardson", label, "points", len(t) - 1, "LU", r.factorizations, out, "NR its", its,
              "||g(1)||", n[-1])
    ax.set(xlabel="homotopy parameter t", ylabel=r"$\log_{10}\|g(x(t))\|_\infty$",
           title="Mismatch along the path")
    fig.legend(*ax.get_legend_handles_labels(), fontsize=6.8, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.2))
    ax2.set(xlabel="t at the start of the step", ylabel="accepted step size", yscale="log",
            title="Accepted step sizes")
    save(fig, "richardson_path")


def fig_feasibility():
    """(K, dt0) maps of the paper's method on four cases, in the report's style (cached)."""
    from improvements.feasibility import DT0_VALUES, K_VALUES, feasibility_grid
    cases = ["case18482", "case36964", "case6024", "case_ACTIVSg2000limit"]
    cache = GEN / "feasibility.json"
    grids = json.loads(cache.read_text()) if cache.exists() else {}
    for c in cases:
        if c not in grids:
            grids[c] = feasibility_grid(load(c)).tolist()
            cache.write_text(json.dumps(grids))
            print("  feasibility", c, flush=True)
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#F4D6D2", "#FBE3B5", SKY])
    fig, axes = plt.subplots(1, 4, figsize=(7.5, 2.75), sharey=True, gridspec_kw=dict(wspace=0.08))
    lk = np.log10(K_VALUES)
    for ax, c in zip(axes, cases):
        g = np.array(grids[c])
        ax.imshow(g, origin="lower", cmap=cmap, vmin=0, vmax=2, aspect="auto",
                  extent=(lk[0] - 0.25, lk[-1] + 0.25, -0.5, len(DT0_VALUES) - 0.5))
        for d in (0.02,):
            xs = np.linspace(lk[0], lk[-1], 60)
            rows = np.interp(xs - np.log10(d), np.log10(DT0_VALUES), np.arange(len(DT0_VALUES)),
                             left=np.nan, right=np.nan)
            ax.plot(xs, rows, "--", color=INK, lw=0.9)
        for dt0, K in [(0.005, 1e-4), (0.05, 1e-3), (0.1, 2e-3)]:
            ax.plot(np.log10(K), DT0_VALUES.index(dt0), "o", color=INK, ms=4.5, mfc="white", mew=1.2)
        ax.set_title(f"{SHORT.get(c, c)}  ({int((g == 2).sum())}/{g.size})", fontsize=8.5)
        ax.set_xlabel(r"$\log_{10} K$")
        ax.set_xticks([-6, -4, -2])
        ax.grid(False)
    axes[0].set_yticks(range(len(DT0_VALUES)), [f"{d:g}" for d in DT0_VALUES], fontsize=7)
    axes[0].set_ylabel(r"$\Delta t_0$")
    fig.legend(handles=[Patch(color=SKY, label="reference solution"), Patch(color="#FBE3B5", label="another root"),
                        Patch(color="#F4D6D2", label="fails"),
                        Line2D([], [], color=INK, ls="--", label=r"$\delta=K/\Delta t_0=0.02$"),
                        Line2D([], [], color=INK, marker="o", mfc="white", ls="", label="settings S1, S2, S3")],
               ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.2), fontsize=7)
    save(fig, "feasibility")


def fig_adaptive():
    """Log every step the adaptive rule tries on case6748 (S1)."""
    name, dt0, K = "case6748", 0.005, 1e-4
    pf = load(name)
    attempts = []
    original = imp_path.try_step

    def logging_step(step, h, x, t, dt):
        x_new, norm = original(step, h, x, t, dt)
        attempts.append((t, dt, norm))
        return x_new, norm

    imp_path.try_step = logging_step
    try:
        r = solve(pf, pf.flat_start(), K, [0.0, dt0, 1.0], Options(multiplier=True, adaptive=True))
    finally:
        imp_path.try_step = original
    fixed = solve(pf, pf.flat_start(), K, [0.0, dt0, 1.0], Options())
    acc_t = r.trajectory.times
    acc_n = r.trajectory.norms
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    last_from = {}
    for i, (t, dt, n) in enumerate(attempts):
        last_from[round(t, 12)] = i
    start_norm = {round(t, 12): n for t, n in zip(acc_t, acc_n)}
    for i, (t, dt, n) in enumerate(attempts):
        ok = last_from[round(t, 12)] == i
        v = n if np.isfinite(n) else 1e6
        ax.annotate("", (t + dt, np.log10(v)), (t, np.log10(start_norm[round(t, 12)])),
                    arrowprops=dict(arrowstyle="-|>", color=GREEN if ok else CORAL, lw=0.9,
                                    ls="-" if ok else "--", shrinkA=2, shrinkB=2))
        if not ok:
            ax.scatter(t + dt, np.log10(v), marker="x", color=CORAL, s=22, zorder=4)
    ax.plot(acc_t, np.log10(acc_n), "o", color=GREEN, ms=5, zorder=5, label="accepted point")
    ax.plot(fixed.trajectory.times, np.log10(fixed.trajectory.norms), "s--", color=MAROON, ms=4,
            lw=1, label="paper's fixed path {0, 0.005, 1}")
    ax.plot([], [], "x", color=CORAL, label="rejected trial (norm grew)")
    ax.set(xlabel="homotopy parameter t", ylabel=r"$\log_{10}\|g(x(t))\|_\infty$",
           title="Adaptive path on case6748 (S1)")
    ax.legend(fontsize=7)
    save(fig, "adaptive_path")
    print("  adaptive accepted times", [round(t, 5) for t in acc_t], "LUs", r.factorizations,
          "trials", len(attempts), "converged", r.converged, "NR its",
          r.refine.iterations if r.refine else None, "| fixed converged", fixed.converged,
          "fixed x(1) norm", fixed.trajectory.norms[-1])


if __name__ == "__main__":
    st = systems()
    table_systems(st)
    fig_systems(st)
    fig_spy()
    fig_tutorial()
    fig_rescues()
    fig_adaptive()
    fig_richardson()
    fig_feasibility()
