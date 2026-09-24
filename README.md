# Dynamic homotopy for ill-conditioned power flow: a reproduction

This repository is a Python reimplementation of

> A. Lima-Silva, F. D. Freitas, "Exploring a Dynamic Homotopy Technique to Enhance the
> Convergence of Classical Power Flow Iterative Solvers in Ill-Conditioned Power System
> Models," *Energies* 17(18):4642, 2024. doi:10.3390/en17184642 (`energies-17-04642-v3.pdf`)

Every table and figure in the paper is regenerated from scratch. The pieces are a
MATPOWER-equivalent power-flow model, Newton-Raphson, fast-decoupled XB, the
dynamic homotopy with FE / RK2 / BE, and the GSH-NR baseline. Each result is
printed next to the value the paper published.

## Quick start

```bash
pip install -r requirements.txt
python run_all.py                 # all experiments, in the paper's order
python run_all.py --only 3 4      # just Tables 2-5
python -m pytest tests            # sanity and regression tests
```

The nine test systems (Table 1) are downloaded from Zenodo into `data/cases/` on
first use. The sources are records [3514739](https://zenodo.org/records/3514739)
(ill-conditioned) and [3491654](https://zenodo.org/records/3491654) (limit cases).
Outputs go to `results/tables/` (CSV + Markdown) and `results/figures/` (PNG).

## The method in one picture

```
flat start x0 ──► dynamic homotopy (integrate the ODE along γ(t), 0 ≤ t ≤ 1) ──► x̂(0) = x(1) ──► NR or FDXB ──► x*
                    G(x,t)  = t·g(x) + (1−t)·K(x − x0)                  (12), (22)
                    Gx      = t·J(x) + (1−t)·K·I                          (23)
                    Gt      = g(x) − K(x − x0)                            (24)
                    dx/dt   = −Gx⁻¹·Gt                                    (15)
                    first step t1 = Δt0 by BE:  Δx(0) = −[J + ((1−Δt0)/Δt0)·K·I]⁻¹ g(x0)   (25)/(27)
                    later steps: FE (16), BE with one fixed-point iteration (19), RK2 (20)-(21)
```

## Code map (paper → code)

| Paper | Module |
|---|---|
| Test systems, Table 1 (Sec. 4.1) | `dynhomotopy/datasets.py` (registry + Zenodo download) |
| MATPOWER case, `ext2int`, `bustypes` | `dynhomotopy/powerflow/case.py` |
| Ybus, Sbus, B′/B″ (`makeYbus`, `makeSbus`, `makeB`) | `dynhomotopy/powerflow/network.py` |
| Power-balance equations (1)-(2), Jacobian (3) | `dynhomotopy/powerflow/model.py` |
| Generic g(x) = 0 interface; tutorial system (Sec. 3.3) | `dynhomotopy/problem.py` |
| Newton-Raphson (4) | `dynhomotopy/solvers/newton.py` |
| Fast decoupled XB (5)-(8) | `dynhomotopy/solvers/fdxb.py` |
| GSH-NR static homotopy [12] (Tables 7-8) | `dynhomotopy/solvers/gsh.py` |
| FPV homotopy G, Gx, Gt (12), (22)-(24) | `dynhomotopy/homotopy/fpv.py` |
| FE (16), BE (17)-(19), RK2 (20)-(21), linear first step (26)-(27) | `dynhomotopy/homotopy/integrators.py` |
| Pathway instants t_k (Sec. 4.2) | `dynhomotopy/homotopy/pathway.py` |
| Integrating the IVP along the pathway (Sec. 4.3) | `dynhomotopy/homotopy/dynamic.py` |
| Homotopy followed by NR / FDXB | `dynhomotopy/hybrid.py` |
| LU factorisation (counted) | `dynhomotopy/linalg.py` |

Every solver works on the abstract `NonlinearProblem` (`g`, `jacobian`), so the
2×2 tutorial and a 109,272-bus grid go through the same code.

## Experiment map (paper result → script → output)

| # | Paper result | Script (`experiments/`) | Output |
|---|---|---|---|
| 1 | Fig. 1(a)-(d), Fig. 2 (tutorial) | `exp01_tutorial.py` | `fig1*.png`, `fig2*.png`, `fig1_fig2_tutorial_*.md` |
| 2 | Sec. 4.2.1 first time step | `exp02_first_step.py` | `sec421_first_step.md` |
| 3 | Tables 2, 3, 4 | `exp03_tables2to4.py` | `table2_FE.md`, `table3_RK2.md`, `table4_BE.md` |
| 4 | Table 5 | `exp04_table5.py` | `table5.md` |
| 5 | Sec. 4.4 pathway study (text) | `exp05_sec44_pathways.py` | `sec44_pathways.md` |
| 6 | Fig. 3(a)-(c), Table 6 | `exp06_fig3_table6.py` | `fig3*.png`, `table6.md` |
| 7 | Fig. 4 | `exp07_fig4.py` | `fig4_case109272_states.png`, `fig4_states.md` |
| 8 | Fig. 5 | `exp08_fig5.py` | `fig5_timestep_sensitivity.png`, `fig5_timestep.md` |
| 9 | Table 7 | `exp09_table7.py` | `table7.md`, `table7_native_base.md` |
| 10 | Table 8 | `exp10_table8.py` | `table8.md` |

Shared pieces: `common.py` (paths, formatting, timing), `paper_values.py` (the
published numbers), `norm_tables.py` (Tables 2-5), `setups.py` (the Table 7/8
configurations), `plots.py` (figure style).

## How closely the results match

Each generated table prints `ours / paper` in every cell.

| Result | Agreement |
|---|---|
| Flat-start mismatches (Table 2-5, t = 0) | All nine cases match (532, 200, 13, 54, 0.012, 0.006). |
| **Table 4** (BE + NR) | Essentially digit for digit on all nine cases, including every NR iteration (e.g. case18482: 532, 65, 33, 7.2, 3.7 → 0.04, 1.5e-5, 1e-11). |
| **Table 5** | Essentially digit for digit for all nine case/solver rows. |
| Tables 2-3 (FE, RK2) | Every converging trajectory matches. Divergent runs diverge in the same places, with different blow-up magnitudes (these are chaotic). One cell differs: FE case18482 at t = 0.01 is 3.7 vs 9.7. The paper's own case36964 (the same network ×4) gives 3.5. |
| Sec. 4.2.1 | Explicit first steps explode (‖g‖ up to 1e12), BE/LIN give 65 / 25, and BE ≡ linear approximation (27). |
| **Table 7** | NR-MAT, NR-flat, BE(NR), RK2(NR), BE(FDXB), RK2(FDXB) match on all ill-conditioned cases. Exception: case109272 BE(FDXB) is 29 vs 19; case54636 (the same network ×4) is 29 in both. GSH-NR matches 18482 (10) and 36964 (28 vs 27) and fails where the paper fails, but also fails on case27318 (paper: 25). See note 3. |
| Fig. 4 | Same trajectories (bus-6 angle ≈ 60° at t = 1 → 21.1°; NR in 4 iterations). |
| Fig. 5 | Same behaviour. dt = 0.1 hands NR ≈ 1.005 pu / 28°, and every dt converges to 1.0045 pu / 21.1°. |
| Fig. 3 | Same shapes. The paper's single failure (FE, t3 = 0.25, case36964) is reproduced. Here RK2 with t3 = 0.20 on case36964 also fails. |
| Sec. 4.4 text | Holds for BE, which works with 3-4 points almost everywhere (BE+NR fails only on case36964 with 3 points; FDXB still converges there). Unlike the paper's statement, FE/RK2 also succeed on some 4-point pathways. |
| Table 6 (CPU %) | case109272 within a few points of the paper (FE 145/145/164 vs 146/145/170; RK2 ≈228 vs 226-237; BE ≈166 vs 169-171). |
| Table 8 (CPU %) | Headline reproduced: BE(FDXB) is cheaper than the NR-MAT reference on every case (66-91 % vs the paper's 57-85 %), and BE(NR) costs about the same (102-123 % vs 92-128 %). RK2 is costlier here than in MATLAB (≈250-330 %). |
| Figs. 1-2 (tutorial) | Fig. 1(a) matches for ε = 0.005 / 0.01. The hybrid runs keep the conclusion (BE gives the fastest NR), but some counts differ by 1-2. See note 4. |

## Reproduction notes and deviations

These are the places where the paper leaves something unstated, and what this
implementation does instead.

1. **Base MVA of `case69limit` and `case141limit`.** Both files store
   `baseMVA = 10`. With that base, their flat-start mismatch is 0.124 / 0.064
   p.u., but the paper reports 0.012 / 0.006. On a 100 MVA base, every entry of
   Tables 2-5 and BE(NR) = 1 in Table 7 reproduce. The registry therefore uses
   100 MVA for these two cases (`paper_base_mva`). The paper's NR-flat / NR-MAT
   counts for them (4 and 3) only match the native 10 MVA base, which
   `table7_native_base.md` shows. The paper's own numbers for these two cases
   are therefore mutually inconsistent.
2. **FDXB iteration count.** The paper counts P- and Q-half-iterations
   separately, the way MATPOWER's `fdpf` reports them (e.g. 15 P + 14 Q = 29).
   `fast_decoupled_xb` counts the same way.
3. **GSH-NR** (ref. [12]) is described in this paper only in prose. It is
   reconstructed as a fictitious-shunt compensation, which makes the flat start
   exact at h = 0, plus impedance scaling of the slack's branches by δ (see the
   docstring of `solvers/gsh.py`). With Δh = 1, δ = 1 it reduces to NR-flat, as
   the limit-case rows of Table 7 require.
4. **Tutorial (Figs. 1-2).** NR alone matches Fig. 1(a) for ε = 0.005 and 0.01
   (11 and 10 iterations). The ε = 0.05 curve of the paper is not consistent
   with the equations: one NR step from (1, 1.05) gives ‖g‖ ≈ 10², not the ≈10³
   plotted. The hybrid curves keep the paper's conclusion (BE gives the best NR
   start). A few per-curve iteration counts differ by 1-2, because the paper
   does not state every detail of its tutorial runs. Some runs reach the
   system's other root (0.5, 1.5), which is equally valid.
5. **CPU times** (Tables 6, 8) come from Python/SciPy (SuperLU), not MATLAB
   (UMFPACK), so only the ratios to NR-MAT are comparable. The paper averages
   100 runs. Here the methods run round-robin for `--reps` rounds and the
   median is reported, because single-run laptop timings drifted by up to 4×.
   On small systems (the 2000-bus case) Python's fixed per-step overhead
   dominates, so the ratios there say little about the method.
6. **Fig. 4** is drawn as stacked panels on a shared x-axis rather than with the
   paper's second y-axis for the mismatch. The data are the same.
7. The flat start sets every angle to zero, including the slack's. NR-MAT
   shifts the case file's angles so the slack is at 0. Both are exact
   re-referencings and leave g unchanged.
8. The NR iteration limit is MATPOWER's default of 10. A run that does not
   reach 1e-8 within it, or that hits a singular Jacobian or a non-finite
   state, counts as a fail.
