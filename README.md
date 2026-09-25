# Dynamic Homotopy for Ill-Conditioned Power Flow

CSE 402 project. This is a Python implementation of

A. Lima-Silva and F. D. Freitas, "Exploring a Dynamic Homotopy Technique to Enhance the
Convergence of Classical Power Flow Iterative Solvers in Ill-Conditioned Power System
Models", Energies 17(18):4642, 2024. doi:10.3390/en17184642

The paper is in `energies-17-04642-v3.pdf`. We rerun every table and figure of the
paper. Most results match closely; the ones that do not are listed under "Differences".
On top of that we test changes to the method (`improvements/`) and implement the
numerical parts of our project proposal (`extensions/`).

## Running

```
pip install -r requirements.txt
python run_all.py                           # all experiments
python run_all.py --only 3 4                # only some of them
python -m experiments.exp03_tables2to4      # a single experiment
python -m pytest tests
```

The test systems are downloaded from Zenodo into `data/cases/` the first time they are
needed (records 3514739 and 3491654, the same files the paper uses). Results are written
to `results/tables` and `results/figures`. Where the paper has a value, a table cell is
printed as `ours / paper`.

Experiments 1 to 10 (the paper) take about an hour. Experiment 11 (our improvements,
17 configurations x 170 runs) takes several hours on a laptop.

## The method

Power flow means solving g(x) = 0, where x holds the bus voltage angles and magnitudes.
Newton-Raphson (NR) often fails on large, badly conditioned grids when it starts from a
flat start (all voltages 1 pu, all angles 0). The paper builds the homotopy

    G(x, t) = t g(x) + (1 - t) K (x - x0)

which is easy to solve at t = 0 (the answer is x0) and equals the real problem at t = 1.
Differentiating G(x(t), t) = 0 gives the ODE

    dx/dt = -[t J(x) + (1 - t) K I]^-1 (g(x) - K (x - x0))

This ODE is integrated from t = 0 to 1 with only a few steps of forward Euler, RK2 or
backward Euler, which gives a rough answer. That answer is then used as the starting
point for NR or the fast decoupled method (FDXB). The first step is always taken with
backward Euler, because it adds delta = K / dt0 to the diagonal of the Jacobian and
keeps the step small.

## Code layout

```
dynhomotopy/            the paper's method, nothing else
  problem.py            g(x) = 0 interface, and the 2x2 example from Section 3.3
  datasets.py           the nine test systems and their download
  linalg.py             sparse LU with a counter for the number of factorizations
  hybrid.py             homotopy followed by NR or FDXB
  powerflow/
    case.py             reads MATPOWER .mat files, renumbers buses, finds PV/PQ buses
    network.py          Ybus, Sbus, and B' / B'' for FDXB
    model.py            power flow equations and Jacobian
  solvers/
    newton.py           Newton-Raphson
    fdxb.py             fast decoupled load flow, XB version
    gsh.py              GSH-NR, the static homotopy the paper compares against
  homotopy/
    fpv.py              G, Gx and Gt
    integrators.py      FE, BE, RK2 and the linearised first step
    pathway.py          the time points t_k
    dynamic.py          runs the integration along the path
improvements/           changes to the method, each behind a switch (see "Improvements")
  solve.py              solve() and the Options switches
  multiplier.py         NR with Iwamoto's optimal multiplier
  steps.py              step rules (BE-chord, FE, RK2, RK4) and the Newton corrector
  homotopies.py         the scaled and the Newton homotopy
  path.py               fixed, adaptive and Richardson-controlled time points
  cases.py              the wider test bed of 34 cases and their reference solutions
extensions/             the numerical parts of our proposal (see "Proposal features")
  gauss.py              Gauss elimination and LU with partial pivoting, from scratch
  scratch.py            BE path, NR and FDXB using that LU
  spectrum.py           power method and shifted inverse iteration
  tuning.py             golden-section search for delta and dt0
experiments/            one script per result
tests/
run_all.py
```

The solvers only use `g(x)` and `jacobian(x)`, so the 2x2 example and the 109,272-bus
system go through the same code.

| Result | Script | Output |
|---|---|---|
| Figures 1 and 2 | exp01_tutorial.py | fig1*.png, fig2*.png |
| Section 4.2.1 | exp02_first_step.py | sec421_first_step.md |
| Tables 2, 3, 4 | exp03_tables2to4.py | table2_FE.md, table3_RK2.md, table4_BE.md |
| Table 5 | exp04_table5.py | table5.md |
| Section 4.4 | exp05_sec44_pathways.py | sec44_pathways.md |
| Figure 3, Table 6 | exp06_fig3_table6.py | fig3*.png, table6.md |
| Figure 4 | exp07_fig4.py | fig4_case109272_states.png |
| Figure 5 | exp08_fig5.py | fig5_timestep_sensitivity.png |
| Table 7 | exp09_table7.py | table7.md, table7_native_base.md |
| Table 8 | exp10_table8.py | table8.md |
| Our improvements | exp11_improvements.py | improvements_*.md |
| Proposal: spectrum | exp12_spectrum.py | spectrum.md |
| Proposal: (K, dt0) maps | exp13_feasibility.py | feasibility.md, feasibility_maps.png |
| Proposal: tuning delta | exp14_tuning.py | tuning_per_case.md, tuning_rules.md |
| Proposal: own LU | exp15_scratch_lu.py | scratch_lu.md, gauss_vs_lu.md |

## Results

What matches the paper:

- Tables 4 and 5 (backward Euler) match almost digit for digit, including the NR
  iterations after the homotopy. For case18482 we get 532, 65, 33, 7.2, 3.7 along the
  path and then 0.04, 1.5e-5, 1e-11, the same as the paper.
- Tables 2 and 3 match wherever FE and RK2 do not blow up. Where they do blow up, they
  do so at the same place, but the size of the blow-up differs.
- Table 7 matches on the five large ill-conditioned cases in every column except two
  (GSH-NR on case27318 and BE(FDXB) on case109272, see below), and on case500limit and
  case2000limit in every column. Every converged run in Table 7 reaches the reference
  operating point (see "Checking the solution").
- Figures 4 and 5 show the same voltage paths. Bus 6 of case109272 ends at 1.0045 pu
  and 21.1 degrees for every time step.
- Figure 3 has the same shapes. The paper's one failure (FE with t3 = 0.25 on
  case36964) also fails here.
- Run times (median of 7 runs, relative to NR from the case file): BE followed by FDXB
  takes 66-93% of the NR time on the five large cases (paper: 57-85%), and BE followed
  by NR takes 102-122% (paper: 92-128%). In Table 6, case109272 is within 13
  percentage points of the paper in every row, and within 8 in most (for example BE
  163% against 169-171%).

Differences:

- case69limit and case141limit store baseMVA = 10. The paper's numbers for them in
  Tables 2-5 only come out on a 100 MVA base, so every experiment uses 100 MVA for these
  two cases (the Table 7 caption says so). On that base NR from the case file and from
  a flat start take 2 iterations where Table 7 says 4 and 3. On the 10 MVA base those
  two columns match, but then BE(NR) takes 2 instead of 1 (table7_native_base.md).
  The paper's numbers for these two cases do not fit a single base.
- FDXB: we report the number of P updates plus Q updates (for example 15 + 14 = 29),
  because that is what Table 7 counts. MATPOWER's fdpf returns one count per P/Q cycle
  and prints the P and Q counts separately. With this count every FDXB entry of
  Table 7 matches except the next one.
- Table 7 gives BE(FDXB) = 19 for case109272. We get 29, the same as for case54636,
  which is the same network repeated half as many times. In Table 2, FE on case18482 at
  t = 0.01 gives 3.7 against 9.7 in the paper, while every other cell in that column
  matches.
- GSH-NR is only described in words in the paper, so solvers/gsh.py is our own version
  of it. It matches case18482 (10 iterations) and case36964 (28 against 27), but fails
  on case27318 where the paper reports 25.
- Figure 3 and Table 6: RK2 with t3 = 0.20 fails on case36964, where the paper's plot
  converges.
- Section 4.4: BE works on almost every path, as the paper says. Unlike the paper, FE
  and RK2 also work on some 4 point paths, and BE followed by NR fails on case36964
  with only 3 points (FDXB still converges there).
- Tutorial: Figure 1(a) matches for eps = 0.005 and 0.01. The paper's curve for
  eps = 0.05 does not agree with the equations (one NR step from (1, 1.05) gives
  ||g|| of about 100, not 1000). In Figures 1(b)-(d) BE is still the best start for
  NR, but some iteration counts differ by 1 or 2. Some runs reach the other root
  (0.5, 1.5) of the system.
- In Table 8, RK2 costs 252-331% of plain NR here against 84-218% in the paper. Our
  timings come from Python and SciPy and the paper's from MATLAB, so only percentages
  are compared.

## Checking the solution

Power flow equations have more than one solution, and a mismatch below 1e-8 does not
tell them apart. For every case we take as reference the solution NR reaches from the
initial guess stored in the case file (`improvements.cases.reference_voltage`). A run
counts as solved only if its voltages are within 1e-4 pu of the reference.

This matters. With the paper's own method on our 34-case test bed, 14 of the 141 runs
that converge reach another root, for example case27318, case54636 and case109272 with
dt0 = 0.1 or with K = 1e-5 or 1e-3. With the paper's settings (dt0 = 0.005, K = 1e-4)
all of the paper's nine cases reach the reference.

## Improvements

The `improvements` package changes the paper's method without touching `dynhomotopy`.
Every change is a switch in `Options`, and with all switches off `improvements.solve`
gives exactly the result of the paper's method (a test checks this).

```python
from improvements import Options, solve

options = Options(multiplier=True, corrector=True)
result = solve(pf, pf.flat_start(), K=1e-4, times=[0, 0.005, 1], options=options)
```

| Switch | What it does | Where the idea comes from |
|---|---|---|
| `corrector=True` | after every path step, one Newton step on G(x, t) = 0 | predictor-corrector, Sec. 2.1.2, refs. [36, 37] |
| `multiplier=True` | scales every step of the final NR by Iwamoto's optimal multiplier | ref. [10] |
| `adaptive=True` | after t1, tries to jump to t = 1 and halves the step while the mismatch grows | future work (Sec. 5) |
| `richardson=True` | chooses the steps from a Richardson estimate of the local error (stops after 150 LUs) | our proposal |
| `step="FE"`, `"RK2"`, `"RK4"`, `"BE-chord"` | other step rules (BE-chord: 3 fixed point iterations on one LU, remark after eq. (19)) | future work, our proposal |
| `homotopy="scaled"`, `"newton"` | K D (x - x0) with D the Jacobian diagonal, or the Newton homotopy | ref. [29], ref. [36] |

### How we tested them

`experiments/exp11_improvements.py` runs each configuration on 34 cases: the paper's 9
plus 25 more from the same Zenodo records. Each case runs under 5 settings of (dt0, K)
with the path {0, dt0, 1}: the paper's default (S1), the two larger first steps of
Section 4.4 (S2, S3), a weaker K = 1e-5 (S4) and a stronger K = 1e-3 (S5). That is 170
runs per configuration. A run is solved if NR converges within 10 iterations and reaches
the reference operating point. Cost is the number of LU factorizations. Run times are
not compared, because each configuration runs only once.

| Configuration | Solved (of 170) | Other root | Lost vs paper | Mean LUs, S1 (paper: 8.04) |
|---|---|---|---|---|
| paper | 127 | 14 | - | 8.04 |
| PC | 152 | 2 | 0 | 8.43 |
| OM + PC | 152 | 3 | 0 | 7.57 |
| OM + PC + adaptive | 152 | 2 | 0 | 7.86 |
| PC + adaptive | 152 | 1 | 0 | 8.71 |
| OM | 143 | 12 | 0 | 7.14 |
| OM + adaptive | 148 | 10 | 1 | 7.46 |
| adaptive | 132 | 14 | 1 | 8.36 |
| richardson (BE) | 151 | 9 | 1 | 10.61 |
| OM + richardson | 156 | 8 | 0 | 9.75 |
| richardson + FE | 157 | 3 | 1 | 37.14 |
| richardson + RK2 | 158 | 4 | 1 | 47.50 |
| richardson + RK4 | 145 | 10 | 9 | 87.15 |
| scaled homotopy | 129 | 6 | 1 | 8.04 |
| BE-chord | 89 | 24 | 47 | 8.89 |
| newton homotopy | 75 | 10 | 56 | 9.00 |
| RK4 (fixed steps) | 13 | 1 | 114 | 7.00 |

OM is the multiplier, PC the corrector. The LU column only counts S1 cases solved by
both the configuration and the paper's method; for the last four rows these are fewer
cases, so their paper mean differs from 8.04. The full tables, including all five
settings, are in `results/tables/improvements_*.md`.

What we found:

- The corrector is the most useful change. It solves 25 more runs than the paper's
  method, never loses one, and almost removes convergence to other roots (2 instead of
  14), for about 0.4 more LUs per case. It also makes the result nearly independent of
  dt0 (32 of 34 cases in S1, S2 and S3).
- The multiplier lowers the cost (7.14 LUs instead of 8.04, about 11% fewer) but does
  not stop the method from landing on other roots (12). Together with the corrector
  it keeps the robustness of the corrector and most of the saving (7.57 LUs), which
  makes OM + PC the best overall choice.
- Richardson step control is the only change that solves case6024 at the paper's own
  setting, and case7092, case10595 and case12110 with K = 1e-5, where the paper's
  method diverges.
  With BE it costs about 30% more LUs than the paper's method, and 20% more together
  with the multiplier (OM + richardson, 156 solved). With FE or RK2 it solves slightly
  more runs (157, 158), but needs 4-6 times as many LUs.
- The paper's claim that BE is the best integrator still holds under error control:
  with Richardson control BE needs 10.6 LUs, FE 37, RK2 48 and RK4 87, and RK4 also
  solves fewer runs.
- Adaptive halving alone helps little (132). Combined with the corrector it gets very
  expensive when K is large (about 20 LUs in S5), because each rejected step costs two
  LUs.
- Extra fixed point iterations on one LU (the paper's own suggestion) make the start
  worse. Fixed-step RK4 and the Newton homotopy fail on most cases; the Newton homotopy
  has no K I term, which shows that this shift is what makes the method work. Scaling K
  by the Jacobian diagonal changes almost nothing.

`improvements_gains.md` lists every case the paper's method does not solve and says
whether it diverges, reaches another root, or is only slow (for case14limit and
case2736splimit it converges in 11 NR iterations if allowed more than 10).

## Proposal features

Our proposal promised a spectral study of the Jacobian, maps of where the method works
in the (K, dt0) plane, golden-section tuning, Richardson-controlled steps (see
"Improvements") and linear solvers written from scratch. All of them are implemented;
this is what they showed.

**Spectrum** (`extensions/spectrum.py`, `exp12`). At the flat start the eigenvalue of J
closest to zero is 4e-4 to 6e-3 on the twelve ill-conditioned cases, and 0.01 to 1.9 on
most limit cases (the two limit cases built from the 9241- and 13659-bus grids are
lower). The largest is 3e4 to 5e4 on the large grids, so the ratio
lambda_max / lambda_min reaches 1.2e8 (case109272). The paper's shift delta = 0.02 lifts
the smallest eigenvalue to about 0.02 and lowers the ratio to 1e6-2e6 on the
ill-conditioned cases, but hardly changes most limit cases, whose smallest eigenvalue is
already of that size or larger. That explains why the method behaves like plain NR on
those cases. Inverse iteration
matches ARPACK to 4 digits on 33 of 34 cases. The power method is off by up to 19% on
the networks built from copies of the same grid, where the largest eigenvalues come in
near-equal clusters and it converges very slowly.

**Feasible region** (`exp13`, `feasibility_maps.png`). The method works in a band along
the diagonals delta = K / dt0 = constant. For case18482 and case2000limit almost the
whole grid works. For case36964 only 29 of 88 cells work, and no cell with dt0 = 0.1
or 0.2 works. This is why the paper found that dt0 = 0.1 with K = 0.002 fails even though
delta is unchanged. For case6024 the band lies at larger delta (0.03 to 32), and the
paper's delta = 0.02 line runs along its edge, where runs reach other roots.

**Tuning delta** (`extensions/tuning.py`, `exp14`). Golden-section search over
log10(delta), then over log10(dt0), minimising total LUs, needs 21 solves per case.
The tuned values show that the paper's delta suits some of its own cases (0.025 on
case27318, case54636 and case109272) but is far too small for the other ill-conditioned
cases (tuned delta 1 to 5). We fitted two rules on the paper's 9 cases and tested them on the other 25:

| How delta is chosen | Solved, 9 training cases | Solved, 25 test cases | Mean LUs, test |
|---|---|---|---|
| paper, delta = 0.02 | 9 | 19 | 9.47 |
| fixed tuned value, delta = 0.011 | 9 | 19 | 9.53 |
| spectral rule, delta = 1.86 lambda_min | 5 | 20 | 11.00 |
| golden search on each case | 9 | 23 | 9.13 |

The spectral rule does not work. On the ill-conditioned cases the tuned delta is 50 to
800 times lambda_min, so no single factor fits, and on the limit cases delta hardly
matters at all. The rule loses four training cases. Searching
each case solves four more test cases, but the search itself costs about 21 solves.

**Own linear solvers** (`extensions/gauss.py`, `extensions/scratch.py`, `exp15`). The
paper's method with our dense LU with partial pivoting gives the same iterations as
with SciPy's sparse LU and the same solution to 5e-12, on the tutorial and the 69- to
2000-bus cases. It is much slower (424 s against 0.04 s for NR on the 2000-bus case),
because a dense LU costs O(n^3) and ignores the sparsity of the Jacobian. Factorising
once and reusing the LU is 14 times faster than repeating Gauss elimination for 20
right-hand sides, which is why FDXB factorises B' and B'' only once.

## Other details

- The flat start sets all angles, including the slack angle, to 0. For NR from the case
  file we shift the stored angles so the slack is at 0. Both give the same solution.
- NR stops after 10 iterations (MATPOWER's default). A run counts as failed if it has not
  reached ||g|| < 1e-8 by then.
- In Figure 4 the mismatch is drawn in its own panel instead of on a second y axis.
