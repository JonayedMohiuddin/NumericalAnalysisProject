# Dynamic Homotopy for Ill-Conditioned Power Flow

CSE 402 project. This is a Python implementation of

A. Lima-Silva and F. D. Freitas, "Exploring a Dynamic Homotopy Technique to Enhance the
Convergence of Classical Power Flow Iterative Solvers in Ill-Conditioned Power System
Models", Energies 17(18):4642, 2024. doi:10.3390/en17184642

The paper is in `energies-17-04642-v3.pdf`. We reproduce all of its tables and figures.

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
to `results/tables` and `results/figures`. Every table cell is printed as `ours / paper`.
Running everything takes around 30 minutes, mostly for the timing tables.

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
backward Euler, because it adds K / dt0 to the diagonal of the Jacobian and keeps the
step small.

## Code layout

```
dynhomotopy/
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
improvements/           our additions, built on top of dynhomotopy (see "Improvements")
  solve.py              solve() and the Options switches
  multiplier.py         NR with Iwamoto's optimal multiplier
  steps.py              step rules (BE-chord, RK4) and the Newton corrector
  homotopies.py         the scaled and the Newton homotopy
  path.py               fixed and adaptive time points
  cases.py              the wider test bed of 34 cases
experiments/            one script per result in the paper, plus exp11 for our improvements
tests/
run_all.py
```

The solvers only use `g(x)` and `jacobian(x)`, so the 2x2 example and the 109,272-bus
system go through the same code.

| Paper result | Script | Output |
|---|---|---|
| Figures 1 and 2 | exp01_tutorial.py | fig1*.png, fig2*.png |
| Section 4.2.1 | exp02_first_step.py | sec421_first_step.md |
| Tables 2, 3, 4 | exp03_tables2to4.py | table2_FE.md, table3_RK2.md, table4_BE.md |
| Table 5 | exp04_table5.py | table5.md |
| Section 4.4 | exp05_sec44_pathways.py | sec44_pathways.md |
| Figure 3, Table 6 | exp06_fig3_table6.py | fig3*.png, table6.md |
| Figure 4 | exp07_fig4.py | fig4_case109272_states.png |
| Figure 5 | exp08_fig5.py | fig5_timestep_sensitivity.png |
| Table 7 | exp09_table7.py | table7.md |
| Table 8 | exp10_table8.py | table8.md |
| Our improvements | exp11_improvements.py | improvements_*.md |

## Results

Most results match the paper closely:

- Tables 4 and 5 (backward Euler) match almost digit for digit, including the NR
  iterations after the homotopy. For case18482 we get 532, 65, 33, 7.2, 3.7 along the
  path and then 0.04, 1.5e-5, 1e-11, the same as the paper.
- Tables 2 and 3 match wherever FE and RK2 do not blow up. Where they do blow up, they
  do so at the same place, but the size of the blow-up differs.
- Table 7 matches for NR from the case file, NR from a flat start (fails on all large
  cases), and every BE / RK2 run with NR or FDXB.
- Figures 4 and 5 show the same voltage paths. Bus 6 of case109272 ends at 1.0045 pu
  and 21.1 degrees for every time step.
- Figure 3 has the same shapes, and the one failure mentioned in the paper (FE with
  t3 = 0.25 on case36964) also fails here.
- For run times, BE followed by FDXB is faster than plain NR on all large cases (66-91%
  of the NR time, the paper reports 57-85%). BE followed by NR costs about the same as
  plain NR. Table 6 for case109272 is within a few percent of the paper.

Differences we found:

- case69limit and case141limit store baseMVA = 10. The paper's numbers for them in
  Tables 2-5 only come out if 100 MVA is used, so we use 100 MVA for these two cases.
  But the paper's plain NR iteration counts for them in Table 7 match the 10 MVA base
  (see table7_native_base.md), so the paper is not consistent here.
- FDXB iterations are counted the way MATPOWER prints them, as P steps plus Q steps
  (for example 15 + 14 = 29). With that count all FDXB numbers in Table 7 match.
- Table 7 gives BE(FDXB) = 19 for case109272. We get 29, the same as for case54636,
  which is the same network repeated half as many times. We think 19 is a typo. In
  Table 2, FE on case18482 at t = 0.01 gives 3.7 against 9.7 in the paper, while every
  other cell in that column matches.
- GSH-NR is only described in words in the paper, so solvers/gsh.py is our own version
  of it. It matches case18482 (10 iterations) and case36964 (28 against 27), but fails
  on case27318 where the paper reports 25.
- Section 4.4: BE works on almost every path, as the paper says. Unlike the paper, FE
  and RK2 also work on some 4 point paths, and BE followed by NR fails on case36964
  with only 3 points (FDXB still converges there).
- Tutorial: Figure 1(a) matches for eps = 0.005 and 0.01. The paper's curve for
  eps = 0.05 does not agree with the equations (one NR step from (1, 1.05) gives
  ||g|| of about 100, not 1000). In Figures 1(b)-(d) BE is still the best start for
  NR, but some iteration counts differ by 1 or 2. Some runs reach the other root
  (0.5, 1.5) of the system.
- RK2 is slower here than in the paper (250-330% of plain NR). Our timings come from
  Python and SciPy and the paper's from MATLAB, so only the percentages can be compared.
  We report the median of several runs because single runs varied a lot on a laptop.

## Improvements

The `improvements` package adds changes to the paper's method without touching
`dynhomotopy`. Every change is a switch in `Options`, and with all switches off
`improvements.solve` gives exactly the result of the paper's method (a test checks this).

```python
from improvements import Options, solve

options = Options(multiplier=True, adaptive=True)
result = solve(pf, pf.flat_start(), K=1e-4, times=[0, 0.005, 1], options=options)
```

| Switch | What it does | Where the idea comes from | Code |
|---|---|---|---|
| `multiplier=True` | scales every step of the final NR by Iwamoto's optimal multiplier | ref. [10] | multiplier.py |
| `corrector=True` | after every path step, one Newton step on G(x, t) = 0 | predictor-corrector, Sec. 2.1.2, refs. [36, 37] | steps.py |
| `adaptive=True` | after t1, tries to jump to t = 1 and halves the step while the mismatch grows | future work (Sec. 5) | path.py |
| `step="BE-chord"` | BE with 3 fixed point iterations on one LU | remark after eq. (19) | steps.py |
| `step="RK4"` | classical RK4 after a BE first step | future work (Sec. 5) | steps.py |
| `homotopy="scaled"` | K D (x - x0) with D the Jacobian diagonal | ref. [29] | homotopies.py |
| `homotopy="newton"` | Newton homotopy g(x) - (1 - t) g(x0) | future work, ref. [36] | homotopies.py |

The first three helped in our tests; the last four did not and are kept for comparison.

### How we tested them

`experiments/exp11_improvements.py` runs every configuration on 34 cases: the paper's 9
plus 25 more from the same Zenodo records (`improvements/cases.py`). Each case runs under
5 settings of (dt0, K) with the path {0, dt0, 1}: the paper's default (S1), the two larger
first steps of Section 4.4 (S2, S3), a weaker K = 1e-5 (S4) and a stronger K = 1e-3
(S5). That is 170 runs per configuration. A run is solved if NR then converges within
10 iterations, as in the paper. Cost is the number of LU factorizations, which does not
depend on the machine.

| Configuration | Solved (of 170) | Lost vs paper | Mean LUs, S1 | Time vs paper, S1 |
|---|---|---|---|---|
| paper | 141 | - | 8.30 | 1.00 |
| OM | 155 | 0 | 7.37 | 0.92 |
| OM + adaptive | 158 | 1 | 7.90 | 0.92 |
| OM + PC | 155 | 0 | 7.70 | 0.98 |
| PC | 154 | 0 | 8.60 | 1.00 |
| OM + PC + adaptive | 154 | 1 | 7.97 | 0.98 |
| PC + adaptive | 153 | 1 | 8.87 | 1.01 |
| adaptive | 146 | 1 | 8.83 | 0.99 |
| scaled homotopy | 135 | 6 | 8.30 | 1.04 |
| BE-chord | 113 | 30 | 8.75 | 1.06 |
| newton homotopy | 85 | 56 | 8.94 | 1.10 |
| RK4 | 14 | 127 | 7.00 | 2.30 |

OM is the multiplier, PC the corrector. "Mean LUs" only counts cases solved by both the
configuration and the paper's method (the paper's own mean on the same cases is 8.30,
except for BE-chord and RK4, which solve fewer cases). The full tables are in
`results/tables/improvements_*.md`.

What we found:

- The multiplier is the safest improvement. It solves 14 more runs, never loses one the
  paper's method solves, and needs about 11% fewer LU factorizations. Most of the saving
  is on the stressed limit cases, where NR needs 2 or 3 fewer iterations.
- The multiplier together with adaptive steps solves the most runs (158). It is the only
  configuration that solves case6024 (with the larger first steps S2 and S3) and case6748,
  where the paper's method diverges. It loses one run: case2383wplimit with strong K.
- The corrector rescues the same divergent case36964 runs as the multiplier, and makes
  the result almost independent of dt0, but it costs one extra LU per path point.
- Adding the corrector to the multiplier gives nothing extra. Adding the corrector to
  adaptive steps can get expensive (19 LUs on average with strong K), because every
  rejected step now costs two LUs.
- Extra fixed point iterations on one LU (the paper's own suggestion) make the start worse.
  RK4 and the Newton homotopy fail on most cases. The Newton homotopy has no K I term,
  which shows that this shift is what makes the method work. Scaling K by the Jacobian
  diagonal changes nothing at the default setting and loses cases with larger steps.

Real rescues and speed-ups are different things. For case14limit and case2736splimit the
paper's method does converge if NR is allowed 30 iterations (it needs 11), so there the
improvements only make it faster. `improvements_gains.md` marks every such case. The
paper's method really diverges on case36964 (S2, S3, S4), case6024 (S2, S3), case6748
(S1, S2), case10595 and case12110 (S4), and each of these is solved by at least one
configuration above.

## Other details

- The flat start sets all angles, including the slack angle, to 0. For NR from the case
  file we shift the stored angles so the slack is at 0. Both give the same solution.
- NR stops after 10 iterations (MATPOWER's default). A run counts as failed if it has not
  reached ||g|| < 1e-8 by then.
- In Figure 4 the mismatch is drawn in its own panel instead of on a second y axis.
