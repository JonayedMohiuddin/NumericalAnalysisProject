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
experiments/            one script per result in the paper
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

## Other details

- The flat start sets all angles, including the slack angle, to 0. For NR from the case
  file we shift the stored angles so the slack is at 0. Both give the same solution.
- NR stops after 10 iterations (MATPOWER's default). A run counts as failed if it has not
  reached ||g|| < 1e-8 by then.
- In Figure 4 the mismatch is drawn in its own panel instead of on a second y axis.
