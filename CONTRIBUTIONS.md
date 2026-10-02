# Team Contributions

CSE 402 Numerical Analysis Sessional, project on

A. Lima-Silva and F. D. Freitas, "Exploring a Dynamic Homotopy Technique to Enhance the
Convergence of Classical Power Flow Iterative Solvers in Ill-Conditioned Power System
Models", Energies 17(18):4642, 2024.

Team: Jonayed, Hasbi, Apurbo, Sizan and Sadman.

## What the team did

The paper deals with a hard problem in power system analysis. Power flow means solving
the nonlinear equations g(x) = 0, where x holds the voltage angle and magnitude of every
bus. The usual method, Newton-Raphson (NR), starts from a "flat start" (all voltages
1 pu, all angles 0). On large, badly conditioned grids it diverges from there, even
though a solution exists.

The paper's idea is to reach the solution in two stages:

1. Build the homotopy G(x, t) = t g(x) + (1 - t) K (x - x0). At t = 0 its solution is
   simply the flat start x0, and at t = 1 it is the real power flow problem.
   Differentiating G(x(t), t) = 0 turns the problem into an ODE, which is integrated
   from t = 0 to t = 1 with a few steps of forward Euler, RK2 or backward Euler.
2. The rough point reached at t = 1 is handed to NR, or to the fast decoupled method
   (FDXB), which then converges quickly.

The paper tests this on grids of up to 109,272 buses. Its main claim is that backward
Euler with only two path points is enough to make NR converge where it fails from a
flat start.

As a team we did five things:

1. **Reproduced the paper.** We wrote the whole method from scratch in Python: power
   flow model, solvers, homotopy and integrators. We reran every table and figure of the
   paper on the same public test systems. Most results match the paper digit for digit,
   and every difference is documented.
2. **Checked the answers properly.** Power flow has several solutions, and a small
   mismatch does not show which one was found. We compare every result with a reference
   operating point, and found that some "converged" runs, including some of the paper's
   own method, reach the wrong solution.
3. **Improved the method.** We tried ten changes suggested by the paper and its
   references on a larger set of 34 test grids. Three clearly help: a Newton corrector,
   Iwamoto's optimal multiplier and Richardson step control. Each change is a separate
   switch that can be turned on alone or combined with others.
4. **Delivered our proposal.** Our proposal presentation stated a spectral study of the Jacobian,
   maps of where the method works, golden-section tuning, error-controlled steps, and
   linear solvers written by hand. All of these are implemented, and each has a
   measured result, including the ones that did not work.
5. **Generated the report.** Two scripts in `report/scripts/` build every figure and
   table of the project report (`report/main.tex`): one runs the project code, the
   other reads the saved results. So the report can be rebuilt after any rerun.

## Who did what

All five members contributed equally. Every member owns one part of the core
code, one part of the paper's results, one improvement and one numerical study from the
proposal, the report figures and tables built from that work, and the tests for their own
code. The table gives an overview, by roll number. The sections below explain each
member's work in full.

| Roll | Member | Core code | Paper results reproduced | Improvements and studies | Report |
|---|---|---|---|---|---|
| 2105039 | Sizan | table output, timing benchmark | Tables 2 to 6, Figure 3 | FE, RK2, RK4 and BE-chord step rules; golden-section tuning of delta | Tables 2 to 5 in LaTeX, fidelity plot, tuning figure |
| 2105040 | Sadman | sparse LU counter, result types, figure style | Section 4.4, Figure 5 | Richardson step control, scaled homotopy; Gauss elimination and LU from scratch; command-line tool | Jacobian sparsity, Section 4.4 grid, Figure 5 cost, Richardson paths |
| 2105057 | Apurbo | homotopy, integrators, hybrid method | Figures 1, 2, 4, Section 4.2.1 | `solve()` with switches, Newton corrector, adaptive steps; spectrum of the Jacobian | tutorial geometry, first step, adaptive path, spectrum figures |
| 2105058 | Hasbi | NR, FDXB, GSH-NR | Tables 7, 8 | optimal multiplier, Newton homotopy; NR and FDXB with our own LU | timing, cost against robustness, per-case map, rescue histories |
| 2105060 | Jonayed | power flow model, test systems, experiment runner | flat-start check and base MVA finding, root check of Table 7, audit of the paper | wider test bed, solution (root) check, the comparison study of all changes; feasibility maps | test-systems table and figure, audit table, feasibility figure, study summary figure and table |

---

## Jonayed

Jonayed built the power flow model that every other part runs on, the checks that decide
whether a result is right, and the study that compares every member's changes. His root
check is behind the project's main finding: that 14 of the paper's 141 converged runs on
the wider test bed reach the wrong solution. He also set up and maintains the project's
public GitHub repository, where the code of all five members comes together, and his
`g(x)` / `jacobian(x)` interface is what every other member's solver, integrator and
modification is written against.

### Core code: the power flow model

- `dynhomotopy/datasets.py` lists the paper's nine test systems (five ill-conditioned
  grids of 18,482 to 109,272 buses and four stressed "limit" cases). It downloads them
  from the two Zenodo records the paper cites the first time they are needed.
- `dynhomotopy/powerflow/case.py` reads the MATPOWER `.mat` files and passes them to
  PYPOWER's `ext2int` and `bustypes`, which remove isolated buses and switched-off
  generators and branches, renumber the buses 0 to n-1, and mark each bus as slack,
  PV or PQ.
- `dynhomotopy/powerflow/network.py` wraps PYPOWER's `makeYbus`, `makeSbus` and `makeB`
  for the bus admittance matrix Ybus, the specified injections Sbus, and the B' and B''
  matrices of the fast decoupled method.
- `dynhomotopy/powerflow/model.py` defines the power flow equations (equations 1 and 2 of
  the paper). The state is x = [angles of PV and PQ buses, magnitudes of PQ buses], the
  mismatch is g(x) = V conj(Ybus V) - Sbus, and the Jacobian (equation 3) is assembled
  in sparse form from the derivatives of the bus injections given by PYPOWER's `dSbus_dV`.
- `dynhomotopy/problem.py` defines the interface `g(x)`, `jacobian(x)` that every solver
  uses, and the 2x2 example of Section 3.3. Because all solvers are written against this
  interface, the same code runs on the 2x2 example and on the 109,272-bus grid.

### Paper results: checking the model

- A finite-difference test checks the Jacobian column by column.
- The flat-start mismatch of all nine cases was compared with the paper (532, 200, 13,
  54, 0.012 and 0.006 pu). Seven matched at once. The two small limit cases came out
  exactly 10 times too large. Jonayed traced this to their files storing a 10 MVA base,
  while the paper's numbers for Tables 2 to 5 only come out on a 100 MVA base. The
  paper's own numbers for these two cases do not fit a single base, so the base used is
  stated wherever it matters.
- `experiments/exp09_table7.py` marks every converged run of Table 7 that reaches a
  different root, using his solution check; none does, so the Table 7 counts are not
  inflated by wrong roots.
- The audit of the paper: Jonayed collected the eight places where the paper is
  inconsistent or under-specified (his base MVA finding, Hasbi's FDXB count for
  case109272, Apurbo's eps = 0.05 curve, Sizan's Table 2 cell and Table 6 run, Sadman's
  Section 4.4 paths, the GSH-NR description and the case2000limit mismatch), checked each
  against our own code first, and wrote the audit table of the report with the evidence
  and our handling of each.

### Improvements: test bed, solution check, comparison study and feasibility maps

- `improvements/cases.py` adds 25 more test grids from the same Zenodo records (34 in
  total), so that changes to the method are tested on more than the paper's nine cases.
- The same file holds the solution check. `reference_voltage` computes the operating
  point NR reaches from the initial guess stored in each case file, and `on_reference`
  tests whether a result is within 1e-4 pu of it. Power flow equations have several
  solutions, and this check showed that 14 of the 141 runs where the paper's method
  converged on the wider test bed had reached a different root. Every entry of Table 7
  was confirmed to reach the correct solution.
- `experiments/exp11_improvements.py` is the study that compares everyone's changes. It
  runs each configuration on the 34 grids under 5 settings of dt0 and K (170 runs), and
  counts a run as solved only if NR converges within 10 iterations and reaches the
  reference solution. For every case the paper's method does not solve, it reruns the
  paper's method with 30 NR iterations, to tell apart a real divergence, a wrong root
  and a merely slow run. The main results:

  | Configuration | Solved (of 170) | Wrong root | Mean LUs |
  |---|---|---|---|
  | paper's method | 127 | 14 | 8.04 |
  | corrector (Apurbo) | 152 | 2 | 8.43 |
  | multiplier (Hasbi) + corrector (Apurbo) | 152 | 3 | 7.57 |
  | multiplier (Hasbi) | 143 | 12 | 7.14 |
  | multiplier + Richardson (Sadman) | 156 | 8 | 9.75 |

- `run_all.py` runs all 15 experiments, or any selection of them, so the whole project
  can be reproduced with one command.

`improvements/feasibility.py` runs the paper's method on a grid of K and dt0 values and
marks every cell as solved, wrong root or failed, using the solution check.
`experiments/exp13_feasibility.py` draws the maps for four grids. The method works in a
band along delta = K / dt0. For case36964 no setting with dt0 of 0.1 or more works, which
explains the paper's unexplained failure at dt0 = 0.1, K = 0.002. For case6024 the
working band lies at a much larger delta than the paper's 0.02.

### Report: test systems and the study summary

- In `report/scripts/figures_from_code.py`, the test-system statistics. For each of the
  34 grids it records the number of buses, PV buses and unknowns, the non-zeros of the
  Jacobian, the flat-start mismatch, and whether NR converges from a flat start and from
  the case file's own start (cached in `report/generated/systems.json`). From these it
  writes the test-systems table of the report (`systems_rows.tex`, which also lists the
  networks each large grid was built from) and a figure comparing the size of every grid
  with its flat-start mismatch (`systems.pdf`).
- The report version of the feasibility maps (`feasibility.pdf`). The grids are cached in
  `feasibility.json`, so the report can be rebuilt without running them again.
- In `report/scripts/figures_from_results.py`, the summary of the improvements study.
  From `improvements_raw.csv` it counts, for every configuration, the runs solved in each
  setting, the wrong roots, and the runs gained and lost against the paper's method,
  splitting the gains into rescues, root fixes and speed-ups. It draws the headline
  figure of runs solved per configuration (`improve_solved.pdf`) and writes the main
  results table of the report (`improve_table.tex`).

---

## Apurbo

Apurbo implemented the core of the paper, the dynamic homotopy, and the framework that
lets the improvements be switched on and off.

### Core code: the homotopy and the integrators

- `dynhomotopy/homotopy/fpv.py` builds the homotopy of Section 3.1:
  G(x, t) = t g(x) + (1 - t) K (x - x0), its derivative in x,
  Gx = t J(x) + (1 - t) K I (equation 23), and in t, Gt = g(x) - K (x - x0)
  (equation 24). The path follows the ODE dx/dt = -Gx^-1 Gt (equation 15), and every
  evaluation of it costs one LU factorization.
- `dynhomotopy/homotopy/integrators.py` has one function per scheme: forward Euler
  (equation 16), backward Euler with one fixed-point iteration (equation 19), RK2
  (equations 20 and 21) and the linearised first step (equations 26 and 27). Apurbo
  showed, and a test checks, that the backward Euler first step and the linear
  approximation give the same point. Their shift K / dt0 on the Jacobian diagonal is
  what tames the ill-conditioning (equation 25).
- `dynhomotopy/homotopy/pathway.py` builds the time points (for example
  {0, 0.005, 0.01, 0.51, 1}), `dynamic.py` follows the path, and
  `dynhomotopy/hybrid.py` joins the homotopy with the NR or FDXB refinement.

### Paper results

- **Figures 1 and 2** (`exp01_tutorial.py`), the 2x2 example: NR alone needs 11 and 10
  iterations for eps = 0.005 and 0.01, exactly as in the paper. Apurbo also showed that
  the paper's eps = 0.05 curve does not fit its own equations.
- **Section 4.2.1** (`exp02_first_step.py`): explicit first steps blow up (mismatch up to
  1e12), while backward Euler brings it down to 65 and 25. This confirms why the paper
  always takes the first step with backward Euler.
- **Figure 4** (`exp07_fig4.py`): on the 109,272-bus grid, two homotopy points are enough
  for NR to converge in 4 iterations. Bus 6 ends at 1.0045 pu and 21.1 degrees, as in the
  paper.

### Improvements: the switch framework, the corrector, adaptive steps and the spectrum of the Jacobian

- `improvements/solve.py` holds `Options` and `solve()`. `solve()` follows the path with
  any step rule and any homotopy function, then refines with NR or FDXB. The other
  members' changes plug into it as switches. With all switches off it gives exactly the
  paper's result, and a test checks this.
- `improvements/path.py`: the basic path follower and the **adaptive** switch, which
  tries to jump straight to t = 1 and halves the step while the mismatch grows. It helps
  little on its own (132 of 170 runs solved, against 127 for the paper's method).
- `improvements/steps.py`: the **corrector**. After every path step, one Newton step on
  G(x, t) = 0 pulls the point back onto the path. This is the most useful change of the
  study: 152 runs solved, none lost, and wrong roots drop from 14 to 2, for about 0.4
  more LUs per case. It also makes the result almost independent of dt0.

`improvements/spectrum.py` estimates the largest eigenvalue modulus of J with the power
method, and the smallest with inverse iteration (one LU of J, then the power method on
J^-1). Shifted inverse iteration at -delta gives the smallest eigenvalue of J + delta I,
the matrix used in the first homotopy step. `experiments/exp12_spectrum.py` runs this on
all 34 grids and checks it against SciPy's ARPACK.

- The eigenvalue ratio lambda_max / lambda_min at the flat start reaches 1.2e8 on the
  109,272-bus grid. The paper's delta = 0.02 lowers it to about 1e6 to 2e6 on the
  ill-conditioned grids, but hardly changes most limit cases, whose smallest eigenvalue
  is already that large. This explains why the method behaves like plain NR on those
  cases.
- Inverse iteration matches ARPACK on 33 of 34 grids. The power method is up to 19% off
  on grids built from copies of one network, where the largest eigenvalues come in
  near-equal clusters.

### Report: the tutorial, the first step, adaptive paths and the spectrum

- `tutorial_geometry.pdf` (`figures_from_code.py`): the 2x2 example drawn in the
  (x1, x2) plane over a map of the mismatch, for K = 0.05 and K = 0.005. It shows the
  exact homotopy path (traced with many small corrected steps), the three-step FE, RK2
  and BE paths with the NR iterations that follow them, the line where the Jacobian is
  singular, and the two roots of the system.
- `adaptive_path.pdf`: every step the adaptive rule tries on case6748 at the paper's
  setting, accepted and rejected, next to the paper's fixed path.
- `first_step.pdf` (`figures_from_results.py`): for every case of the paper, the
  flat-start mismatch next to the mismatch after one FE, RK2 or BE step to t1 = 0.005.
  The explicit steps explode while backward Euler contracts.
- `spectrum.pdf`: for all 34 grids, the smallest eigenvalue modulus of J and of
  J + 0.02 I, and the eigenvalue ratio before and after the shift.

---

## Hasbi

Hasbi wrote the classical solvers that the homotopy hands its result to, and everything
in the project that changes Newton's method itself.

### Core code: the solvers

- `dynhomotopy/solvers/newton.py` implements equation 4, x_(i+1) = x_i - J(x_i)^-1 g(x_i).
  It stops when the largest mismatch is below 1e-8 pu, or after 10 iterations
  (MATPOWER's default), and records the mismatch after every iteration. From a flat start
  it fails on all five large grids, exactly as the paper reports.
- `dynhomotopy/solvers/fdxb.py` implements equations 5 to 8. Each iteration is a P-angle
  update with the constant matrix B' and then a Q-magnitude update with B''. Both are
  factorized only once, which makes an iteration much cheaper than an NR iteration.
  Hasbi found that the paper counts P updates and Q updates separately (for example
  15 + 14 = 29), while MATPOWER returns one count per cycle. With the paper's counting
  every FDXB number in Table 7 matches except one.
- `dynhomotopy/solvers/gsh.py` reconstructs GSH-NR, the static homotopy the paper
  compares against (ref. 12), which the paper only describes in words. Fictitious shunts
  on the Ybus diagonal make the flat start the exact solution at h = 0, the impedances of
  the branches at the slack bus are scaled by delta, and NR solves the problem at
  h = dh, 2 dh, ..., 1. It matches the paper on case18482 (10 iterations) and case36964
  (28 against 27), and fails on case27318, where the paper reports 25.

### Paper results

- **Table 7** (`exp09_table7.py`, with the method list in `setups.py`): the iteration
  counts of all seven methods match the paper on the five large grids in every column
  except two.
- **Table 8** (`exp10_table8.py`): run times relative to NR from the case file. Backward
  Euler followed by FDXB takes 66% to 93% of the NR time (paper: 57% to 85%), and
  followed by NR 102% to 122% (paper: 92% to 128%).

### Improvements: the optimal multiplier, the Newton homotopy and NR and FDXB with our own LU

- `improvements/multiplier.py` scales every NR step by Iwamoto and Tamura's optimal
  multiplier (ref. 10). With a = g(x) and c = g(x + dx), the mismatch along the step is
  modelled as (1 - mu) a + mu^2 c, and mu is the minimiser of its squared norm, the root
  of a cubic. It costs no extra LU. It saves about 11% of the LUs (7.14 against 8.04) but
  does not stop the method from reaching wrong roots. Together with Apurbo's corrector it
  gives the best balance of the study: 152 runs solved with 7.57 LUs.
- `NewtonHomotopy` in `improvements/homotopies.py` replaces the paper's homotopy with the
  classical Newton homotopy g(x) - (1 - t) g(x0). It solves only 75 of 170 runs. It has no
  K I term, which shows that this shift is what makes the paper's method work.

`improvements/scratch.py` runs the paper's method (backward Euler path, then NR or FDXB)
with every linear system solved by Sadman's hand-written LU, and
`experiments/exp15_scratch_lu.py` compares it with SciPy's sparse LU. On the 2x2 example
and the 69- to 2000-bus grids it gives the same iterations and the same solution to
5e-12. It is far slower (424 s against 0.04 s on the 2000-bus grid), because a dense LU
costs O(n^3) and ignores sparsity.

### Report: cost and rescue figures

- `timing.pdf` (`figures_from_results.py`): the Table 8 run times of BE and RK2 followed
  by NR or FDXB, as a percentage of NR from the case file, with the paper's MATLAB value
  marked on each bar.
- `improve_pareto.pdf`: every configuration of the improvements study placed by runs
  solved against LU cost relative to the paper's method, so the robust and cheap ones
  stand out.
- `improve_cases_s1.pdf`: a map of all 34 cases against all configurations at the
  paper's setting, with the LU count of every solved run and marks for wrong roots and
  failures.
- `rescues.pdf` (`figures_from_code.py`): the mismatch history of the paper's method,
  the multiplier, the corrector and both together on three typical runs: a rescue
  (case36964 with dt0 = 0.05), a root fix (case3012wplimit) and a speed-up
  (case14limit).

---

## Sizan

Sizan reproduced most of the paper's convergence tables, built the output that puts our
numbers next to the paper's, and studied the choice of integrator and of delta.

### Core code: tables and timing

- `experiments/common.py` writes every table both as CSV (with proper quoting) and as
  Markdown, formats numbers the way the paper prints them, and holds the timing
  benchmark. Single timings on a laptop varied by up to 4 times, so the benchmark runs
  the methods in turn for several rounds and takes the median.
- `experiments/paper_values.py` holds the numbers copied from the paper, so every table
  prints "ours / paper" side by side. `norm_tables.py` holds the shared code of
  Tables 2 to 5.

### Paper results

- **Tables 2, 3, 4 and 5** (`exp03_tables2to4.py`, `exp04_table5.py`): the backward Euler
  tables match almost digit for digit. For case18482 we get 532, 65, 33, 7.2, 3.7 along
  the path and 0.04, 1.5e-5, 1e-11 in the NR iterations, the same as the paper. Forward
  Euler and RK2 blow up in the same places as in the paper.
- **Figure 3 and Table 6** (`exp06_fig3_table6.py`): five-point paths with a varying third
  point. The paper's one failure (FE with t3 = 0.25 on case36964) is reproduced, and the
  run times for the 109,272-bus grid are within 13 percentage points of the paper in
  every row. RK2 with t3 = 0.20 fails on case36964, where the paper converges.

### Improvements: other step rules and golden-section tuning of delta

In `improvements/steps.py` Sizan added forward Euler, RK2 and RK4 as step rules, and
BE-chord (backward Euler with three fixed-point iterations on one LU, the refinement the
paper suggests after equation 19).

- With fixed steps, RK4 fails almost everywhere (13 of 170 runs), and BE-chord makes the
  start worse (89 of 170), because the extra iterations drift when the step is large.
- Together with Sadman's Richardson step control, Sizan re-tested the paper's claim that
  backward Euler is the best integrator. Under error control it still is the cheapest:
  10.6 LUs per case, against 37 for FE, 48 for RK2 and 87 for RK4.

`improvements/tuning.py` finds delta, then dt0, by golden-section search, minimising the
total number of LU factorizations. `experiments/exp14_tuning.py` tunes every grid, then
fits two rules for choosing delta without a search on the paper's 9 grids and tests them
on the other 25:

| How delta is chosen | Solved, training (9) | Solved, test (25) |
|---|---|---|
| paper, delta = 0.02 | 9 | 19 |
| spectral rule, delta = 1.86 lambda_min | 5 | 20 |
| golden search on each case | 9 | 23 |

The paper's delta suits some of its own grids, but the other ill-conditioned grids need a
delta of 1 to 5. The spectral rule, built on Apurbo's estimate of lambda_min, does not
work, because the best delta is not proportional to the smallest eigenvalue. Searching
each case solves four more test grids, but the search costs about 21 solves.

### Report: the reproduction tables, fidelity and tuning

- In `figures_from_results.py`, the LaTeX rows of Tables 2 to 5 (`table2_rows.tex` to
  `table5_rows.tex`), with our value and the paper's value on alternate lines in the
  report's number format. Tables 4 and 5 appear in the main text and Tables 2 and 3 in
  the appendix. The same script's readers for the saved result tables are also Sizan's.
- `fidelity.pdf`: every numeric cell of Tables 2 to 5 plotted as our value against the
  paper's value, with the number of cells that agree within a factor of 2.
- `tuning.pdf`: the golden-section tuned delta of every grid against its smallest
  eigenvalue, with the paper's delta = 0.02 and the fitted spectral rule drawn in. It
  shows why no single factor of lambda_min fits.

---

## Sadman

Sadman worked on the numerical kernels: the linear algebra, the step-size control, and
the tools that run any part of the project from the command line.

### Core code: linear algebra and results

- `dynhomotopy/linalg.py` wraps SciPy's sparse LU and counts every factorization. The
  paper's cost argument is "one LU per path point", so these counts are the cost measure
  used across the whole project.
- `dynhomotopy/results.py` defines the result types (solver result, homotopy path,
  combined result) that all solvers return.
- `experiments/plots.py` holds the common figure style of the experiment plots.

### Paper results

- **Section 4.4** (`exp05_sec44_pathways.py`): larger first steps (dt0 = 0.05 and 0.1)
  with paths of 3 to 5 points, on all nine cases, with NR and with FDXB. Backward Euler
  works on almost every path, as the paper says. Unlike the paper, FE and RK2 also work on
  some 4-point paths, and backward Euler with NR fails on case36964 with only 3 points
  (FDXB still converges there).
- **Figure 5** (`exp08_fig5.py`): the effect of the time step on the 109,272-bus grid.
  With dt = 0.1 the homotopy hands NR a point of about 1.005 pu and 28 degrees at bus 6,
  and for every time step NR converges to 1.0045 pu and 21.1 degrees, as in the paper.

### Improvements: Richardson step control, scaled homotopy, Gauss elimination and LU from scratch, and the command-line tool

- `follow_path_richardson` in `improvements/path.py` chooses the step size by comparing
  one step of size dt with two steps of size dt / 2. Their difference estimates the local
  error, the kept point is the Richardson extrapolation of the two results, and the next
  step is scaled from the error. A path that needs more than 150 LUs is stopped. It
  solves 151 of 170 runs, and 156 together with Hasbi's multiplier. It is the only change
  that solves case6024 at the paper's own setting, and case7092, case10595 and case12110
  with a weak K, where the paper's method diverges.
- `ScaledHomotopy` in `improvements/homotopies.py` scales the shift K by the Jacobian
  diagonal, in the spirit of the conditioning step of ref. 29. It changes almost nothing
  (129 of 170 runs).
- `improvements/__main__.py` runs any combination of switches on any grid, for example
  `python -m improvements case36964 --corrector --multiplier`, and prints the result next
  to the paper's method. It also runs each study on one grid, for example
  `python -m improvements spectrum case18482`.

- `improvements/gauss.py` implements Gauss elimination and LU factorization with partial
  pivoting, with forward and back substitution. Both agree with NumPy to about 1e-13.
  Hasbi's `scratch.py` uses this LU inside NR and FDXB.
- The comparison at the end of `exp15_scratch_lu.py` shows that factorizing once and
  reusing the LU is 14 times faster than repeating Gauss elimination for 20 right-hand
  sides, which is why FDXB factorizes its matrices only once.

### Report: figure style, sparsity, Section 4.4, Figure 5 and Richardson paths

- `report/scripts/style.py`: the colours, fonts and plot settings shared by every report
  figure, and the helpers that save the figures and the generated LaTeX tables.
- `spy.pdf` (`figures_from_code.py`): the sparsity pattern of the Jacobian for the 2000-,
  18,482- and 109,272-bus grids, with the number of unknowns and non-zeros. It shows why
  a sparse LU is needed at this size.
- `richardson_path.pdf`: the Richardson-controlled paths on case6024 with BE, FE and RK2,
  next to the paper's fixed path: the mismatch along the path and the accepted step
  sizes.
- `sec44_grid.pdf` (`figures_from_results.py`): the Section 4.4 re-run as a grid of paths
  against integration schemes, coloured by whether NR, FDXB, both or neither converge,
  with the NR iteration counts.
- `fig5_tradeoff.pdf`: the homotopy and NR LU factorizations for each time step of
  Figure 5, showing that smaller steps cost more without a better result.

---

## Shared work

Everyone wrote the tests for their own code (24 tests in total):

- `tests/test_core.py`: Jacobian check (Jonayed); pathway, first step and 2x2 example
  (Apurbo); NR and FDXB iteration counts (Hasbi); a regression row of Table 4 (Sizan).
- `tests/test_improvements.py`: all switches off equals the paper's method, corrector and
  adaptive steps, power method and inverse iteration (Apurbo); optimal multiplier and NR
  with our own LU (Hasbi); explicit first step and golden-section search (Sizan); scaled
  homotopy, Gauss elimination, LU, pivoting and Richardson control (Sadman); Newton
  homotopy (Hasbi).
- `tests/test_cli.py`: the command-line tool (Sadman).
