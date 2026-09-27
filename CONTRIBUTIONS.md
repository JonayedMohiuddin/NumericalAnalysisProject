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

As a team we did four things:

1. **Reproduced the paper.** We wrote the whole method from scratch in Python: power
   flow model, solvers, homotopy and integrators. We reran every table and figure of the
   paper on the same public test systems. Most results match the paper digit for digit,
   and every difference is documented.
2. **Checked the answers properly.** Power flow has several solutions, and a small
   mismatch does not show which one was found. We compare every result with a reference
   operating point, and found that some "converged" runs, including some of the
   paper's own method, reach the wrong solution.
3. **Improved the method.** We tried ten changes suggested by the paper and its
   references on a larger set of 34 test grids. Three clearly help: a Newton corrector,
   Iwamoto's optimal multiplier and Richardson step control. Each change is a separate
   switch that can be turned on alone or combined with others.
4. **Delivered our proposal.** Our proposal promised a spectral study of the Jacobian,
   maps of where the method works, golden-section tuning, error-controlled steps, and
   linear solvers written by hand. All of these are implemented, and each has a
   measured result, including the ones that did not work.

The work was split into five parts, one per member. Only the code counts as project
work here (the solver packages, the experiment scripts and the tests), not the
documentation. The parts build on each other in this order:

```
Jonayed: power flow model, data, validation, integration
   |
   +--> Hasbi: classical solvers (NR, FDXB, GSH-NR), Tables 6-8
   |
   +--> Sizan: dynamic homotopy and integrators, Figures 1-5, Tables 2-5
            |
            +--> Apurbo: improvements (corrector, multiplier, Richardson, ...)
                     |
                     +--> Sadman: proposal features (spectrum, maps, tuning, own LU)
```

| Member | Part | Main code | Main results |
|---|---|---|---|
| Jonayed | Power flow model, data, validation and integration | `dynhomotopy/powerflow/`, `datasets.py`, `problem.py`, `linalg.py`, `improvements/cases.py`, `run_all.py`, `experiments/common.py` | flat-start mismatches, base MVA finding, solution check, all table output |
| Hasbi | Classical solvers and the paper's cost comparison | `dynhomotopy/solvers/`, `experiments/setups.py` | Tables 6, 7, 8, Figure 3 |
| Sizan | Dynamic homotopy and the paper's convergence results | `dynhomotopy/homotopy/`, `hybrid.py` | Figures 1, 2, 4, 5, Tables 2-5, Section 4.2.1, Section 4.4 |
| Apurbo | Improvements to the method | `improvements/` (except `cases.py`) | 17-configuration study on 34 grids |
| Sadman | Proposal features | `extensions/` | spectrum, feasibility maps, tuning, own LU |

---

## Part 1: Jonayed, power flow model, data, validation and integration

Jonayed built the foundation that every other part runs on. This
covers turning the published grid files into equations, checking those equations, and
running the whole project from one command.

### Loading the test systems

- `dynhomotopy/datasets.py` lists the nine test systems of the paper (five
  ill-conditioned grids of 18,482 to 109,272 buses and four stressed "limit" cases). It
  downloads them from the two Zenodo records the paper cites the first time they are
  needed.
- `dynhomotopy/powerflow/case.py` reads the MATPOWER `.mat` files. It repeats what
  MATPOWER does before a power flow: it removes isolated buses and switched-off
  generators and branches, renumbers buses 0 to n-1, and classifies every bus as slack,
  PV or PQ. It also accepts both field names used in the files (`mpc` and `mpc_m`).
- `dynhomotopy/powerflow/idx.py` holds the MATPOWER column indices.

### The power flow equations

- `dynhomotopy/powerflow/network.py` builds the bus admittance matrix Ybus with the full
  MATPOWER branch model (series impedance, line charging, tap ratios, phase shifters,
  bus shunts). It also builds the specified injections Sbus and the B' and B'' matrices
  used by the fast decoupled method.
- `dynhomotopy/powerflow/model.py` defines the power flow problem (equations 1 and 2 of
  the paper). The state is x = [angles of PV and PQ buses, magnitudes of PQ buses], the
  mismatch is g(x) = V conj(Ybus V) - Sbus, and the Jacobian (equation 3) is computed
  in sparse form from the derivatives of the bus injections. It also provides the flat
  start and the initial guess stored in each case file.
- `dynhomotopy/problem.py` defines the common interface `g(x)`, `jacobian(x)` that all
  solvers use, and the 2x2 example of Section 3.3 of the paper. Because everything is
  written against this interface, the same solver code runs on the 2x2 example and on
  the 109,272-bus grid.
- `dynhomotopy/linalg.py` wraps SciPy's sparse LU and counts every factorization. The
  paper's cost argument is "one LU per path point", so LU counts are used as the cost
  measure throughout the project.

### Validation

- A finite-difference test checks every column of the Jacobian.
- The flat-start mismatch of all nine cases was compared with the paper (532, 200, 13,
  54, 0.012 and 0.006 pu). Seven matched at once. The two small limit cases came out
  exactly 10 times too large. Jonayed traced this to their files storing a 10 MVA base,
  while the paper's numbers for Tables 2-5 only come out on a 100 MVA base. The paper's
  own numbers for these two cases do not fit a single base, so this is documented and
  the base used is stated in every table.
- **Solution check.** Power flow equations have several solutions. Jonayed added
  `reference_voltage` and `on_reference` (`improvements/cases.py`), which compute the
  operating point NR reaches from the case file's own starting point and test whether
  a result is within 1e-4 pu of it. The check showed that 14 of the 141 runs where the
  paper's method converged on our wider test bed had reached a different root. From
  then on, a run only counts as solved if it reaches the reference. Every entry of
  Table 7 was confirmed to reach the correct solution.

### Integration

- `run_all.py` runs all 15 experiments in order, or any selection of them.
- `experiments/common.py` holds the shared helpers: writing every table both as CSV
  (with proper quoting) and as Markdown (with `|` escaped), number formatting close to
  the paper's, and saving figures.
- `experiments/paper_values.py` holds the numbers copied from the paper, so every table
  prints "ours / paper" side by side.
- The tests in `tests/test_core.py` for the model, the solvers' reference iteration
  counts and the table output.

---

## Part 2: Hasbi, classical solvers and the paper's cost comparison

Hasbi wrote the solvers that the homotopy hands its result to, and the method the paper
compares against. Hasbi also produced the paper's iteration and CPU time tables.

### Newton-Raphson

`dynhomotopy/solvers/newton.py` implements equation 4 of the paper:
x_(i+1) = x_i - J(x_i)^-1 g(x_i). It stops when the largest mismatch is below 1e-8 pu,
or after 10 iterations (MATPOWER's default). It also records the mismatch after every
iteration, since the paper's tables print these. From a flat start it fails on all five
large grids, exactly as the paper reports.

### Fast decoupled load flow (FDXB)

`dynhomotopy/solvers/fdxb.py` implements equations 5 to 8. The coupling between angles
and magnitudes is ignored, so each iteration is a P-angle update with the constant matrix
B' and then a Q-magnitude update with B''. Both matrices are factorized only once, which
makes an FDXB iteration much cheaper than an NR iteration. Hasbi found that the paper
counts P updates and Q updates separately (for example 15 + 14 = 29), unlike MATPOWER,
which returns one count per cycle. With that counting every FDXB number in Table 7 matches
except one.

### GSH-NR

The paper compares against a static homotopy from an earlier paper (ref. 12), but only
describes it in words. Hasbi reconstructed it in `dynhomotopy/solvers/gsh.py`:

- fictitious shunts are added to the diagonal of Ybus so that the flat start is the exact
  solution of the modified network at h = 0, and they are removed as h goes to 1;
- the impedances of the branches connected to the slack bus are scaled by delta at h = 0
  and restored at h = 1;
- NR solves the problem at h = dh, 2 dh, ..., 1, each time starting from the previous
  solution, and the iterations are added up.

The reconstruction matches the paper on case18482 (10 iterations) and case36964 (28
against 27), and reduces to plain NR when dh = 1 and delta = 1, as the paper's limit-case
rows require. It fails on case27318, where the paper reports 25 iterations, which shows
that the original method differs in details the paper does not give.

### Tables and timing

- `experiments/setups.py` defines the seven methods of Table 7 (NR from the case file,
  NR from a flat start, GSH-NR, BE or RK2 followed by NR or FDXB).
- **Table 7** (`exp09_table7.py`): the iteration counts match the paper on the five large
  grids in every column except two.
- **Figure 3 and Table 6** (`exp06_fig3_table6.py`) and **Table 8** (`exp10_table8.py`):
  run times relative to NR from the case file. Single timings on a laptop varied by up to
  4 times, so Hasbi wrote an interleaved benchmark: the methods run in turn, 7 rounds,
  and the median is taken. With it, BE followed by FDXB takes 66-93% of the NR time
  (paper: 57-85%), BE followed by NR 102-122% (paper: 92-128%), and Table 6 for the
  109,272-bus grid is within 13 percentage points of the paper in every row.

---

## Part 3: Sizan, the dynamic homotopy and the paper's convergence results

Sizan implemented the core of the paper and reproduced its convergence results.

### The homotopy and the integrators

- `dynhomotopy/homotopy/fpv.py` builds the homotopy of Section 3.1:
  G(x, t) = t g(x) + (1 - t) K (x - x0), its derivative in x,
  Gx = t J(x) + (1 - t) K I (equation 23), and in t, Gt = g(x) - K (x - x0)
  (equation 24). The path then follows the ODE dx/dt = -Gx^-1 Gt (equation 15), and
  every evaluation of it costs one LU factorization.
- `dynhomotopy/homotopy/integrators.py` has one function per scheme: forward Euler
  (equation 16), backward Euler solved with one fixed-point iteration (equation 19),
  RK2 (equations 20 and 21), and the linearised first step (equations 26 and 27).
  Sizan showed, and a test checks, that the backward Euler first step and the linear
  approximation give the same point. Their shift K/dt0 on the Jacobian diagonal is what
  tames the ill-conditioning (equation 25).
- `dynhomotopy/homotopy/pathway.py` builds the time points (for example
  {0, 0.005, 0.01, 0.51, 1}), `dynhomotopy/homotopy/dynamic.py` follows the path, and
  `dynhomotopy/hybrid.py` joins the homotopy with the NR or FDXB refinement.

### Results reproduced

- **Figures 1 and 2** (`exp01_tutorial.py`), the 2x2 example: NR alone needs 11 and 10
  iterations for eps = 0.005 and 0.01, exactly as in the paper. Sizan also showed that
  the paper's eps = 0.05 curve does not fit its own equations.
- **Section 4.2.1** (`exp02_first_step.py`): explicit first steps blow up (mismatch up
  to 1e12), while backward Euler gives 65 and 25, confirming why the paper always takes
  the first step with backward Euler.
- **Tables 2, 3, 4 and 5** (`exp03_tables2to4.py`, `exp04_table5.py`): the backward Euler
  tables match almost digit for digit. For case18482 we get 532, 65, 33, 7.2, 3.7 along
  the path and 0.04, 1.5e-5, 1e-11 in the NR iterations, the same as the paper.
  Forward Euler and RK2 blow up in the same places as in the paper.
- **Section 4.4** (`exp05_sec44_pathways.py`): larger first steps with paths of 3 to 5
  points.
- **Figures 4 and 5** (`exp07_fig4.py`, `exp08_fig5.py`): voltage paths on the
  109,272-bus grid. Bus 6 ends at 1.0045 pu and 21.1 degrees for every time step, and NR
  needs 4 iterations after only two homotopy points.

---

## Part 4: Apurbo, improvements to the method

Apurbo built the `improvements` package, which changes the paper's method without touching
the paper's code. Every change is a switch in `Options`. With all switches off, the
result is exactly the paper's method, and a test checks this.

### The changes

| Switch | Idea | Source |
|---|---|---|
| `corrector` | after every path step, one Newton step on G(x, t) = 0 pulls the point back onto the path | predictor-corrector, refs. 36, 37 |
| `multiplier` | every NR step is scaled by the step length that minimises a quadratic model of the mismatch | Iwamoto and Tamura, ref. 10 |
| `richardson` | the step size is chosen by comparing one step with two half steps (Richardson error estimate) | our proposal |
| `adaptive` | try to jump straight to t = 1, halve the step while the mismatch grows | paper's future work |
| `step` | FE, RK2, RK4 or BE with extra fixed-point iterations instead of plain BE | paper's future work |
| `homotopy` | a Jacobian-scaled shift, or the classical Newton homotopy | refs. 29, 36 |

The code is split by kind of change: `multiplier.py` (the optimal multiplier and NR with
it), `steps.py` (step rules and the corrector), `homotopies.py` (other homotopy
functions), `path.py` (fixed, adaptive and Richardson-controlled time points) and
`solve.py` (the switches). Each change can be run on any case from the command line, for
example `python -m improvements case36964 --corrector --multiplier`.

### The study

`experiments/exp11_improvements.py` runs 17 configurations on 34 grids (the paper's 9
plus 25 more from the same records) under 5 settings of dt0 and K: 170 runs per
configuration. A run counts only if NR converges within 10 iterations and reaches the
reference solution. Cost is counted in LU factorizations.

| Configuration | Solved (of 170) | Wrong root | Mean LUs |
|---|---|---|---|
| paper's method | 127 | 14 | 8.04 |
| corrector | 152 | 2 | 8.43 |
| multiplier + corrector | 152 | 3 | 7.57 |
| multiplier | 143 | 12 | 7.14 |
| multiplier + Richardson | 156 | 8 | 9.75 |

What Apurbo found:

- The corrector is the most useful change: 25 more runs solved, none lost, and wrong
  roots drop from 14 to 2.
- The multiplier saves about 11% of the LUs but does not prevent wrong roots. Together
  with the corrector it gives the best balance (multiplier + corrector).
- Richardson step control is the only change that solves case6024 at the paper's own
  setting, and case7092, case10595 and case12110 with a weak K, where the paper's
  method diverges.
- Under error control backward Euler is still the cheapest integrator (10.6 LUs against
  37 for FE, 48 for RK2 and 87 for RK4), which confirms the paper's main claim.
- Extra fixed-point iterations (the paper's own suggestion), fixed-step RK4 and the
  Newton homotopy make things worse. The Newton homotopy has no K I term, which shows
  that this shift is what makes the method work.

An earlier version of this study counted any small mismatch as success. After the
solution check (Part 1) was added, the study was rerun, and the recommendation changed
from "multiplier + adaptive" to "multiplier + corrector".

---

## Part 5: Sadman, the proposal features

Sadman implemented the numerical parts promised in our proposal deck, in the `extensions`
package, and measured what each one gives.

### Spectrum of the Jacobian

`extensions/spectrum.py` estimates the largest eigenvalue modulus of J with the power
method, and the smallest with inverse iteration: one LU of J, then the power method on
J^-1. Shifted inverse iteration at -delta gives the smallest eigenvalue of the shifted
matrix J + delta I used in the first homotopy step. `exp12_spectrum.py` runs this on all
34 grids and checks it against SciPy's ARPACK.

- The eigenvalue ratio lambda_max / lambda_min at the flat start reaches 1.2e8 on the
  109,272-bus grid. The paper's delta = 0.02 lowers it to 1e6-2e6 on the ill-conditioned
  grids but hardly changes most limit cases, whose smallest eigenvalue is already that
  large. This explains why the method behaves like plain NR on those cases.
- Inverse iteration matches ARPACK on 33 of 34 grids. The power method is up to 19% off
  on grids built from copies of one network, where the largest eigenvalues come in
  near-equal clusters.

### Feasible region maps

`extensions/feasibility.py` and `exp13_feasibility.py` run the paper's method on a grid
of K and dt0 and mark each cell as solved, wrong root or failed. The method works in a
band along delta = K / dt0. For case36964 no setting with dt0 of 0.1 or more works, which
explains the paper's unexplained failure at dt0 = 0.1, K = 0.002. For case6024 the
working band lies at a much larger delta than the paper's 0.02.

### Golden-section tuning

`extensions/tuning.py` finds delta, then dt0, by golden-section search, minimising the
total LU factorizations. The paper's delta suits some of its own grids, but the other
ill-conditioned grids need delta between 1 and 5. Two rules for choosing delta without a
search were fitted on the paper's 9 grids and tested on the other 25 (`exp14_tuning.py`):

| How delta is chosen | Solved, training (9) | Solved, test (25) |
|---|---|---|
| paper, delta = 0.02 | 9 | 19 |
| spectral rule, delta = 1.86 lambda_min | 5 | 20 |
| golden search on each case | 9 | 23 |

The spectral rule does not work, because the best delta is not proportional to the
smallest eigenvalue. Searching each case solves four more test grids, but the search
costs about 21 solves.

### Linear solvers written from scratch

`extensions/gauss.py` implements Gauss elimination and LU factorization with partial
pivoting, and `extensions/scratch.py` runs the paper's method (backward Euler path, then
NR or FDXB) with them. On the 2x2 example and on the 69- to 2000-bus grids, they give the
same iterations as SciPy's sparse LU and the same solution to 5e-12. They are far slower
(424 s against 0.04 s on the 2000-bus grid), because a dense LU costs O(n^3) and ignores
sparsity. Reusing one LU is 14 times faster than repeating Gauss elimination for 20
right-hand sides, which is why FDXB factorizes its matrices only once.

Each feature runs on its own, for example `python -m extensions spectrum case18482`.
Richardson step control, the fifth feature of the proposal, was built together with
Apurbo as a switch of the improvements package.

---

## Shared work

Everyone wrote tests for their own part: `tests/test_core.py` (Jonayed, Hasbi, Sizan),
`tests/test_improvements.py` (Apurbo), `tests/test_extensions.py` (Sadman) and
`tests/test_cli.py` (Apurbo, Sadman), 24 tests in total.
