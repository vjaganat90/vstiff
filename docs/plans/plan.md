# vstiff plan: a verified, property-tested stiff solver

Status (2026-10-03): accepted. D15 to D18 and D20 follow their recommendations, and D19 is
decided as section 9 states. This is the top-level plan. It sets the order of work for two
accepted documents:

- the [roadmap](roadmap.md), with milestones M1 to M13 and decisions D1 to D14;
- the [formal verification report](formal-verification.md), called the FV report here, with
  theorems T1 to T12.

It also adds what they leave open:

- a definition of done;
- a property-testing program and a robustness program;
- global error estimation;
- a comparison of correctness with other solvers.

Section 8 lists where it changes either document. Claims from memory of the literature are
marked **[M]**. Facts about the local toolchain were checked on 2026-10-03.

---

## 0. Summary

vstiff today is a variable-step BDF2 integrator with:

- a halve/double controller;
- full Newton with a finite-difference Jacobian;
- dense Gaussian elimination.

47 corpus lines and 4 soak lines pin its behaviour. One lemma about its coefficients is checked
in Rocq, in `verif/` on the verification branch.

The target is a general stiff solver whose correctness is shown in three independent ways:

- **Tested.** Every promise in a public interface is covered by one of: an expect line, a
  property run on thousands of generated problems, or a theorem. No non-equivalent mutant of the
  core modules survives.
- **Proved.** Rocq checks these, against mirrors that the build keeps in sync with the code:
  - the method coefficients;
  - the zero-stability of every step sequence the controller can produce;
  - the controller's invariants and termination;
  - the linear solver;
  - the conservation of linear invariants.
- **Measured.** A reproducible scorecard runs vstiff, scipy and SUNDIALS on the same problems. It
  reports cost, and how often each solver returns a wrong answer without saying so.

The work runs in four tracks: solver, properties, proofs and benchmarks. They proceed in eight
waves. Each wave leaves the solver green and better, so the plan can stop after any of them. The
estimate is 121 to 184 person-days of engineering, plus 12 to 30 weeks of proof work. These
figures are judgement, not data.

No solver can be right on every input, because it samples a right-hand side at finitely many
points. The plan therefore promises:

- machine-checked properties of the method and its control logic;
- measured behaviour on published suites;
- error estimates whose quality is tested;
- rigorous bounds only where section 5's research step shows they can be reached.

On speed it promises no more than the roadmap's criterion 5: within 2x of scipy's BDF at equal
error.

---

## 1. Definition of done

The plan is done when every item below holds. Each one is checkable.

### 1.1 Capability

This is the roadmap's "general" (its section 1.3), plus four items:

| Capability | Milestone |
|---|---|
| rtol and atol vectors, WRMS norm | M2 |
| Variable-order BDF 1 to 5; Radau IIA of order 5 | M5, M10 |
| Finite-difference Jacobians (dense, grouped for banded, coloured for sparse), or a user Jacobian | M6, M6b |
| Dense, banded and sparse LU; Krylov only if a target problem needs it | M3, M6, M6b, M11 |
| Index-1 DAEs with a mass matrix | M9 |
| Dense output, output at requested times, events | M7 |
| New: integration backward in time (`t_end < t0`) | M2 (D17) |
| New: on failure, the last accepted point, the statistics and the reason | M2 (D16) |
| New: the solver state is a value, so a run stopped after any step resumes with the same steps | M7 (with D7's `Seq.t`) |
| New: an optional global error estimate, and a certified mode that tightens the tolerance until the estimate meets the request | G1, G2 (section 5) |

### 1.2 Robustness

- **R1.** The roadmap's robustness criteria 1 to 4 (its section 1.3) hold on every corpus
  problem, the adversarial suite of section 4 included. They hold at rtol 1e-3 to 1e-8, and at
  1e-10 for problems with n <= 20.
- **R2.** Every entry point returns `Ok` or a typed `Error`, for every input:
  - Invalid options are programmer errors and raise `Invalid_argument`. So is a right-hand side
    that returns a vector of the wrong length.
  - An exception raised by the user's right-hand side propagates unchanged (D15).
  - The solver raises nothing else.
- **R3.** No silent wrong answer under fault injection. Such a run's right-hand side returns NaN
  or infinity at seeded random calls, or adds noise below the tolerance. The run must return
  either a typed `Error`, or an `Ok` within criterion 2's bound of the fault-free run. The nightly
  campaign (3.5) must find zero violations.
- **R4.** Bounded work. Every run ends within a number of right-hand-side calls fixed by its
  options (steps, rejections, budget). A property checks this, and T4 proves it for the driver
  model.
- **R5.** Determinism within a build. The same binary, given the same inputs, returns the same
  outputs, across processes and across OCaml 5 domains running solves concurrently. This rules
  out hidden state, randomness and races. Platforms, compilers and refactors may change the last
  bits, and every check holds on each of them (D19).
- **R6.** Scale. These finish: a banded heat problem with n = 1e5, the Brusselator with n = 2000,
  and a sparse circuit or network problem with n from 1e3 to 1e4. Their times are reported, not
  gated.

### 1.3 Property coverage

- **P1.** `docs/contracts.md` lists every promise made in a public `.mli`. For each one it lists
  the checks that cover it: expect lines, properties and theorems. A test fails when a listed
  check does not exist. Each wave audits the interfaces for promises the list is missing.
- **P2.** CI runs every property with a fixed seed, within a time budget of about a minute. The
  nightly run uses a fresh seed and a hundred times the cases. A failure prints its seed and its
  shrunk case. Its fix lands with a pinned line that reproduces it.
- **P3.** Mutation testing (3.4) leaves no surviving non-equivalent mutant in the kernel, the
  methods, the controllers or the drivers. Equivalent mutants are listed, each with the reason.
- **P4.** A statistical property, such as tolerance proportionality or estimate quality, prints
  its exact count. Its threshold is a correctness bound. The count may move with an algorithm
  (re-pinned as D13 says), but the threshold never loosens.

### 1.4 Verification

These theorems are machine-checked in CI, on Rocq mirrors. The FV report's tripwire and
differential pins keep the mirrors in sync with the code (its section 4.5). Theorem IDs T1 to T12
are the FV report's; T13 to T17 are new.

| Group | Theorems |
|---|---|
| Method | BDF coefficients of every order 1 to 5 are exact on polynomials of degree <= k, through any distinct nodes (T1, generalized as T13). Radau IIA's order conditions and algebraic stability hold, with its exact coefficients (T16, T17). |
| Stability | Variable-step BDF2 is zero-stable for ratios below 1 + sqrt 2 (T2). For every order, a certificate shows that every ratio sequence the controller can produce is zero-stable (T14). BDF1 and BDF2 are A-stable, and the A(alpha) angles of BDF3 to BDF5 have proved lower bounds (T5, T15). |
| Control | The controller never attempts a ratio outside its order's certified set, nor a step below the clock's floor. The driver terminates within its attempt bound. (T3, restated for the M2 controller; T4.) |
| Kernel | LU is correct over an ordered field (T6a). With modified Newton and finite-difference Jacobians, the stage solve conserves linear invariants exactly (T9). |
| Floats | Every accepted state is finite (T10). The computed coefficients have rounding bounds (T1f). |

Each theorem comes with three things:

- a note that maps it to the OCaml lines it covers;
- a property that checks its conclusion on the compiled code;
- its written-out assumptions: properties of the user's right-hand side, IEEE arithmetic and the
  compiler.

### 1.5 The scorecard

The scorecard (section 7) is a published, reproducible comparison of vstiff with other solvers:

- scipy's BDF, Radau and LSODA;
- SUNDIALS's CVODE and IDA (D18).

It runs on:

- the IVP Test Set;
- the roadmap's corpus (its section 4.1);
- ten thousand generated stiff problems with exact solutions;
- the adversarial suite;
- the fault-injection runs.

vstiff is built to earn these claims, each reported next to the other solvers' numbers:

1. **No silent wrong answer on any suite.** Every `Ok` is within criterion 2's bound, and every
   failure is typed.
2. **Certified mode meets its tolerance.** The achieved error is at most the requested tolerance
   on at least 99 % of (problem, tolerance) pairs.
3. **The reliability profile is at least as high as every compared solver's at K = 100.** The
   profile is the share of runs with error/tol <= K.
4. **The guarantees of 1.4 are machine-checked.** None of the compared solvers ships such proofs
   **[M]**.

At equal achieved error, the cost is within 2x of scipy's BDF (the roadmap's criterion 5).

---

## 2. Tracks and gates

| Track | Scope | Branch | Gate |
|---|---|---|---|
| S, solver | The roadmap's milestones, with M6b added | `feat/solver` | The milestone's corpus lines and properties; mutation on the modules it touches; the work-precision gate (roadmap 4.3) |
| P, properties and robustness | Harness, generators, properties, mutation tool, adversarial suite, fault injection, nightly campaigns | `feat/solver` | Its runs are part of CI |
| V, proofs | Theorems, mirrors, tripwire, differential pins | `feat/verification`, rebased onto `feat/solver` after each solver milestone | `dune build --root verif` in CI; the tripwire is green; the theorem's shadow property exists |
| B, benchmarks | `bench/` and the scorecard | `feat/solver` (`bench/compare/` is outside the build) | Regenerated from scripts; references follow roadmap 4.2 |

Five rules hold across tracks:

1. A defect can be found anywhere: by a property, a nightly run, a proof attempt or the
   scorecard. It lands together with its fix and a pinned line that would have shown it. The
   commit message records the old output.
2. Proofs trail code by one wave. A theorem is stated against a module's mirror once that
   module's contract has settled. No proof starts about a module that the next milestone
   replaces.
3. Pins move only as D13 says, and correctness bounds never loosen.
   [`test/refs.ml`](../../test/refs.ml) is never edited; new references go into new modules.
4. Checks are about correctness, not bits (D19). An expect line prints a verdict, a typed
   outcome, or a value to the digits its bound makes meaningful, so the same files pass on every
   platform of the CI matrix. A change that moves a line either changed behaviour, or exposed a
   line printed more precisely than its check needs; the commit message says which.
5. Documentation is updated once per wave.

---

## 3. Property-testing program

### 3.1 Harness

The harness uses the standard library only (D2) and lives in `test/prop/`.

- **Generators.** A generator draws from a `Random.State.t`. It returns its value together with a
  lazy tree (`Seq.t`) of smaller candidates. Shrinking therefore comes with the generator, and
  nobody writes a shrinker per type.
- **Properties.** A property is a module, like the soak test's `Case`
  ([`test/soak.ml`](../../test/soak.ml)):

  ```ocaml
  module type Property = sig
    type t

    val name : string
    val gen : t Gen.t
    val show : t -> string
    val holds : t -> bool
  end
  ```

  `Prop.check (module P : Property)` runs it. A statistical property also has a threshold and
  reports its count.
- **Seeds.** Each property draws from its own stream, split from the run seed and its name. Adding
  a property therefore changes no other property's cases. CI fixes the run seed; the nightly job
  reads it from the environment.
- **Effects stay where they are.** The harness returns values and
  [`Report`](../../test/report.ml) prints them. [`Guard`](../../test/guard.ml) bounds each case's
  right-hand-side calls, so a hanging case fails instead of stalling the run.
- **Output.** A passing property prints `prop <name>: ok <count>`, and no generated value. Its
  line survives a change to the generators' internals or to the standard library's `Random`.

### 3.2 Generators and their oracles

| Generator | Oracle |
|---|---|
| Floats over all magnitudes, and the special values (+-0, subnormals, +-infinity, NaN, `max_float`, `epsilon`) | The contracts: `Ok` or a typed `Error`, and finiteness |
| Matrices with prescribed singular values (random Householder reflections around a chosen diagonal); banded, sparse, and exactly singular ones | Residual and backward-error bounds; singularity detected; dense, banded and sparse LU agree |
| Linear ODEs y' = A y with A = Q D Q^T. D is block diagonal: real eigenvalues up to a stiffness ratio of 1e10, and 2 x 2 rotation blocks for complex pairs placed in chosen sectors. Non-normal variants use a well-conditioned V in place of Q | The exact solution, in closed form per block |
| Manufactured stiff problems f(t, y) = y*'(t) + L (y - y*(t)) + N(y - y*(t)). Here y* is a random smooth function, L = -Q diag(lambda) Q^T, and N is quadratic with N(0) = 0 | y* itself |
| Problems with conserved quantities: random mass-action reaction networks, and fields with c . f = 0 by construction | The conserved quantities, and positivity for the networks |
| Index-1 DAEs with a manufactured solution and a nonsingular algebraic Jacobian (M9) | The solution and the constraint residual |
| Problems with known event times (M7) | The event times |
| Adversarial wrappers around any of the above: NaN or infinity at seeded calls or in regions, exceptions, noise below the tolerance, discontinuities, finite-time blow-up | The fault-free run, and the contracts R2 and R3 |
| Step-ratio sequences inside the controller's set and at its edge | T2 and T14 |

### 3.3 Properties by layer

The contract map lists every property. These examples show the range.

- **Kernel.**
  - LU's residual stays within the backward-error bound of partial pivoting, with the observed
    growth factor.
  - `solve (P A) (P b)` agrees with `solve A b` within that bound.
  - Banded and sparse LU agree with dense LU within that bound.
  - `Jac.forward` is exact on affine maps up to its rounding bound. It keeps the orientation
    J_ij = df_i/dy_j; an asymmetric map catches a transpose.
  - Coloured differences agree with plain ones within the finite-difference rounding bound.
  - `Newton.solve` never returns `Ok` with a non-finite vector. An `Ok` passes its own stopping
    test when that test is evaluated again.
- **Methods.**
  - The coefficients are exact on polynomials of degree <= k at random nodes. The oracle is
    double-double arithmetic built on `Float.fma`.
  - Fixed-step runs of each order show their order on at least 95 % of manufactured problems:
    halving the step divides the error by between 0.8 x 2^k and 1.25 x 2^k.
  - BDF and Radau agree within the sum of their tolerances (M10).
- **Controllers.**
  - Every proposal lies within the clamp and above the clock's floor.
  - Every attempted ratio lies inside its order's certified set (the shadow of T3 and T14).
  - After a rejection, the next proposal is smaller.
  - The attempt count stays within T4's bound.
- **Drivers.**
  - Requested output times are hit exactly.
  - Dense output reproduces the accepted states to rounding error.
  - These transformations act on the solution in the same way, within the tolerance:
    - permuting the components;
    - scaling y and atol by a power of two;
    - shifting t0.
  - A nonstiff problem integrated forward and then back returns to its start within criterion
    2's bound.
  - Tolerance proportionality holds on at least 99 % of manufactured problems: achieved error
    over tolerance stays within [1e-2, 1e2], for rtol from 1e-3 to 1e-8.
- **Robustness.** R2 to R5 are properties: fault injection, budgets, partial results, and the
  same results from concurrently running domains as from a sequential run.
- **Global error.**
  - The estimate is within a factor 10 of the true error on at least 95 % of manufactured
    problems.
  - Certified mode meets the tolerance on at least 99 % of them.

### 3.4 Mutation testing

`mutaml` does not install on OCaml 5.5 (roadmap 4.5). This plan replaces the roadmap's sed
script with a small executable built on `compiler-libs`. That library ships with the compiler, and
the FV report's planned tripwire uses it as well.

The executable parses each module and enumerates its mutants:

- arithmetic and comparison operators swapped;
- conditions negated;
- constants perturbed;
- integer offsets of one;
- arguments of the same type swapped;
- one branch's body replacing another's.

Each mutant is built in a scratch copy and run against the corpus, the soak test and the
properties, under `Guard`'s budget. Its outcome is one of: killed (by a test, or by the budget),
survived, or does not compile. A survivor is resolved in one of two ways:

- a new line or property kills it;
- it joins the list of equivalent mutants, with the reason.

The tool runs on the touched modules for every milestone, on changed modules nightly, and on
everything once per wave.

### 3.5 Nightly campaigns

A scheduled CI job runs outside the per-push budget. It covers:

- every property, with a fresh seed and a hundred times the cases;
- differential runs: generated problems through every combination of method, controller and
  linear solver, reporting any disagreement beyond the combined tolerances;
- the fault-injection campaign of R3;
- mutation on the modules changed since the last run;
- the large problems and the work-precision gate (roadmap D11).

---

## 4. Robustness program

The adversarial suite extends the failure modes that the roadmap pins (its section 4.4). Each
case is a corpus line that prints `Ok` or the typed `Error`:

- finite-time blow-up: y' = y^2, from y(0) = 1 to t = 2. The run ends in an `Error` near t = 1
  with the last accepted point, never in `Ok`;
- non-Lipschitz y' = -sqrt |y| through 0; discontinuous forcing; a right-hand side that branches
  on the state;
- a Newton matrix that is exactly singular at the proposed step (h beta lambda = 1), and a nearly
  singular Jacobian;
- extreme scales: components near 1e-300 and 1e300, subnormal states, t0 = 1e15, spans across
  t = 0, negative times, backward integration (D17);
- stiffness from 1 to 1e12 (van der Pol up to mu = 1e6), and eigenvalues near the imaginary axis
  (the roadmap's 60 and 80 degree canaries);
- NaN walls (logarithm and square-root domains), and noise below the tolerance. The run either
  finishes with bounded rejections or fails with a typed error.

The other robustness checks:

- **Fault injection** wraps the right-hand side in `Guard`, the test module allowed to raise. At
  seeded calls the wrapper returns NaN or infinity, raises, or adds noise. Its properties are R2
  and R3.
- **Partial results and budgets** (D16) arrive with M2's contract change, so every later failure
  reports where it stopped.
- **Hidden shared state.** A property runs the same solves in several domains at once and
  compares the bits with a sequential run.
- **Options.** Every new option gets its `Check` rule and a pinned line, as d709e81 did for `tol`
  and for the right-hand side's length.

---

## 5. Global error estimation

Local error control does not bound the global error, and none of the solvers compared in section 7
reports one **[M]**. This plan adds global error in three steps; the third is conditional.

- **G1, certified mode** (wave 4; it needs dense output, M7).
  - Solve at the requested tolerance and again at a tighter one.
  - Estimate the error of the coarser run from the difference on a common grid of output times.
  - Return the finer run, with that estimate as a conservative bound on its error.
  - If the estimate exceeds the requested global tolerance, tighten and repeat a bounded number of
    times. Then return `Error` with the best run and its estimate.

  After M10 the second run can use Radau IIA instead, whose errors arise differently from BDF's.
  The cost is two to three solves.
- **G2, a running estimate** (wave 6). Integrate the linearized error equation
  e' = J e + delta(t) alongside the solution. Here delta is the defect of the dense output. The
  equation reuses the Newton matrix already factored for the step, at the cost of one extra solve
  and a few right-hand-side calls per step. The target is under 40 % extra cost **[M: Lang and
  Verwer, SIAM J. Sci. Comput. 29 (2007)]**.
- **G3, a bound** (wave 7, research, with a go/no-go after a five-day prototype).
  - For dissipative problems, the bound is
    ||e(t)|| <= int_t0^t exp(int_s^t mu) ||delta(s)|| ds.
  - delta is the defect, and mu bounds the logarithmic norm of the Jacobian over a tube around the
    computed solution. mu is cheap to compute in a weighted infinity norm.
  - For small systems the bound becomes rigorous by evaluating the defect and the Jacobian with
    interval arithmetic. Outward rounding uses `Float.pred` and `Float.succ`, so this needs only
    the standard library.
  - The decision is go if the bound is within 100x of the true error on HIRES and Robertson.

The properties of 3.3 test the estimates of G1 and G2 against exact solutions, and the scorecard
reports them.

---

## 6. Verification, re-sequenced

| Step | Theorems | Starts after | Notes |
|---|---|---|---|
| V0 | T1 with the cubic residual; T2; the tripwire for `Bdf2.coeffs` and `Clock`; the proofs job in CI | now | T3 waits: M2 replaces the controller contract, so a proof about `Halving` would be discarded |
| V1 | T3 for the M2 controller and driver (ratio set, floor, `dt_max`); T4 | M2's contract change | Covers both directions of time (D17) from the start |
| V2 | T6a for M3's LU algorithm; T9 with modified Newton | M3 | |
| V3 | T13: exactness of BDF1 to BDF5 on arbitrary nodes; T14: zero-stability certificates; T15: lower bounds on the A(alpha) angles | M5 | T14 carries research risk; see below |
| V4 | T16: Radau IIA's order conditions B(5), C(3), D(2), in Q(sqrt 6); T17: algebraic stability; coefficient rounding | M10 | |
| V5 | T10; T1f; T3 in binary64; T7 if wanted | M8 | |

**T14** is the step with research risk.

- At order k, the step-difference recurrence depends on the last k - 1 step ratios.
- A quadratic Lyapunov function for it is found numerically.
- Rocq checks that the function decreases over the box of ratios the controller allows. The tool
  is interval bisection with CoqInterval. ValidSDP, the obvious alternative, requires MathComp
  below 2.6 (checked).
- The known variable-step limits for BDF3 and above are small (FV report 4.8). If no certificate
  covers the set the controller uses, the controller is restricted to one that can be certified.
  For instance, it can hold the step for k steps after a change at order k. That restriction is a
  D13 re-pin, and its cost shows in the work-precision gate.

**Differential pins** (FV report 4.5, item 4) run in CI throughout. Rocq's primitive-float mirror
computes coefficient tables and controller decision traces, and CI compares them with the OCaml
output within stated rounding bounds, not bit for bit. Rocq's primitive floats have no fused
multiply-add, and the OCaml code may use it (D19). Traces are compared only on inputs that stay
clear of each decision threshold by more than those bounds.

**Not planned:**

- extracting shipped code from Rocq (FV option B);
- end-to-end binary64 correctness (FV section 3: impractical);
- the stability of variable-step BDF2 at ratio 2 for dissipative problems (T12: research).

---

## 7. The scorecard

- **Where.** `bench/` holds the vstiff driver (roadmap M1). `bench/compare/` holds the other
  solvers' drivers. It is outside the dune build, as `verif/` is, and pins its versions in a lock
  file.
- **Solvers.**
  - scipy's BDF, Radau and LSODA: scipy 1.18.1 is installed (checked).
  - CVODE and IDA, through the OCaml bindings `sundialsml` (D18). The opam release 6.1.1p1 needs
    the SUNDIALS C library, which is not installed (checked). Whether the bindings build against
    the current SUNDIALS release has not been checked.
  - Julia's solvers and Fortran RADAU5 are not planned.
- **Transcription check.** Every problem exists in OCaml and in Python.
  - Generated problems are exported as their parameters, and each side builds the same right-hand
    side from them.
  - Before any comparison, both right-hand sides are evaluated at a hundred seeded points and must
    agree to 1e-15 relative.
  - Without this check, a mistranscribed problem would look like a solver failure.
- **Fairness.** Each solver runs twice: at its defaults, and at its documented best settings.
  - Every solver gets the same tolerance definition (mixed WRMS).
  - Every solver gets the same Jacobian information: finite differences for all, or the analytic
    Jacobian for all.
  - CPU time is reported, not compared.
- **Metrics.** Each run is one problem, tolerance and solver, and reports:
  - its status;
  - its achieved error, as scd and mescd (roadmap 4.2), and that error over the tolerance;
  - its right-hand-side calls, Jacobians, factorizations and steps;
  - three flags: an `Ok` beyond criterion 2's bound, a non-finite state, and a conserved quantity
    drifting beyond its rounding bound.

  Each solver also gets a reliability profile: the share of runs with error/tol <= K, for K from
  1 to 1e3.
- **Output.** CSV files under `bench/results/`, and markdown tables generated from them. Both are
  regenerated each wave. From M13 on, the README carries the latest table.

---

## 8. Order of work

| Wave | Solver | Properties and robustness | Proofs | Benchmarks | Exit |
|---|---|---|---|---|---|
| 0 | M1: bench, references, CI matrix (Linux x86_64 and arm64, macOS arm64) | Harness; kernel properties; contract map; mutation tool and a baseline run | V0 | The work-precision baseline (roadmap Appendix B) | CI green: build, corpus, soak, properties, proofs, tripwire |
| 1 | M2: tolerances, controllers, the contract change with D16 and D17 | Controller and driver properties; budgets and partial results | V1 | Table re-pinned | Roadmap M2 acceptance; T3 and T4 checked |
| 2 | M3: in-place LU, modified Newton | LU differential and invariant properties | V2 | Gate | Roadmap M3 acceptance |
| 3 | M4, M5: variable order | Exactness, order, A(alpha) canaries | V3 | Gate | Roadmap M5 acceptance; a certificate per order, or a certified controller |
| 4 | M6, M6b, M7: Jacobians, banded and sparse LU, outputs and events | Adversarial suite, fault injection, domains; G1 | Catch-up | Gate | R1 to R5 on the corpus so far; G1's targets |
| 5 | M8: corpus, nightly gate | Nightly campaigns | — | Scorecard v1 | Roadmap 1.3 on the full ODE set; scorecard v1 published |
| 6 | M10, then M9: Radau IIA, DAEs | Radau and DAE properties; G2 | V4 | Scorecard v2, with IDA, and BDF against Radau | Roadmap M9 and M10 acceptance |
| 7 | M11 if needed; M13: release | Full mutation run; G3 go/no-go | V5 | The scorecard in the README | Section 1 holds |

What changes in the two accepted documents:

- M8's property tests and mutation script move to wave 0. The sed script becomes the
  `compiler-libs` tool (3.4).
- M6b, for sparse Jacobians and sparse LU, is new (D20).
- M10 comes before M9. Of D12's two conditions, the cross-check is the one a correctness goal
  needs.
- M12, the verification hooks, is spread over V0 to V3.
- The roadmap's bit-identity requirements become agreement within correctness bounds (D19). These
  are 4.4's bit-identical refactors and hex-dump detector, M3's bit-identical LU, M6's
  bit-for-bit run with a user Jacobian, and D3's bit-identical PRs.
- The FV pilot's T3 is restated for the M2 controller and moves to V1.

Effort is judgement, not calibrated against this code base. It is in person-days for one
developer who knows the code.

| Item | Effort |
|---|---|
| Roadmap M1 to M13, its revised range | 69-103 d |
| Property harness, generators, per-wave properties, contract map | 15-24 d |
| Mutation tool | 3-4 d |
| Adversarial suite and fault injection | 5-8 d |
| M6b | 8-12 d |
| G1 and G2. G3 adds 10-20 d, stopped after 5 on a no-go | 12-18 d |
| Scorecard v1 and v2 | 9-15 d |
| **Engineering total** | **121-184 d** |
| Proofs V0 to V5, for someone fluent in MathComp (the FV report's newcomer factor applies) | 12-30 weeks |

---

## 9. Decisions

These follow the roadmap's format; D1 to D14 stand. D15 to D18 and D20 are accepted as
recommended. D19 is decided.

| ID | Question | Recommendation |
|---|---|---|
| **D15** | What happens to an exception raised by the user's right-hand side? | It propagates unchanged, and the documentation says so. A non-finite value stays the recoverable signal: the step is rejected and retried, as today. Catching every exception would also catch `Out_of_memory`, `Stack_overflow` and `Sys.Break`. |
| **D16** | What does an `Error` carry? | The reason, the last accepted point and the statistics. This lands in M2, with D3's history on failure. The reasons gain `Budget`, and later event and DAE-initialization failures, instead of folding them into `StepRejected`. |
| **D17** | Integrate backward in time? | Yes, in M2. That is before T3 and T4 are stated, so both cover the two directions once. |
| **D18** | Which other solvers does the scorecard run? | scipy now. CVODE and IDA through `sundialsml`, once the SUNDIALS C library is installed outside the build. Nothing else, unless a gap appears. |
| **D19** | Bit-identical results across platforms, compilers or refactors? | No. The goal is correct results, not bug compatibility. Every check is a correctness bound, a typed outcome, or a value printed to the digits its bound makes meaningful, so the last bits may change between arm64 and x86_64, between compilers, and in a refactor. Determinism within a build stays (R5): it costs nothing in pure code and catches hidden state. `Float.fma` and compensated sums are used wherever they make a result more accurate. |
| **D20** | A sparse direct solver? | Yes, as M6b: after the banded solver and before Krylov. Circuits and networks (roadmap 1.1) are sparse without a band. |

---

## 10. Risks

- **T14 may not close for BDF4 and BDF5.** The fallback restricts the controller, and the
  work-precision gate shows what that costs.
- **Toolchain churn** (FV report section 5). The proof switch and the CI image are pinned. Expect
  one migration a year.
- **Certified mode costs two to three solves.** It is opt-in, and the scorecard reports it
  separately from the default mode.
- **The scorecard's fairness.** Scripts, versions and settings are published. Problems pass the
  transcription check, and every solver runs at its defaults and at its best settings.
- **References for new problems** follow roadmap 4.2: two solvers must agree. After M10, vstiff's
  two methods add a cross-check.
- **CI time.** The per-push run stays within its budget, with about a minute for properties.
  Everything heavier runs nightly.
- **No bit-level tripwire.** Without bit-identical refactors, a change of behaviour can hide below
  the printed digits. Properties, mutation testing and the work-precision gate carry that load.
- **Estimates** are judgement. The waves are ordered so that each one leaves the solver better and
  green, and stopping after any wave leaves a coherent state.
