# vstiff roadmap: from a small BDF2 integrator to a general stiff solver

Status: proposal. Code references and measurements are against commit f6d9b4e unless stated
otherwise (OCaml 5.5.0, dune 3.24.2, standard library only, no flambda, flat float arrays). Two
later commits fixed defects found by running the code: 052fdc7 (`Newton.solve` reports a
converged step that overflows as `Error Nan`) and dc3bb67 (`Adaptive.integrate` snaps steps to
the clock and rejects steps that cannot move `t`; the floor 16 eps |t| is now held by
[`Clock.resolution`](../../src/clock.ml)). The text states the current behaviour where they
matter and the baseline behaviour where a measurement predates them. The companion
[formal verification report](formal-verification.md), called the FV report here, covers Rocq and
MathComp.

How claims are marked:

- **[E]** measured with throwaway prototype programs built on a copy of the repository; the
  programs are not kept in the repository (Appendix A).
- **[D]** read in a primary document: the CVODE user guide (mathematical considerations) and
  sources (`cvode.c`, `cvode_impl.h`, `cvode_ls.c`), the DDASSL Fortran source at netlib, PETSc's
  `adaptdsp.c`, the Akrivis–Katsoprinakis paper on BDF stability angles, the Li–Liao arXiv
  abstract on variable-step BDF stability, the deTestSet package index.
- **[M]** from memory of the literature, not re-verified; check before relying on a number.
- **[R]** checked or re-measured a second time, independently of the first measurement or
  reading, with the same toolchain. **[?]** could not be verified.

---

## 0. Summary

vstiff is a readable, pure-functional variable-step BDF2 integrator with a halve/double step policy,
full Newton with a finite-difference Jacobian rebuilt at every iteration, and a dense recursive
Gaussian elimination. The expect-test discipline is its best asset. One defect was found while
measuring, and is fixed: at f6d9b4e [`Adaptive.integrate`](../../src/adaptive.ml) never returned
when every proposed step was below half an ulp of `t`, because the 16 eps |t| floor was applied on
rejection only, and a step of a few ulps advanced the state by its nominal length but the clock by
less (Sections 1.3 and 2.2) **[R]**; dc3bb67 fixed both. Measured against scipy's stiff solvers on
Robertson, HIRES, van der Pol (mu = 1000) and a Brusselator method-of-lines problem with n = 80, it
spends 36 to about 1250 times more right-hand-side calls than scipy's BDF at equal achieved error
(120x, 220x and 215x at rtol = 1e-6 on Robertson, HIRES and van der Pol; 650x on the Brusselator at
1e-4), counting scipy's finite-difference Jacobian evaluations as vstiff's counter does (scipy's
reported `nfev` leaves them out, see 2.0). At n = 80 about 80 % of a Newton iteration is the dense
solve (the finite-difference Jacobian takes most of the rest), 98 % once the Jacobian is reused
**[E, R]**.

Three prototypes (about 150 lines of OCaml; the controller prototype needs its own driver because
the current [`Ode.Controller`](../../src/ode.mli) cannot express it) show where the payoff is
**[E]**. The figures below were recomputed from the raw tables **[R]**; "equal achieved error" means
log-log interpolation of the work-precision curve over rtol 1e-3 ... 1e-8.

| change (BDF2 kept) | rhs calls | rejections |
|---|---|---|
| continuous error-based controller (I, safety 0.9, clamp [0.2, 5]) instead of halve/double | -27 % to -43 % at the same tolerance (-35 % to -43 % for rtol <= 1e-4), but it is 1.1-1.8x less accurate there; at equal achieved error -16 % to -32 % (-25 % to -31 % at rtol 1e-6) | from a quarter of all attempts (a third of the accepted steps) to 0-10 % of accepted steps (under 1 % at rtol <= 1e-5) |
| the same with DASSL's discrete rule (hold / double / shrink), the rule recommended for BDF in 2.2 | +14 % to -32 % at the same tolerance (1.5-3x more accurate); at equal achieved error -7 % to -29 % (+7 % on HIRES at rtol <= 1e-4) | 0-1.4 % of accepted steps |
| modified Newton, Jacobian carried in the history (same controller, same estimate) | n = 3: /2.1 to /2.8; n = 8: /2.0 to /6.7 (/5.9 at rtol 1e-6); n = 80: /42 (/65 if compared with Halving, which adds the controller's gain). CPU at n = 80 does not follow until the factorization is reused: 4.2 s against 3.1 s back to back, 11 against about 6 Newton iterations per step, each a full O(n^3) solve | unchanged |
| one Newton solve per step with a divided-difference error estimate | about /2 to /4.3 at equal achieved error for rtol <= 1e-6 (at looser tolerances the error ranges do not overlap); 6-10x fewer steps at rtol 1e-6, 14-21x at 1e-8, for an achieved error of 6x to 890x the tolerance (criterion 2 of 1.3 allows 100x) | 0-24 % of accepted steps (van der Pol at rtol 1e-3 and 1e-4: 24 %, 16 %) |

Taken together (the last row on top of modified Newton) the three changes cut the rhs calls at the
error Halving reaches at rtol 1e-6 by about 10x, 36x and 13x on Robertson, HIRES and van der Pol
(extrapolated just outside the measured error range). What remains is the order: scipy's
variable-order BDF (orders 1-5) still needs 4-9x fewer calls at the errors it reaches at rtol
1e-4 ... 1e-6 and 10-23x at 1e-8 (extrapolated). Variable order is the single largest item, and
it is where the representation decision (Section 2.1) matters most.

The roadmap keeps what is good: values instead of exceptions, effects quarantined in named
modules, contracts as module types composed with modular explicits, pins that are never
weakened. It changes the contracts where the measurements say the current ones cannot express
the needed behaviour (the controller never sees the error ratio; a failed step loses the
method's state; nothing counts work without a mutable cell).

---

## 1. Goals and non-goals

### 1.1 Who would use this

- OCaml projects that need a dependable stiff initial-value solver without C bindings: chemical
  kinetics (n = 3 to 100), circuit models (n = 10 to 1000), method-of-lines PDEs in one space
  dimension (n = 100 to 100 000, banded or sparse Jacobian), control and biology models.
- People who read numerical code to learn or to verify it: the code should remain the clearest
  BDF implementation available, and a target for the Rocq/MathComp investigation.
- Not a competitor to SUNDIALS or RADAU5 on speed for n > 1e5, and not a DAE solver for index
  above 1.

### 1.2 Accuracy range

Tolerances from 1e-3 to 1e-10 in a weighted RMS norm with rtol and atol vectors; achieved
accuracy at least that of scipy's BDF at the same tolerance on the corpus below. The criteria of
1.3 are stated down to 1e-8 only: add 1e-10 for the small problems (scipy's BDF already fails on
van der Pol at 1e-11, `success = False` **[R]**), and note that the divided-difference prototype
is less accurate than scipy's BDF at equal tolerance on van der Pol (2.1e-4 against 4.2e-5 at 1e-6).

### 1.3 Robust and general, in measurable terms

Robust means, for every problem of the corpus in Section 4.1 and every rtol in
{1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8} with atol as the problem prescribes:

1. The integration finishes with `Ok` (no give-up) and the state is finite. Appendix B uses
   rtol = atol, which for Robertson at rtol 1e-3 and 1e-4 is larger than y2 (at most 3.7e-5): scipy
   BDF and Radau report failure at rtol = atol = 1e-3 and LSODA returns NaN **[R]**, and the
   divided-difference prototype fails the same way (Section 2.3). The
   corpus table therefore has to prescribe atol per problem (Robertson: atol_2 <= 1e-8 is typical,
   <= 1e-4 sufficed at rtol 1e-3 in the prototype), and this criterion holds for those values; at
   rtol = atol on Robertson `Ok` or `Error` are both acceptable, a loop or a garbage `Ok` is not.
2. The achieved error at the end, `max_i |y_i - ref_i| / (rtol |ref_i| + atol_i)`, is at most
   100 (tolerance proportionality within two decades). This is scipy's level, not a margin: scipy
   BDF reaches 2.6x to 71x the tolerance at rtol <= 1e-6 on Robertson, HIRES and van der Pol (71x
   on van der Pol at 1e-8; LSODA 84x), the divided-difference prototype 6x to 890x **[R]**.
   Appendix B reports `max_i |y_i - ref_i| / (1 + |ref_i|)`, which is this measure divided by tol
   only when rtol = atol, and is blind to components far below 1 (Robertson's y2 is 4.8e-7 at
   t = 1e4).
3. Rejections (error test plus convergence failures) are under 10 % of accepted steps (the
   divided-difference prototype exceeds this on van der Pol at rtol 1e-3 and 1e-4: 24 %, 16 %
   **[R]**).
4. The solver never loops without advancing `t`, and the state advances by exactly what the clock
   does. f6d9b4e failed this: the 16 eps |t| floor was applied on rejection only, and
   [`Adaptive.integrate`](../../src/adaptive.ml) did not return when every proposed step was below
   half an ulp of `t` (t0 = 1, t_end = 1 + ulp; t0 = 1e10 with dt_max = 5e-7; `t` never moved in 3e6
   rhs calls) **[R]**; and a step of a few ulps advanced the state by its nominal length but the
   clock by less (t0 = 1e15, a span of 100, dt_max = 0.19: `Ok` with y = 76.36 instead of 100).
   Since dc3bb67 it holds: every non-final step is snapped to the floats, h = (t + dt) - t; a step
   that cannot move `t` (h <= 0) is rejected as `Too_small` without calling the method; and a
   remainder at or below [`Clock.resolution t`](../../src/clock.ml) = 16 eps |t| is taken as the
   last step. The controller decides what a rejection does: [`Halving`](../../src/halving.ml) halves
   the step, which is below the floor, so it gives up with `StepRejected`. A step that does move `t`
   is judged by the error test alone, however short, so the floor binds only on rejections.
   Termination follows from two facts: every accepted non-final step strictly increases `t`, and the
   controller ends any run of rejections (FV report, T4); a new controller has to keep the second.
   Corpus lines 25 to 27 pin the three cases;
   [numerics/05-step-control.md](../numerics/05-step-control.md) explains the rules.
5. Work-precision: at equal achieved error, right-hand-side calls are within 2x of scipy's BDF
   with a finite-difference Jacobian (its Jacobian evaluations counted, see 2.0), and within
   1.5x once a user Jacobian is supplied.

General means: rtol and atol vectors; variable order 1-5; user-supplied or finite-difference
Jacobian, dense or banded; mass matrix and index-1 DAEs; dense output, output at requested
times, event location; two methods (BDF and Radau IIA) cross-checking each other; every
statistic a paper would report.

### 1.4 Non-goals

Partitioned or IMEX methods, Krylov linear solvers before a user needs them (kept as a late
option), parallelism, GPU, automatic differentiation, index-2 and index-3 DAEs, extended
precision, Fortran-level throughput.

---

## 2. Capabilities and their design

### 2.0 Where the current design stands (measured baseline at f6d9b4e) [E]

Corpus test: 1.4 s CPU; soak test: 12.4 s CPU on an idle machine (direct runs of the test
executables; 1.8 s and 16.6 s were measured earlier) **[R]**. The corpus then had 23 lines (the
9 original cases and 14 regression pins); it has 27 now, 18 of them regression pins, and the
soak test has 4. CPU times on the measuring machine vary by up to 2.5x with load (Apple
silicon, other jobs running): ratios between runs made back to back are reliable, absolute
times are not, so acceptance criteria below are stated in counts.

Linear solver: the pure recursive elimination costs about 4x an in-place LU with partial pivoting at
every size tried (3.5-4.5x over repeated runs; idle: n = 10: 1.9 vs 0.5 us; n = 100: 0.97 vs 0.22
ms; n = 300: 24 vs 6.5 ms). The prototype LU agrees with the recursion to 1e-17 but not bit for bit,
because its back substitution sums in a different order. An in-place LU that performs the same
floating-point operations in the same order (elimination applied to the right-hand-side column stage
by stage, back substitution as `(b_i - s) / pivot` with `s` the products rounded one by one and
added left to right from 0) is bit-identical to [`Linalg.solve`](../../src/linalg.ml): 20 119
random, badly scaled, singular and Laplacian systems of sizes 1-120 gave no mismatch on an arm64
build **[R]**. That is a property of this compiler's instruction selection (no fused multiply-add in
either loop), not a language guarantee.

Work-precision at rtol = atol = 1e-6 (full tables in Appendix B):

| problem | vstiff BDF2 + Halving: steps (rejected), rhs calls, achieved error | scipy BDF (FD Jacobian): steps, rhs calls counted with a wrapper (reported `nfev`), error |
|---|---|---|
| Robertson, t = 1e4 | 2531 (831), 84 267, 6.5e-7 | 176, 589 (548), 2.6e-6 |
| HIRES, n = 8 | 3212 (1059), 187 292, 1.0e-6 | 164, 649 (450), 1.0e-5 |
| van der Pol, mu = 1000, t = 2000 | 39 041 (12 984), 978 933, 2.8e-6 | 847, 2791 (2620), 4.2e-5 |
| Brusselator MOL, n = 80, rtol 1e-4 | 737 (230), 453 298, 6.0e-5, 3.0 s idle (7.9 s measured earlier, probably under load) | 94, 688 (283), 1.7e-4, 0.03 s |

In the scipy column the first numbers are counted with a wrapper around the right-hand side,
because scipy's reported `nfev` for BDF and Radau leaves out the finite-difference Jacobian
columns (njev x n, plus a few extra evaluations), while vstiff's counter includes them; LSODA's
`nfev` includes them **[R]**. Comparing vstiff's count with the reported `nfev` would overstate
the gap to scipy by 7-60 % on n <= 8 and by 2.4x at n = 80.

Two readings. First, the halving controller rejects a quarter to a third of all steps because
every rejection halves and every three accepts double, so the step oscillates around the
largest acceptable one. Second, the finite-difference Jacobian at every Newton iteration costs
(n + 1) calls each, which is why the n = 80 problem needs 450 000 calls for 737 steps.

### 2.1 Variable-order BDF, orders 1 to 5

**Representation.** Four formulations; the last three are in production use
**[M, with [D] where marked]**:

| | variable-coefficient (what Bdf2 does today) | fixed-leading-coefficient, Nordsieck array (CVODE **[D]**, VODE) | fixed-leading-coefficient, modified divided differences (DASSL **[D]**) | quasi-constant step, backward differences (ode15s, scipy BDF) |
|---|---|---|---|---|
| history carried | the last k+1 points (t_i, y_i) | k+1 scaled derivative vectors z_j = h^j y^(j)/j! plus the last k step sizes | phi_j vectors plus psi_j step sums | k+1 backward differences, re-interpolated when h changes |
| coefficients | from the k step ratios, each step | from the ratios (l, tq vectors) | from psi (alpha, beta, gamma, sigma) | fixed per order |
| gamma in the stage equation x = psi + gamma f | beta(ratios) h: changes whenever a ratio changes | h / l_1: depends on h only | h / alpha_s: depends on h only | h / gamma_k: depends on h only |
| step change | free (zero-stable while ratios stay under a bound) | rescale z | nothing | interpolate the differences |
| readability | highest: "the polynomial through the last k+1 points" | medium | lowest | medium |
| error and order estimates | divided differences of the points | from z and the last correction | from phi | from the differences |

Recommendation: **stay variable-coefficient and generalise what Bdf2 already does**, with the
history as the last k+2 old points (k+1 at the maximum order). Reasons:

1. It is the existing formulation ([`Bdf2.coeffs omega`](../../src/bdf2.ml)), so order 3-5 is a
   generalisation, not a rewrite, and no coefficient tables appear anywhere: the stage equation
   is derived by evaluating the interpolating polynomial, Section 2.1.1.
2. Its zero-stability under step changes is the best of the families compared (Hairer–
   Wanner II, V.5 **[M]**), though not firmly, because the fixed-leading-coefficient figures
   below come from a reconstruction. Largest constant step ratio for which the homogeneous
   recurrence satisfies the root condition, variable-coefficient formulas **[R]**: BDF2 2.414
   (= 1 + sqrt 2, Grigorieff 1983, sharp **[D]**: confirmed in the literature search, not in the
   Li–Liao abstract), BDF3 1.618, BDF4 1.28, BDF5 1.127. A fixed-leading-coefficient corrector
   built on the degree-k predictor polynomial (a reconstruction of the Jackson–Sacks-Davis form,
   not checked against their paper **[?]**) gives 1.92, 1.41, 1.21, 1.12 **[R]**. At a sustained ratio 2 the parasitic roots have modulus 1.43, 2.83 and 5.63 per step for BDF3,
   BDF4 and BDF5 (variable-coefficient), so a ratio above the limit is acceptable only as a single
   transient followed by constant steps. That is why DASSL changes the order only after k+1 steps
   at constant step size and order (`NS`, `KP1.GE.NS` in `DDASTP`) and at most doubles the step
   **[D]**. Li–Liao (arXiv 2201.00527: "the variable-step BDF3 scheme is stable if almost all
   adjacent step ratios are less than 2.553" **[D]**) concerns stability for dissipative
   parabolic problems on arbitrary grids, a different notion; it does not license a sustained
   ratio above the zero-stability limits.
3. The cost of the fixed-leading-coefficient trick (gamma depends on h alone, so the Newton
   matrix survives ratio changes) is small once the step policy is discrete (Section 2.2): with
   DASSL's rule the ratio is exactly 1 on most steps, so gamma is constant on most steps in
   either family. A step change itself moves gamma by about the step ratio in both families
   (variable-coefficient: x1.6-1.8 at a doubling, x0.56-0.64 at a halving, orders 2-5), so the
   matrix is refactored there either way. In the k steps after it the variable-coefficient gamma
   keeps drifting: by at most 21 % after a doubling and 28 % after a halving (order 5), 58 % after
   a ratio of 0.25 **[R]**, against CVODE's 30 % threshold and DASSL's [0.6, 5/3].
4. It is the formulation closest to a formal statement: "exact on polynomials of degree <= k
   through these nodes" is one lemma about Lagrange interpolation (Section 4.6).

The Nordsieck fixed-leading-coefficient form is the fallback if, after M5, factorization reuse
turns out to be limited by gamma drift; the decision is D1 in Section 6.

#### 2.1.1 The stage equation for order k, derived not tabulated

With nodes t_{n+1} > t_n > ... > t_{n+1-k} and the unknown y_{n+1}, the BDF step is: the
polynomial P of degree k through (t_{n+1}, y_{n+1}), (t_n, y_n), ..., (t_{n+1-k}, y_{n+1-k})
satisfies P'(t_{n+1}) = f(t_{n+1}, y_{n+1}). Write P in Newton form with the new node first;
every divided difference is affine in y_{n+1}, so P'(t_{n+1}) = a y_{n+1} + b with a > 0 a
scalar and b a vector computed from the old points. Then

    y_{n+1} = psi + gamma f(t_{n+1}, y_{n+1}),  gamma = 1 / a,  psi = -b / a,

which is exactly [`Stage.equation`](../../src/stage.ml). For k = 2 and ratio omega this reproduces
[`Bdf2.coeffs`](../../src/bdf2.ml) (a1, a0, beta), which is the first unit test of the
generalisation. The computation of (a, b) is O(k^2) scalar operations plus k vector operations per
step.

#### 2.1.2 Error estimate and order selection from divided differences

For the exact solution, the derivative of the interpolation error at the node t_{n+1} is
omega'(t_{n+1}) y^(k+1)(xi) / (k+1)! with omega'(t_{n+1}) = prod_{i=1..k} (t_{n+1} - t_{n+1-i}),
so the local error of the step is

    LTE_k ~ (omega'_k(t_{n+1}) / a_k) * D_{k+1},
    D_{k+1} = [y_{n+1}, y_n, ..., y_{n-k}]   (the (k+1)-th divided difference, including the new point).

Check for k = 2 with constant h (omega' = 2 h^2, D_3 = y'''/6 + O(h)): the formula gives
LTE = (2/9) h^3 y''', the known BDF2 constant; for k = 3 the same formula gives (3/22) h^4 y''''
**[R: derived and checked]**. This is DASSL's `ERK = SIGMA(K+1) * ENORM`
in another notation **[D]**; the prototype used the equivalent predictor-corrector
form `cc / (cp - cc) * (y - y_pred)` for k = 2.

The same differences give the estimates for the neighbouring orders at no extra solves: D_k
with omega'_{k-1} and a_{k-1} for order k-1, and D_{k+2} (needs one more old point) for order
k+1. Order selection follows DASSL, verified in `DDASTP` **[D, R]**, with TERK_j = (j+1) * est_j:
at every attempt, lower the order when (k = 2) TERK_1 <= 0.5 TERK_2, or (k >= 3)
max(TERK_{k-1}, TERK_{k-2}) <= TERK_k; never below 1. After an accepted step the order is
lowered if that rule said so; kept if k is the maximum, the order was raised on the previous
step, or the step size and order have not been constant for k+1 completed steps; otherwise at
k = 1 raised when TERK_2 < 0.5 TERK_1, and at k >= 2 lowered when TERK_{k-1} <= min(TERK_k,
TERK_{k+1}), else raised when TERK_{k+1} < TERK_k, else kept. The next step size uses the estimate
of the order just chosen. CVODE's rule is the same shape: consider a change only after q+1 steps
at order q, and pick the order whose step ratio is largest **[D]**.

Startup: order 1 (backward Euler) from y0 with half the gap to the explicit Euler predictor as
the estimate (as today), then raise the order and double the step on every accepted step until
the first failure, an order-lowering decision or the order limit (DASSL's phase 0, `IPHASE = 0`
**[D]**). Check the doubling before copying it under the variable-coefficient representation: a
ratio of 2 is above the zero-stability limits of BDF3-5 (1.62, 1.28, 1.13 above; parasitic roots
of modulus 1.43, 2.83 and 5.63 per step at a sustained ratio), so the 3-4 steps of the ramp can
amplify noise by up to about 20x. DASSL ramps the same way with a corrector that is less tolerant
by that reconstruction, so the effect may be harmless; a cap below the limit of the order in use
costs a few steps **[R]**. The current "BDF2 after one step" is the k = 2 special case.

Stability limits: BDF1 and BDF2 are A-stable; BDF3-6 are A(alpha)-stable with maximum angles
86.03, 73.35, 51.84 and 17.84 degrees (Akrivis–Katsoprinakis, BIT 60 (2020) 93-99 **[D]**; the
same values appear in Hairer–Wanner II **[M]**). Default maximum order 5, as CVODE **[D]**; BDF6
is not worth its 17.84 degrees. The corpus gets a canary for this: y' = lambda y with
lambda = |lambda| e^{i (pi - theta)}, an eigenvalue at angle theta from the negative real axis,
|h lambda| large. theta = 60 degrees is outside BDF5's sector (51.84) and inside BDF4's (73.35),
so the order selector must settle at order <= 4; theta = 80 degrees (outside BDF4's, inside
BDF3's 86.03) at order <= 3. Check the orders used at the end of the run, once |h lambda| is
large, not in the first steps.

#### 2.1.3 History as a list of points

```ocaml
(* bdf.mli (new, replaces bdf1.ml/bdf2.ml as the production method; both stay as Methods
   for the fixed-step tests) *)
type history = {
  points : Ode.point list;     (** newest first, at most max_order + 1 of them: order k < max_order needs k + 2 *)
  order : int;                 (** the order of the next step *)
  held : int;                  (** steps taken at this order and step size, for the k+1 rule *)
  cache : Stage.cache;         (** Jacobian, factorization and their age; see 2.3 *)
}
```

Nothing is mutated; a step conses a point and drops the oldest. Variable order fits the
existing abstraction because `history` was always "what the method carries", and the driver
never looks inside.

### 2.2 Error-based step and order controller

**Tolerances and norm.** rtol scalar, atol vector (scalar atol is the vector with equal
entries), weights W_i = 1 / (rtol |y_i| + atol_i) and the WRMS norm
`sqrt (1/n sum (W_i v_i)^2)`, as CVODE **[D]**. The weights are computed from the state at the
start of the step (CVODE) so they are the same for every candidate order. An error estimate
passes when its norm is <= 1. The current `max_i |err_i| / (1 + |y_i|) <= tol` is the special
case rtol = atol = tol with the max norm and the weights taken at the end of the step. RMS <= max
<= sqrt(n) RMS, so the RMS test is up to sqrt n times looser at the same tol (n = 80: up to 9x);
on n <= 8 the RMS prototype took 5-20 % fewer steps for 1.3-1.9x the error, and it failed on
Robertson at rtol = atol = 1e-3 and 1e-4 **[E, R]**. scipy's BDF uses the RMS norm.

**Step-size formula.** With r = ||est||_WRMS for the estimate of order k (error ~ h^(k+1)):

    h_new = h * clamp(fac_min, fac_max, safety * r^(-1/(k+1)))

DASSL's discrete version, verified in the source **[D]** (`est` and `k` are those of the order
chosen for the next step): r' = (2 est + 1e-4)^(-1/(k+1)); if r' >= 2 double h; if r' <= 1 take
r' := max(0.5, min(0.9, r')) and shrink; else keep h. CVODE keeps h and q unless the best
candidate ratio exceeds 1.5 (`ETA_MAX_FX`), and caps growth at 10 (1e4 on the first step)
**[D]**. The prototype used the continuous formula with safety 0.9 and clamp [0.2, 5] (fac_max
= 1 on the step after a rejection); its rhs-call and rejection figures are in the table of
Section 0 **[E, R]**. With the discrete rule on the same (backward Euler gap) estimate the dead
zone r' in (1, 2) makes the step sit 2-8x below the tolerance: rhs calls at rtol 1e-6 are
89.6k / 181k / 828k on Robertson / HIRES / van der Pol, against 49.1k / 110.9k / 592k for the
continuous rule and 84.3k / 187k / 979k for Halving, with 1.5-3x smaller errors **[R; table in
Appendix B]**.

Recommendation for BDF: **DASSL's discrete rule** (hold, double, or shrink), accepting that it
buys stability and factorization reuse, not rhs calls. A ratio of exactly 1 on most steps keeps
the variable-coefficient method inside its stability bounds (2.1), keeps gamma constant so
factorizations are reused (2.3), and keeps the history polynomial well conditioned; the
prototypes do not measure the factorization reuse. The continuous rule is the cheaper one at a
given tolerance, and the M2 acceptance numbers belong to it.

**After failures** **[D, DASSL source]**: first error-test failure: the order drops to the one the
order test above chose (`K = KNEW`) and h *= clamp(0.25, 0.9, 0.9 r') with r' from that order's
estimate; second: same order rule, h *= 0.25; third and later: order 1 and h *= 0.25 (CVODE: eta in
[0.1, 0.2] from the second failure; after the third, the order drops by one per failure, restarting
at order 1, with eta = 0.1; give up after seven **[D]**). Newton convergence failure with a current
matrix, or a singular matrix: h *= 0.25 (DASSL, and CVODE's eta_cf = 0.25 **[D]**); with an
out-of-date matrix both first redo the step at the same h with a fresh matrix. Give up after 10
convergence failures in a row (both) or when h < 16 eps |t|. The floor rule alone was not enough: it
is applied on rejection only, so at f6d9b4e [`Adaptive.integrate`](../../src/adaptive.ml) never
returned when all proposed steps were below half an ulp of `t` (the accept-path cases of 1.3,
criterion 4 **[R]**); DASSL tests its `HMIN` on retries only as well **[D]**. Since dc3bb67 the
driver covers the accept path itself: a step that cannot move `t` (h <= 0 after snapping) is
rejected as `Too_small` without calling the method, and the controller decides.
[`Halving`](../../src/halving.ml) halves it, which is below its floor
[`Clock.resolution t`](../../src/clock.ml) = 16 eps |t| (the same rule as before, now defined once
in `Clock`), so it gives up with `StepRejected`. A new controller has to end every run of rejections
in the same way (FV report, T4).

**PI and digital filters.** Soderlind's filters, with coefficients verified in PETSc's
`adaptdsp.c` **[D]**: step ratio rho = (1/r_n)^(b1/kappa) (1/r_{n-1})^(b2/kappa)
(1/r_{n-2})^(b3/kappa) (rho_{n-1})^(-a2) (rho_{n-2})^(-a3), with (b1, b2, b3 | a2, a3) scaled by
the listed factor and kappa the exponent of h in the local error (PETSc's `order`: kappa = k+1 for
the BDF of order k, 3 for BDF2, not k): PI34 = (7, -4, 0 | 0, 0)/10, PI42 = (3, -1, 0 | 0, 0)/5,
H211b = (1, 1, 0 | 1, 0)/4, H211PI = (1, 1, 0 | 0, 0)/6, H312b = (1, 2, 1 | 3, 1)/8. The PI34
prototype gave smoother step sequences and fewer rejections on van der Pol (63 vs 112 at rtol
1e-4) at 14-46 % more steps (typically 30 %) **[E, R]**. Worth trying for Radau. For BDF2 the
prototype shows no difference in rhs calls between PI34 and the discrete rule (within 10 % at
rtol 1e-6); the case for the discrete rule is stability and factorization reuse, which the
prototype does not measure.

**Contract change (required).** Today `acceptable : t -> y:Vec.t -> err:Vec.t -> bool` (here `t` is
the controller state) returns a bool and `accepted : t -> t` sees no error ratio, so an
error-proportional controller cannot be written against the current
[`Ode.Controller`](../../src/ode.mli). The prototype needed its own driver. Section 3 gives the new
signature: one `decide` that sees the attempt and returns the next state with the verdict; what it
does to the existing pins is in 3.3.

### 2.3 Newton solver

**Modified Newton.** Evaluate J once per Jacobian refresh, form M = I - gamma J (or
M_mass - gamma J), factor once, iterate x <- x - M^-1 G(x) with the factorization reused.
Predictor: the history polynomial extrapolated to t_{n+1} (today: the line through the last
two points; the quadratic extrapolation in the prototype cut iterations and is what the
divided-difference history gives for free).

**Convergence test tied to the tolerance, with a rate estimate.** DASSL **[D, R]**: rate
rho = (||delta_m|| / ||delta_0||)^(1/m), converged when rho/(1 - rho) ||delta_m||_WRMS <= 0.33
(after the first iteration the previous step's factor rho/(1 - rho) is used, 100 after a matrix
update), at most 4 iterations, failure when rho > 0.9; CVODE **[D]**: at most 3 iterations,
divergence declared when ||delta_m|| / ||delta_{m-1}|| > 2 (m > 1), convergence when
R ||delta_m||_WRMS < 0.1 eps, with R a convergence-rate estimate and eps the local-error
acceptance threshold of the same norm. The norm is the WRMS norm with the step's weights: the
corrector only needs to be converged well below the error tolerance, not to 1e-10 in the step.
Its benefit is iterations per step (the prototypes need 2-15 with the fixed 1e-10 test; CVODE and
DASSL cap at 3-4); that it is needed for robustness is **not shown** (next paragraph).

**The prototype's Robertson failure is not evidence for it.** The divided-difference prototype
with the fixed `1e-10 (1 + |x|)` step test failed on Robertson at rtol = atol = 1e-3 (its RMS
variant also at 1e-4): h fell from 1e-9 to 3e-17 with r held in 0.75-1.03 at t = 8.7265e-3 until
the 16 eps |t| rule gave up (a trace of 382 attempts). Printing the state shows what happens
**[R]**: y2 goes from +3.8e-5 (true value 3.6e-5) at t = 3.0e-3 to -4.9e-5 at
t = 7.9e-3 in one accepted step of 4.9e-3 (r = 2.7e-3: the absolute tolerance 1e-3 is 26 times
y2), then runs off (-1.1e-3 at t = 8.70e-3, -1.5e-2 at 8.72e-3) because the -3e7 y2^2 term is
unstable on the negative side; |err| reaches 1e5 while r stays near 1 because the tolerance
tol (1 + |y|) grows with |y|. The cause is an absolute tolerance larger than y2, so the error
test cannot see the sign error; it is not the corrector. With the Newton test unchanged the same
prototype finishes at rtol 1e-3 for every atol from 1e-4 down to 1e-10 (70-127 steps); the
modified-Newton prototypes with the library's estimate finish at every atol; and scipy's BDF and
Radau report failure and LSODA returns NaN at rtol = atol = 1e-3 **[R]**. Of the
four divided-difference variants run at rtol = atol = 1e-3, I(1/3) and its RMS variant fail while
PI34 and the discrete rule finish, which fits a chance excursion below zero, not a systematic
defect. A solver that gives up with `Error` there behaves correctly (DASSL's `NONNEG` option
exists for this: an iterate with a negative part above 0.33 in the WRMS norm counts as a
convergence failure **[D]**). The pins should be Robertson at rtol 1e-3 with the corpus atol
(`Ok`) and at rtol = atol = 1e-3 (`Ok` or `Error`, never a loop, never a garbage `Ok`). A less
conservative estimate (M4) is what makes such excursions possible: all five prototype variants
that use the backward Euler gap finish at rtol = atol = 1e-3.

The same trace shows an unrelated safeguard that is still worth having: after a convergence
failure the prototype re-evaluated J at the same predictor point and grew the step by 5x after
the next accept, so the integration oscillated between "accept, accept, fail, fail". Hence:

**Refresh policy** (CVODE **[D]** unless marked): rebuild M when more than 20 steps have
passed since the last build or when |gamma / gamma_bar - 1| > 0.3, and after a convergence
failure. When the bound is exceeded and J is younger than 50 steps, reuse the saved J and re-form
and refactor M = I - gamma J: a step change costs a factorization but no rhs calls. While gamma
stays within the bound the old factorization is kept and the correction is scaled by
2/(1 + gamma/gamma_bar) for BDF (CVODE `cvLsSolve`; DASSL's analogue is a residual scaling by
2/(1 + cj/cjold), `DDASTP`) **[D, R]**. Re-evaluate J itself after 50 steps, or after a
convergence failure with an out-of-date matrix and |gamma/gamma_bar - 1| < 0.2 (CVODE; DASSL
redoes the step with a fresh matrix after such a failure) **[D]**. If J was already fresh, the
cure is h *= 0.25, not another J. DASSL's equivalent of the gamma test is cj / cjold outside
[0.6, 5/3] **[D]**. No growth of h on the step after a failure (prototype lesson, not a DASSL or
CVODE rule), and DASSL's rule of changing the order only after k+1 steps at constant step size
and order.

**Line search.** Not inside the integrator: the step-size cut plays that role, as in every
production code **[M]**. The present damped Newton with Armijo backtracking stays as
[`Newton.solve`](../../src/newton.ml) for consistent initial conditions (DAEs, Section 2.6) and for
the unit tests.

**User Jacobian.** `problem.jac : (float -> Vec.t -> fy:Vec.t -> Mat.t) option`; `None` means
finite differences. Finite differences with column grouping (Curtis–Powell–Reid **[M]**):
for a banded Jacobian with lower and upper bandwidths l, u, perturb columns j, j + (l+u+1),
j + 2(l+u+1), ... together, so a Jacobian costs l + u + 2 calls (l + u + 1 beyond f(t, y)
itself) instead of n + 1. For general sparsity the grouping is a graph colouring of the column
intersection graph (greedy is fine).

Measured effect of Jacobian reuse alone, BDF2 kept, same controller, same estimate **[E]**:

| problem, rtol 1e-6 | rhs calls, J every iteration | rhs calls, J carried in history | Jacobians evaluated |
|---|---|---|---|
| Robertson (n = 3) | 49 055 | 20 069 | 16 |
| HIRES (n = 8) | 110 851 | 18 939 | 19 |
| van der Pol (n = 2) | 592 084 | 231 437 | 42 |
| Brusselator (n = 80, rtol 1e-4) | 292 495 | 6 932 | 3 |

At n = 80 the factor is 42; comparing the carried-J run with Halving
(453 298 calls) gives 65, which includes the controller's gain. CPU does not follow the rhs counts there:
back to back, the carried-J run takes 4.2 s and the full-Newton run 3.1 s, because the
stale-Jacobian iteration needs 11.1 instead of about 6 Newton iterations per step and every
iteration is a full O(n^3) solve (the saved Jacobian is 18 % of a library iteration). The CPU
gain needs the factorization reuse of M3 **[R]**.

### 2.4 Linear algebra

Measured: pure recursion about 4x slower than in-place LU at every n; at n = 80 the solves are
98 % of CPU time in the Jacobian-reusing prototypes and about 80 % of a library Newton iteration
(solve 0.8 ms, forward-difference Jacobian 0.13-0.16 ms, forming I - gamma J 0.03 ms, idle
machine); idle, the divided-difference prototype takes 0.4-1.3 s per run (1.2-3.6 s measured
earlier) and the carried-J and library runs about 3 s, against 0.03-0.09 s for scipy, whose LU
is LAPACK's **[E, R]**. Three steps, all standard-library:

1. **Dense LU behind a pure interface.** `Linalg.factor : matrix -> factored option` copies
   the matrix once and eliminates in place on the copy; `Linalg.solve : factored -> Vec.t ->
   Vec.t` allocates its result. No argument is ever written; the mutation is on arrays the
   function allocated itself and returns frozen (nothing writes them after the return; the
   factorization escapes into the cache, so "never escape" would be false), so nothing
   observable changes (whether this counts as quarantined is D4). With the operation order of 2.0 the result is
   bit-identical to the recursion, so this step lands as a pure refactor **[R]**. Expected: 4x on
   the solve, and with factorization reuse (2.3) the per-iteration cost drops from O(n^3) to
   O(n^2): at n = 80 the 3144 solves of the divided-difference prototype at rtol 1e-6 (0.5-1.15 ms
   each, depending on load) become one factorization per gamma change plus 3144 triangular
   solves at about 10 us **[M for the 10 us]**.
2. **Banded LU** (LAPACK dgbtrf layout, `Band.t = { lower; upper; data : float array }`),
   O(n (l+u)^2) per factorization. The Brusselator with interleaved u, v has l = u = 2 (central
   three-point stencil, computed from the Jacobian pattern **[R]**), so grouped differences cost 5
   rhs calls per Jacobian beyond f(t, y); with banded storage and grouped differences, n = 2000
   should become a sub-second problem (a target, not a measurement).
3. **Large n**: matrix-free GMRES with a user preconditioner (CVODE's SPGMR route **[M]**) as
   the stdlib-only path. A sparse direct LU (KLU-like) is a project in itself; binding
   SuiteSparse or using `lacaml`/`owl` (both on opam, both drag BLAS/LAPACK) is the
   alternative. Recommendation: no dependency until a user has a sparse problem; `sundialsml`
   (SUNDIALS binding, on opam as 6.1.1p1, which needs a matching system SUNDIALS 6.x through
   `conf-sundials`; whether it builds in the development switch was not checked **[?]**) would be a test-only oracle,
   never a library dependency; scipy remains the practical oracle.

Bigarray: not needed. `float array` is already unboxed and flat (`flat_float_array: true`
in the development switch **[E]**); Bigarray only pays when handing memory to C.

### 2.5 Outputs

- **Dense output** is free with the point history: the interpolating polynomial through the
  last k+1 points, evaluated by Neville or in Newton form, is the method's own continuous
  extension (the same for DASSL: `DDATRP` evaluates the step's own polynomial **[D]**).
  `Bdf.interpolate : history -> float -> Vec.t`.
- **Output at requested times** by interpolation, never by cutting the step (cutting breaks the step
  sequence and the k+1 hold rule; today [`Adaptive`](../../src/adaptive.ml) cuts only the last step,
  which is acceptable, but a list of output times must not). API: `Adaptive.integrate ~at:[t1; t2;
  ...]` returning the points, plus `Adaptive.steps : ... -> step Seq.t` for callers who want every
  step with its interpolant (pure, lazy, no callbacks).
- **Events**: g(t, y) sign changes located on the interpolant by Illinois/regula falsi to a
  tolerance in t (CVODE's rootfinding uses the Illinois algorithm with tolerance
  100 u (|t_n| + |h|) **[D]**; the root of the interpolant can be located to that tolerance, the
  true event time is only as accurate as the solution, about tol relative); the integration
  stops at or restarts from the event with a fresh history (order 1). `events : (float ->
  Vec.t -> float) list` in the problem, `Event of int * point` in the result.

### 2.6 Index-1 DAEs and mass matrices

Whether: yes, the test set's circuit and chemistry problems are index-1 DAEs (transistor
amplifier, two-bit adder, NAND gate, Chemical Akzo Nobel **[D, deTestSet index and Test Set
pages]**) and "general" is hollow without M y' = f. When: after variable order (M5) and banded
storage (M6), because the machinery is the same: the stage equation becomes

    M (y_{n+1} - psi) = gamma f(t_{n+1}, y_{n+1}),   iteration matrix M - gamma J,

and the error test excludes nothing for index 1 (Brenan–Campbell–Petzold **[M]**; for index 2 the
algebraic components would have to be excluded, out of scope). Consistent initial conditions are the
user's job; a helper that solves f_alg(y0) = 0 with [`Newton.solve`](../../src/newton.ml) is cheap.
Radau IIA extends to index 1 with the same change **[M]**.

### 2.7 Radau IIA as a second method

3-stage, order 5, stiffly accurate, L-stable (Hairer–Wanner II, IV.8 **[M]**). The 3n
nonlinear system is solved by simplified Newton after the similarity transform that
diagonalizes A^-1 into one real eigenvalue and one complex-conjugate pair, so a step costs one
real and one complex n x n factorization per Jacobian (`Complex` is in the standard library; a
complex LU still needs its own code, preferably on separate real and imaginary arrays). Error
estimate: the embedded formula with the (I - gamma_0 h J)^-1 filter, with the second filtered
estimate on rejection (RADAU5 **[M]**). Controller: PI (2.2). Dense output: the collocation
polynomial.

Why a second method at all: on van der Pol scipy's Radau needs 5411 calls for 1.4e-7 and 15 003
for 6.9e-10, where BDF needs 5706 for 7.1e-7 (about 7400 for 1.4e-7 and 17 300 for 7e-10 by
log-log extrapolation; true counts) **[E, R]**: Radau is 1.15-1.4x cheaper at equal error there,
not the order of magnitude the unequal comparison suggests. The stronger reason is that two independent
methods agreeing is the best regression oracle the corpus can have (Section 4.2).

### 2.8 Statistics and instrumentation, without mutable state

Replace [`Instrument.count`](../../src/instrument.ml) by cost accounting as values: each step
returns a `Cost.t = { rhs_calls; jac_evals; factorizations; solves; newton_iterations }` (a module
of its own: `Ode` is interface-only and cannot hold `zero` and `add`) and the driver sums them.
Statistics then are pure data: steps, accepted, rejected by cause, cost, order histogram, min and
max h, and they are deterministic and pinnable. The price is that every rhs call must be counted by
hand through [`Stage`](../../src/stage.ml), [`Jac`](../../src/jac.ml) and the methods, where a
wrapped rhs is exact by construction: keep `Instrument` (or Guard's own counter) in the tests as the
oracle, and make "driver-summed `rhs_calls` equals the wrapped count on every corpus problem" a
corpus line. `Instrument` can go from the library once no caller needs it (Guard keeps its own
budget wrapper in the tests; that is the one place a counter is natural).

---

## 3. How the contracts evolve

Signatures are proposals for the next versions; names follow the repository (CamelCase module types,
no functors, no packed first-class modules beyond the modular explicits already used). Labels are
used more liberally below than the repository's rule allows (labels only where same-typed arguments
could be swapped: `weights`, `order`, `gamma`, `mass`, `h` in `decide` and `failed`, and `n` carry
none) and should be pruned when the signatures are implemented. A compile check against OCaml 5.5
showed two constraints that the signatures respect: `val`s cannot live in
[`ode.mli`](../../src/ode.mli), which has no implementation, and a module that holds only a module
type (`linsolve.mli`) must be listed in `modules_without_implementation` too **[R]**.

### 3.1 Problems, tolerances, costs

```ocaml
(* ode.mli *)
type rhs = float -> Vec.t -> Vec.t

(** Jacobian of [rhs] at [(t, y)]; [fy] is [rhs t y], already evaluated. *)
type jacobian = float -> Vec.t -> fy:Vec.t -> Mat.t

type problem = {
  rhs : rhs;
  jac : jacobian option;          (** [None]: finite differences, grouped by [structure] *)
  structure : Mat.structure;      (** [Dense | Banded of { lower : int; upper : int }] *)
  mass : Mat.t option;            (** [M y' = rhs]; [None] is the identity (M9) *)
  t0 : float;
  t_end : float;
  y0 : Vec.t;
}
(* Three new fields break every record literal in test/ and in client code; a
   `Problem.make ?jac ?structure ?mass ~rhs ~t0 ~t_end y0` in a module of its own (`Ode` has no
   implementation) keeps them short. *)

type point = { t : float; y : Vec.t }

(* Work done, summed by the drivers, is `Cost.t` (cost.mli below). `Ode` is interface-only
   (`modules_without_implementation ode`): a `val` declared here compiles and then fails at link
   time with "No implementation provided for ... Ode". *)
```

```ocaml
(* cost.mli *)
type t = { rhs_calls : int; jac_evals : int; factorizations : int; solves : int; newton_iterations : int }
val zero : t
val add : t -> t -> t
(* not `rhs`, `jac`: those would clash with the `problem` labels, and an unannotated `p.rhs`
   resolves to the last record defined *)
```

```ocaml
(* tol.mli *)
type t
val make : rtol:float -> atol:Vec.t -> t
val scalar : rtol:float -> atol:float -> n:int -> t

(** [weights tol y] is [1 / (rtol |y_i| + atol_i)]. *)
val weights : t -> Vec.t -> Vec.t

(** Weighted root-mean-square norm. *)
val wrms : weights:Vec.t -> Vec.t -> float
```

### 3.2 Methods

```ocaml
(* ode.mli, continued *)

(** Why a step could not be taken. [No_convergence] and [Singular] invite a shorter retry;
    [Nan] means the state or the right-hand side stopped being finite. *)
(* the driver maps it to [Fail.t]; this [Nan] and [Fail.Nan] are distinct constructors *)
type failure = No_convergence | Singular | Nan

(** One attempted step at order [order]: the new state, the local error estimate as a
    WRMS norm (acceptable when [<= 1]), the error norms the method would have at the
    neighbouring orders, and the history to continue from. *)
type 'history attempt = {
  y : Vec.t;
  err : float;
  err_vec : Vec.t;             (* the unweighted estimate; only Halving reads it, see 3.3 *)
  order : int;
  alternatives : (int * float) list;
  next : 'history;
  retry : 'history;            (* entry history with the refreshed cache, if the attempt is rejected *)
  cost : Cost.t;
}

module type Method = sig
  type history
  val start : problem -> history
  val order : history -> int

  (** [step p ~weights ~h ~order history at] advances one step. On failure the history comes
      back too, so the method can mark its Jacobian stale or drop its order. *)
  val step :
    problem -> weights:Vec.t -> h:float -> order:int -> history -> point ->
    (history attempt, failure * history) result
end

module type Embedded = sig
  include Method

  (** The method's own continuous extension on the last step, for output and events. *)
  val interpolate : history -> float -> Vec.t
end
```

`Method` carries the estimate in every `attempt`; [`Bdf1`](../../src/bdf1.ml), which has none today,
would need one (half the gap to the explicit Euler predictor, which [`Bdf2`](../../src/bdf2.ml)
already computes on its first step, at the price of one more rhs call per step). `Embedded` adds
interpolation. [`Stepper.fixed`](../../src/stepper.ml) keeps working with any `Method`, ignores
`err`, and has to supply the weights and the order that `step` now takes (a fixed order and weights
that reproduce today's `1e-10 (1 + |x|)` test keep its pins). On an error-test rejection the
attempt's `retry` history is the one to continue from; without it every rejection discards the
Jacobian and factorization built during the attempt.

### 3.3 Controllers

```ocaml
module type Controller = sig
  type t
  type stats

  type verdict = Accept of t | Reject of t | Give_up of Fail.t

  val init : Tol.t -> dt0:float -> dt_max:float -> max_rejects:int -> max_order:int -> t
  (* max_rejects stays: Halving and the public ?max_rejects of the driver need it *)

  (** The step and order to try next. *)
  val proposal : t -> float * int

  (** Judge an attempt of length [h]: accept it and plan the next step and order, reject it
      and plan a shorter retry, or give up. *)
  val decide : t -> h:float -> _ attempt -> verdict

  (** The method failed on a step of length [h]: shrink, and after too many, give up. *)
  val failed : t -> h:float -> failure -> verdict

  val stats : t -> stats
end
```

The controller sees only floats and ints (the method has already reduced vectors to norms), which
makes it the first candidate for the Rocq model in Section 4.6. That cannot be combined with keeping
[`Halving`](../../src/halving.ml)'s pins as they are: `Halving` tests `max_i |err_i| / (1 + |y_i|)
<= tol` on the error vector at the new state, which a method-side WRMS norm with start-of-step
weights does not reproduce, and the pinned lines (`accepted 1109, rejected 374`, `StepRejected 46`,
`max error 8.96e-08`, `|y1 - ref| = 7.2e-07`) depend on it. Options: (a) keep `Halving` and
[`Bdf2`](../../src/bdf2.ml) on the current contracts as frozen baselines beside the new ones (two
driver generations until deleting the old one is agreed); (b) give `attempt` the raw `err_vec` (as
above) so that `Halving` can be ported bit-identically while the new controllers ignore it; (c)
re-pin those lines, which breaks the rule that passing lines are never rewritten. Keeping Halving's
pins and changing the contract in one PR cannot both hold (D3, D13). `Dassl` (discrete rule), a
continuous `I` controller and `Pi` (filters) are added beside `Halving`.

### 3.4 Stage solver, linear solver, Jacobian provider

```ocaml
(* linsolve.mli: a linear solver for the iteration matrix; module types only, so it is listed in
   `modules_without_implementation` next to `ode` *)
module type Linsolve = sig
  type factored
  (** [factor ~gamma ~mass j] factors [mass - gamma j] (identity when [mass] is [None]). *)
  val factor : gamma:float -> mass:Mat.t option -> Mat.t -> (factored, [ `Singular ]) result
  val solve : factored -> Vec.t -> Vec.t
end

(* Dense : Linsolve and Band : Linsolve, chosen by Mat.structure. *)
```

```ocaml
(* stage.mli *)
type equation = { t : float; gamma : float; psi : Vec.t }

(** The Jacobian, the factorization, the gamma they were built with, and their age. *)
type cache
val fresh : cache
val stale : cache -> cache

val solve :
  (module L : Linsolve) -> problem -> weights:Vec.t -> cache -> equation -> guess:Vec.t ->
  (Vec.t * cache * Cost.t, failure * cache * Cost.t) result
```

```ocaml
(* jac.mli *)
(** [evaluate p t y ~fy] uses [p.jac] or finite differences grouped by [p.structure]. *)
val evaluate : problem -> float -> Vec.t -> fy:Vec.t -> Mat.t * Cost.t
```

### 3.5 Purity, mutation and the performance plan

Rule proposed for the kernels: a function may mutate only arrays it allocated itself in the same
call, and only until it returns them; arguments and arrays returned earlier (the factorization kept
in the cache) are never written; interfaces say nothing about it because nothing observable changes.
[`Linalg`](../../src/linalg.ml), `Band` and [`Jac`](../../src/jac.ml) are the three modules that use
the rule, and the rule is stated once in `docs/`. Everything above them stays as it is: immutable
`Vec.t`, records, lists.

Performance plan, in order of expected payoff:

1. Factorization reuse (2.3): O(n^3) -> O(n^2) per Newton iteration. Largest win, pure
   algorithmic.
2. In-place LU on a private copy: 4x on each factorization **[E]**.
3. Banded storage and grouped finite differences for method-of-lines problems: O(n) instead
   of O(n^3) per factorization and l+u+2 instead of n+1 rhs calls per Jacobian.
4. Allocation: [`Vec`](../../src/vec.ml) operations allocate one array each; at n <= 1000 that is
   noise next to a dense factorization (not next to banded or Krylov work, where everything is O(n),
   and every array above 256 words is allocated straight in the major heap). Keep them. Avoid
   closures inside the LU loops (plain `for` loops over `float array`), which the non-flambda
   compiler compiles well.
5. flambda: not in the development switch **[E]**; measure on a flambda switch before deciding
   anything; expect 10-20 % on the vector code, nothing on the LU.
6. `-unsafe` on the kernel modules only (`Linalg`, `Band`), after their property tests exist
   (bounds checks cost 10-20 % in such loops **[M]**). dune has no per-module flag field in the
   documentation found (ocaml/dune#3551 is the feature request **[?]**), so this means a separate library
   or explicit `Array.unsafe_get`.
7. Bigarray: no (2.4).

Bit-identical refactors: an in-place LU with the operation order of 2.0 is bit-identical to the
recursive solve **[R]**, so the LU swap can land as a pure refactor; the prototype LU, whose back
substitution sums differently, differs at the 1e-17 level **[E]**. Modified Newton, the
tolerance-linked test and the new estimates are the numerical changes. Policy in 4.4.

---

## 4. Verification and testing

### 4.1 Problems to add

The IVP Test Set (Mazzia, Magherini, Iavernaro, Bari; mirrored by the R package deTestSet
**[D, index read]**) and Hairer–Wanner II are the sources. Dimensions marked **[M]**.

| problem | type, n | why | reference source |
|---|---|---|---|
| Robertson to t = 1e11 with atol 1e-14 (y2 ~ 1e-10) | ODE, 3 | long-time stiffness, tiny components | to compute: scipy Radau 1e-13 cross-checked with BDF and LSODA (only t = 1e4 exists, Appendix B) |
| HIRES | ODE, 8 | moderately stiff, the standard first test | scipy Radau rtol 1e-13 (BDF and LSODA agree with it to 1.4e-13 absolute, Appendix B); the Test Set's published values were not available, and an agreement to 2e-17 with remembered Test Set digits is not recorded anywhere **[?]** |
| OREGO (Oregonator) | ODE, 3 | relaxation oscillations, large step ratios | Test Set |
| Pollution | ODE, 20 | chemistry, dense Jacobian | Test Set |
| E5 | ODE, 4 | needs atol ~ 1e-24 **[?]**, tests the weights | Test Set, Hairer–Wanner II |
| van der Pol, mu = 1000, t = 2000 | ODE, 2 | in the corpus; keep, add the reference state | Appendix B |
| Ring modulator (ODE form, Cs = 2e-12) | ODE, 15 | extremely stiff, oscillatory | Test Set **[D]** |
| Chemical Akzo Nobel | index-1 DAE, 6 (the Test Set's form **[D]**; y6 = Ks y1 y4 is explicit, so a 5-dimensional ODE after elimination **[M]**) | stiff chemistry, first mass-matrix test | Test Set |
| Brusselator MOL, N cells | ODE, 2N, banded (l = u = 2 interleaved) | larger n, banded Jacobian; n = 80 and 2000 | scipy Radau 1e-10 (Appendix B) for n = 80; n = 2000 to compute |
| Heat or Burgers MOL, n = 1000 | ODE, banded | linear and nonlinear large n | exact (heat), fine-grid reference (Burgers) |
| Medical Akzo Nobel | ODE, 400 | MOL with banded structure | Test Set |
| Prothero–Robinson y' = lambda (y - g(t)) + g'(t) | ODE, 1 | order reduction, stiff accuracy; exact | closed form |
| linear y' = A y with eigenvalues at 60 and 80 degrees from the negative real axis | ODE, 2 | A(alpha) canary for order selection | closed form |
| Transistor amplifier | index-1 DAE, 8 | first DAE (M9) | Test Set |
| Robertson in DAE form (y1 + y2 + y3 = 1) | index-1 DAE, 3 | mass matrix (M9) | same as ODE form |
| discontinuous forcing (step at t = 1) | ODE, 1 | rejection handling, restart | closed form |
| Lotka–Volterra, nonstiff | ODE, 2 | a stiff solver must not fail on nonstiff problems | scipy DOP853 |

### 4.2 References

Policy, with the evidence gathered:

1. Compute every reference with two independent solvers at tolerances at least 1e4 tighter than the
   test needs, and store it only if they agree to the digits used. Robertson at t = 1e4: scipy Radau
   (rtol 1e-13, atol 1e-20), BDF and LSODA (rtol 1e-12) agree to 3e-12 absolute (about; the bound
   recorded in [`test/refs.ml`](../../test/refs.ml) is 4e-12); HIRES to 1.4e-13; van der Pol at t =
   2000: Radau (1e-11) and LSODA agree to 1.6e-9 while **scipy's BDF failed at rtol 1e-11 (success =
   False, state wrong by 600)**, which is the argument for cross-checking **[E]**.
2. The Test Set's published reference values (computed with RADAU5/PSIDE at tight tolerances
   **[M]**) are a third source for its problems. A remembered HIRES value matching scipy Radau to
   2e-17 absolute could not be reproduced: no Test Set digits are stored, and the three scipy
   solvers themselves agree on HIRES only to 1.4e-13 absolute **[?]**.
3. References live in `test/refs.ml` with provenance in the comment (solver, version,
   tolerances, agreement), as now; never edited, only added.
4. Accuracy measure: the Test Set's scd = -log10 (max relative error at the end) and mescd
   (mixed, with atol) **[M]**; print scd in the work-precision tables.
5. Once Radau IIA exists (M10), vstiff's two methods at tight tolerance are a fourth and fifth
   source, and their agreement is a corpus line.

### 4.3 Work-precision diagrams as a quality gate

A `bench/` executable (the prototypes are a starting point) prints, per problem and per
tolerance, steps, rejections by cause, cost, order histogram, CPU time and scd against the
stored reference. A golden table is pinned per milestone; the gate is: scd not lower, and no
cost figure more than 10 % higher than the pinned one, with CPU time informational only (it
is machine dependent). The gate runs in CI on the small problems and nightly on the n >= 400
ones (the repository has no CI configuration today; a workflow is part of M1 **[R]**). Appendix B
is the first such table.

### 4.4 Expect tests and pins going forward

- Keep: stdout diffed against `.expected`, names per case, Guard and Report as the only
  effectful test modules, references never edited, no passing line weakened.
- Pins are per module: [`Halving`](../../src/halving.ml)'s step-control counts stay as they are
  because `Halving` stays; the new controllers get their own lines. A milestone never rewrites
  another module's pins.
- A numerical change (new estimate, new controller, modified Newton) lands in two PRs: the new
  module beside the old with a corpus line asserting agreement or a stated improvement, then the
  switch. The rule that passing lines are never rewritten means the switch cannot re-pin the old
  module's lines: either the old module stays as a frozen baseline or the re-pin is approved
  explicitly (D13). An in-place LU with the operation order of 2.0 is not such a change:
  one PR, no pin moves, checked by a random-system differential test.
- Bit-identical refactors keep the FMA discipline (never move products across a sum). The soak
  test's `identical: true` compares ten runs of the same binary (determinism); it cannot detect a
  refactor that changes bits. The detector is a hex dump of the final states compared between the
  old and the new build **[R]**.
- Pin order sequences: for variable order, a line prints the order histogram and the number
  of order changes on Robertson and OREGO; they move only with the selector.
- Pin the failure modes: Robertson at rtol 1e-3 with the corpus atol (`Ok`) and at
  rtol = atol = 1e-3 (`Ok` or `Error`, never a loop or a garbage `Ok`), the NaN wall, the
  blow-up, the Newton overflow (pinned by corpus line 24), the stuck-step cases of 1.3
  (criterion 4; pinned by corpus lines 25 to 27), all as `Ok`/`Error` lines.

### 4.5 Property tests and mutation testing

- Property tests with the standard library's `Random` (seeded, so the expect output is stable): the
  stage equation for order k reproduces [`Bdf2.coeffs`](../../src/bdf2.ml) at k = 2 and is exact on
  polynomials of degree <= k for random nodes; divided differences of a degree-d polynomial vanish
  above d; LU solves random well-conditioned systems to 1e-13 * cond; the controller's next step
  stays within the controller's clamp ([0.25 h, 2 h] for the DASSL rule) and never below the floor;
  interpolation reproduces the nodes. `qcheck-core` and `qcheck-alcotest` 0.91 are installed in the
  development switch (the `qcheck` meta package is not), but adding them is a dependency decision
  (D2); the hand-rolled generators are a hundred lines.
- Mutation testing: `mutaml` 0.3 is on opam but cannot be installed with OCaml 5.5: it requires
  `ppxlib < 0.36.0`, and every ppxlib below 0.36.0 requires OCaml < 5.4 (`opam show`) **[R]**. The
  plan is therefore a 40-line script that applies a dozen sed mutations to `bdf.ml`,
  [`stage.ml`](../../src/stage.ml) and [`linalg.ml`](../../src/linalg.ml) and requires the corpus to
  fail; it needs no dependency and runs in the CI. Run per milestone, not per commit.

### 4.6 How formal verification could fit

Investigated separately with Rocq and MathComp ([formal-verification.md](formal-verification.md),
the FV report); what this roadmap can do is shape the code so that the investigation has targets. In
order of plausibility:

1. **Coefficient lemmas, exact arithmetic**: the stage equation of order k is exact on
   polynomials of degree <= k through arbitrary distinct nodes; the error constant for k = 2
   is 2/9 at equal steps (3/22 for k = 3); divided differences are symmetric. These are
   statements about `poly` over a field in MathComp; the OCaml side must keep the coefficient
   generators as separate pure functions on the node positions (Section 2.1.1 does), so a
   Rocq mirror is the same function over an ordered field and the float version is tested
   against the mirror evaluated on Rocq's primitive floats (differential testing, not full
   verification).
2. **Stability**: A-stability of BDF1 and BDF2 is the root condition for rho(zeta) - z sigma(zeta)
   with Re z <= 0 (BDF2 has no stability function R(z); BDF1's is 1/(1 - z)), a statement about
   complex numbers (`R[i]` over a real closed field), not a real polynomial inequality.
   Zero-stability of variable-step BDF2 for ratios < 1 + sqrt 2 is shorter than a root-location
   argument: a1 + a0 = 1, so the step differences contract by the parasitic root -a0 per step, and
   -a0 < 1 for ratios below 1 + sqrt 2. The theorem with design value is that
   [`Halving`](../../src/halving.ml), [`Adaptive`](../../src/adaptive.ml) and
   [`Bdf2`](../../src/bdf2.ml) together never attempt a ratio above 2 in real arithmetic, a final
   step within the clock's resolution aside (FV report, T2 and T3; at f6d9b4e the largest observed
   ratio was exactly 2, and since dc3bb67 the snapping of steps allows rounding-level excursions
   above it).
3. **Driver and controller invariants**: with the controller as a function on floats and
   ints (3.3), its transition system can be modelled directly: the next proposal is never above
   dt_max (in real arithmetic; in binary64 a snapped step can exceed it by up to half an ulp of
   the new time) and never grows after a rejection. Termination was not provable at f6d9b4e: the loop did not
   terminate when every proposed step was below half an ulp of `t` (1.3, criterion 4 **[R]**).
   Since dc3bb67 it holds with at most (max_rejects + 1)(N + 1) attempts, N the number of floats
   in (t0, t_end]: every accepted non-final step strictly increases `t`, and the controller ends
   any run of rejections (FV report, T4).
4. **Linear algebra**: LU with partial pivoting solves the system in exact arithmetic when no
   pivot is zero; floating-point error bounds would need Flocq and are optional.
5. Tooling facts from the development switch's opam listing **[E, R]**: `rocq-core` and `coq` 9.2.0 are
   listed, but `rocq-runtime` 9.2.0 requires `dune < 3.24` and the development switch has dune
   3.24.2, so Rocq 9.2.0 cannot be installed there (Rocq 9.3.0 can, per the FV report, but not
   with MathComp). MathComp and Flocq live in the `rocq-released`
   repository, which is not configured in the development switch, and MathComp 2.6 additionally needs dune <= 3.23.1:
   the FV report's solver runs accept a separate switch (OCaml 5.5.0, dune 3.23.1, Rocq 9.2.0,
   MathComp 2.6.0). `coq-of-ocaml` exists as 2.5.3+4.12, +4.13 and +4.14 (`ocaml >= 4.14 & < 4.15`
   for the last), so it cannot read OCaml 5.5 modular explicits, and per the FV report it models
   floats as integers. The FV report recommends a hand-written Rocq mirror tested against the
   OCaml (its Option A) and advises against extraction (its Option B: Rocq's primitive floats have
   no `fma` while `ocamlopt` fuses five vstiff expressions on arm64, so extracted kernels would
   not be bit-identical). `gospel` 0.3.1 requires `ocaml <= 5.3.0` and `ortac-core` 0.8.0 pins
   `gospel = 0.3.1`, so neither installs on OCaml 5.5: specifications in `.mli` comments cost
   nothing, but nothing checks them yet.

---

## 5. Milestones

Effort is person-days for one developer who knows the code, shown as a first estimate / a
revised range; each milestone ends with its corpus lines, its pins and a work-precision table.

| ID | scope | depends on | acceptance (measurable) | risks | effort (first estimate / revised range) |
|---|---|---|---|---|---|
| **M1** | Bench and references: `bench/` executable (rhs calls counted with an [`Instrument`](../../src/instrument.ml)-wrapped rhs; costs as values need the contract change of M2), references for Robertson, HIRES, van der Pol, Brusselator-80 with provenance, first golden table, a CI workflow (none exists today) | — | the `lib + Halving` rows of Appendix B reproduced from the repository (an independent re-run matched them to the digit **[R]**) and the scipy rows regenerated by scripts kept in the repository; `dune test` unchanged and run by CI | none | 2 / 2-3 |
| **M2** | Tolerances and controller: `Tol`, WRMS weights, new `Controller` contract (3.3; D3, D13), [`Halving`](../../src/halving.ml) ported through `err_vec` or kept as a frozen baseline, `Dassl` and a continuous `I` controller added, the clock rules already in the driver (snapped steps, `Ode.Too_small`, [`Clock.resolution`](../../src/clock.ml); dc3bb67) carried over to the new contract, an initial step-size heuristic in place of `dt0 = 1e-6 (t_end - t0)` (Hairer–Nørsett–Wanner I, II.4 **[M]**), driver updated | M1 | every existing corpus and soak line unchanged (D13); the I controller at rtol 1e-6: rhs calls <= 50k / 115k / 600k on Robertson / HIRES / van der Pol (prototype 49.1k / 110.9k / 592k: the margins are 1-4 %, fix them after the first run) and rejections < 2 % of accepted steps at rtol <= 1e-5; the discrete `Dassl` rule is judged at equal achieved error (-7 % to -29 % against Halving; at the same tolerance it costs 89.6k / 181k / 828k); the clock cases of 1.3 keep their pinned results (corpus lines 25 to 27) | the contract change touches every test: split it into three PRs (D3); Halving's pins cannot move | 3 / 4-6 |
| **M3** | Linear algebra and modified Newton: in-place dense LU behind a pure interface (first PR, a bit-identical refactor), `Linsolve`, `Stage.cache`, Jacobian and factorization reuse, DASSL/CVODE convergence and refresh rules, failure returns history | M2 | LU: bit-identical to [`Linalg.solve`](../../src/linalg.ml) on the corpus and on 2e4 random systems, about 4x faster; with the I controller of M2, rhs calls at rtol 1e-6 <= 21k / 20k / 240k (carried-J prototype 20.1k / 18.9k / 231k); Brusselator-80 at rtol 1e-4 <= 7k rhs calls and at most one factorization per accepted step; average Newton iterations per step <= 3 (prototype 11 with the fixed test; the 3-4 iteration limits of CVODE and DASSL make 2-4 plausible, unmeasured **[M]**); Robertson at rtol 1e-3 with the corpus atol finishes, at rtol = atol = 1e-3 it returns `Ok` or `Error`; CPU informational (3.0 s idle at f6d9b4e on the Brusselator) | the tolerance-linked convergence test changes achieved errors slightly; two-PR policy (4.4) | 4 / 6-8 |
| **M4** | BDF2 with the divided-difference error estimate (one solve per step), history as a point list, quadratic predictor | M3 | steps at rtol 1e-6 on Robertson <= 500 (from 2531; prototype 320); rhs calls at achieved error 1e-6 <= 9k / 6k / 120k (extrapolated from the rtol 1e-8 runs, whose errors are 1.9e-6 / 1.6e-6 / 8.9e-6: the prototype never reached 1e-6, and van der Pol extrapolates to about 112k); achieved error <= 100x tol (the prototype gives 6x-890x: a safety factor, or relax the criterion to scipy's level); A-stability canary lines unchanged | estimate constants for variable steps: unit-test against the 2/9 limit; the estimate is less conservative than the backward Euler gap, which is what lets Robertson wander negative at rtol = atol = 1e-3 (2.3) | 3 / 3-5 |
| **M5** | Variable order 1-5: stage equation of order k from the point history, estimates at k-1 and k+1, DASSL order selection, startup phase, max order option, A(alpha) canaries | M4 | rhs calls at achieved error 1e-6 within 2x of scipy BDF's true counts (669 / 844 / 5400 at that error): <= 1.3k / 1.7k / 11k (36x-1250x at f6d9b4e); order histogram pins on Robertson and OREGO; the 60 and 80 degree canaries pass; the startup ramp is capped below the zero-stability limit of the order in use, or shown harmless by a corpus line (2.1.2) | order-selection tuning is the long tail; keep DASSL's rules verbatim (2.1.2) before tuning | 8 / 12-20 |
| **M6** | Jacobian providers and banded storage: `problem.jac`, `Mat.structure`, grouped finite differences, `Band` LU | M3 | Brusselator with n = 2000 at rtol 1e-6: l + u + 1 = 5 rhs calls per Jacobian beyond f(t, y) (pinned) and O(n) work per factorization; a user Jacobian that is [`Jac.forward`](../../src/jac.ml) reproduces the finite-difference run bit for bit, an analytic one agrees to 1e-6 (the finite-difference Jacobian even of a linear problem is not bit-identical to the analytic one); CPU under 1 s informational | none | 4 / 4-6 |
| **M7** | Outputs: `interpolate`, output at times, step sequence API, events with Illinois | M5 | output at 100 times costs no extra steps (pinned step count); an event on van der Pol's zero crossing: the root of the interpolant located to 1e-10 in t (the true event time is only as accurate as the solution); discontinuity test restarts at order 1 | API shape (D7) | 4 / 4-6 |
| **M8** | Corpus expansion and the work-precision gate: OREGO, Pollution, E5, Ring modulator, Medical Akzo, heat/Burgers MOL, Prothero–Robinson, nonstiff (Chemical Akzo Nobel is an index-1 DAE: M9); property tests; mutation script; nightly gate | M5, M6 | every robustness criterion of 1.3 holds on the ODE set; mutation script kills all seeded mutants | the n = 400 problems are slow until M6; every problem needs a reference from two solvers | 3 / 6-9 |
| **M9** | Mass matrix and index-1 DAEs: `problem.mass`, stage equation with M, consistent-initial-condition helper, transistor amplifier, Chemical Akzo Nobel and Robertson-DAE in the corpus | M5, M6 | transistor amplifier scd >= 5 at rtol 1e-7; Robertson-DAE matches the ODE form to 1e-9 at rtol 1e-10 | error test on algebraic variables | 6 / 6-9 |
| **M10** | Radau IIA (order 5) as a second `Embedded` method with the PI controller and collocation dense output | M3, M7 | van der Pol at achieved error 1e-8 within 2x of scipy Radau's true count (8.1k at that error): <= 16k; BDF-vs-Radau agreement line at 1e-10 on Robertson and HIRES | complex LU, the transformation constants (test against RADAU5 numbers) | 10 / 12-18 |
| **M11** | Large n: matrix-free GMRES with user preconditioner, Jacobian-vector products by differences | M6 | heat MOL with n = 1e5 at rtol 1e-6 finishes with pinned iteration counts (target under 10 s, informational) | only if a user needs it | 5 / 5-8 |
| **M12** | Verification hooks: coefficient generators as separable pure functions, contracts as `.mli` comments, a hand-written Rocq mirror (not extraction) with a canonical-source tripwire and differential tests (FV report, Option A) | M5, FV report | the mirror evaluated on Rocq primitive floats agrees with `Bdf` coefficients on 1e4 random node sets to 1e-14 relative; the tripwire is green in CI. The proofs (FV pilot T1-T3: 1-2 weeks for an expert, 3-6 for a newcomer; T4 afterwards, its code fix being in dc3bb67) are not in the 3 days | a separate switch (dune < 3.24); toolchain churn | 3 plus the proof work / 3 plus the proof work |
| **M13** | Release: odoc (not installed in the development switch), opam package (GPL-3.0-only), API freeze, README with the work-precision tables | M7, M8, and whichever of M9-M11 ship in the first release | `opam install vstiff` from the repository; docs build without warnings | none | 2 / 2 |

Ordering rationale: M1-M4 are 12 days (revised: 15-22) for a 9x to 36x reduction in rhs calls at
the error Halving reaches at rtol 1e-6 on the three small problems **[R]**, and for the removal
of the rejection storms, the per-iteration Jacobian and the second solve per step; M5 is the only
large step and is where the design decision D1 bites; M6 is what makes n > 100 usable;
everything after is generality. Total: about 57 person-days as first estimated; the revised
range is 69-103 (M2, M3, M5, M8 and M10 carry most of the difference: pin preservation, the
CVODE/DASSL logic and its tuning, ten problems each needing a reference, and a Radau
implementation with complex linear algebra), plus the proof work. None of the estimates is
calibrated against measured effort on this code base **[?]**.

---

## 6. Open questions and decisions

| ID | question | recommendation |
|---|---|---|
| **D1** | BDF representation: variable-coefficient over the point history (current style), fixed-leading-coefficient Nordsieck (CVODE), or divided differences (DASSL)? | Variable-coefficient over the point history: least distance from the code, best readability, the widest zero-stability limits of the two forms computed (2.1), best fit for the proof work; revisit only if factorization reuse measured after M5 is below 80 % of steps. |
| **D2** | Stay standard-library only, including the tests? | Yes for the library, and yes for the tests (hand-rolled property generators on `Random`); allow `sundialsml` or scipy only as external oracles that produce the stored references, outside the build. |
| **D3** | Change the `Controller` contract to `decide` returning a verdict with the next state, and `Method.step` to return the history on failure? | Yes in substance, not as one PR with the pins moved: the change bundles three separable ones (a verdict that sees the error vector, history on failure, costs as values) and Halving's pins cannot move (3.3). Three PRs, each bit-identical on the existing corpus; Halving ported through `err_vec` or frozen (D13). |
| **D4** | Is local mutation of freshly allocated arrays, frozen once returned, in [`Linalg`](../../src/linalg.ml), `Band` and [`Jac`](../../src/jac.ml) compatible with "effects quarantined"? | Yes, with the rule written down in `docs/` and the three modules named; the interfaces stay pure. The soak test does not guard it (it compares ten runs of one binary): a random-system differential test and a hex dump of the corpus states do. |
| **D5** | Step policy for BDF: discrete (DASSL hold/double/shrink) or smooth (PI filters)? | Discrete for BDF (for the stability limits and factorization reuse; no rhs-call advantage was measured, 2.2), PI34 or H211b for Radau. |
| **D6** | Default and maximum order. | Default 5, user-settable 1-5, no BDF6. |
| **D7** | Output API: closures (`dense : float -> Vec.t`) or data (`Dense.t` with `eval`) and a lazy `Seq.t` of steps? | Data plus `Seq.t`; closures as one-line helpers on top. |
| **D8** | DAEs in scope? | Index 1 and mass matrices after M6; no higher index. |
| **D9** | Cost accounting as values (drop [`Instrument`](../../src/instrument.ml) from the library)? | Yes; the counts become deterministic pins, validated against an `Instrument`-wrapped rhs in the tests. |
| **D10** | Newton inside the integrator: keep the Armijo line search? | No; modified Newton with the tolerance-linked rate test and step cuts, as every production code; [`Newton.solve`](../../src/newton.ml) stays for initial conditions and tests. |
| **D11** | Should the n = 400 and n = 2000 problems run in the default `dune test`? | No: nightly gate. The default `dune test` took about 14 s at f6d9b4e (corpus 1.4 s, soak 12.4 s, idle), so the budget has to be stated for the whole run, not the corpus alone. |
| **D12** | Priority of Radau IIA vs DAEs after M8. | DAEs first (M9) if circuit or chemistry users exist, Radau first (M10) if tight tolerances or the cross-check matter more. |
| **D13** | Passing pins of a superseded algorithm ([`Bdf2`](../../src/bdf2.ml) + [`Halving`](../../src/halving.ml) + `Newton`/[`Stage`](../../src/stage.ml): `accepted 1109, rejected 374`, `StepRejected 46`, `max error 8.96e-08`, the Robertson error 7.2e-07): keep the old modules as a frozen baseline until the new stack beats them on every corpus line, or re-pin them in the PR that supersedes them? | Freeze (a `legacy` sublibrary, no edits), compare against it in differential corpus lines, delete in one PR when we decide the new stack has superseded it; the rule that passing lines are never rewritten then stands. |
| **D14** | Default atol when the caller gives one tolerance (today rtol = atol = tol, which lets Robertson's y2 go negative at 1e-3)? | Require atol explicitly in the new API (`Tol.make`), keep the scalar `~tol` only on the legacy path, and document that atol must be well below the smallest component that matters. |

---

## Appendix A. Provenance of the measurements

The measurements in this document come from throwaway prototype scripts that are not kept in the repository, so the figures are not reproducible from the repository alone.

Caveat on the CPU columns: `Sys.time` syscalls around every solve distort the n <= 8 timings
by up to 2x; rhs and Jacobian counts are the reliable columns there. The n = 80 timings are
unaffected by that, but every CPU cell was probably measured while other jobs were running: idle, the
n = 80 runs take 3.0 s (lib + Halving), 2.9 s (mn + I(1/2); 4.2 s in the back-to-back comparison
with lib + I(1/2) at 3.1 s), 0.42 s and 1.3 s (pc + I(1/3) at rtol 1e-4 and 1e-6), and a solve
costs 0.5-0.8 ms, not 1.2 ms **[R]**.

## Appendix B. Measured work-precision data [E]

Mixed error = max_i |y_i - ref_i| / (1 + |ref_i|) at t_end. vstiff rows: BDF2 in all cases. `lib` =
library Newton (Jacobian every iteration), `mn` = modified Newton with the Jacobian in the history,
`pc` = one solve per step with the divided-difference estimate; [`Halving`](../../src/halving.ml) =
the repository controller, `I(e)` = error-based with exponent e, `PI` = PI34. scipy rows use
finite-difference Jacobians; for BDF and Radau the `nfev` shown does not include the Jacobian
columns (it does for LSODA), so the scipy BDF and Radau rhs cells understate the work by the amounts
in the last table of this appendix.

Robertson, t = 1e4 (reference y1 = 0.1073004285378, agrees with [`test/refs.ml`](../../test/refs.ml)):

| config | rtol 1e-4: steps (rej), rhs, err | rtol 1e-6 | rtol 1e-8 |
|---|---|---|---|
| lib + Halving | 261 (71), 11 195, 6.9e-5 | 2531 (831), 84 267, 6.5e-7 | 25 298 (8423), 684 580, 6.1e-9 |
| lib + I(1/2) | 205 (2), 7063, 9.0e-5 | 1957 (3), 49 055, 9.5e-7 | 19 479 (8), 389 871, 9.6e-9 |
| mn + I(1/2) | 211 (5), 2875, 9.1e-5, 21 J | 1958 (4), 20 069, 9.5e-7, 16 J | 19 472 (8), 145 043, 1.3e-8, 11 J |
| pc + I(1/3) | 87 (5), 651, 6.6e-4 | 320 (4), 1890, 3.8e-5 | 1409 (7), 6129, 1.9e-6 |
| pc + PI34 | 127 (3), 856, 3.2e-4 | 428 (3), 2540, 2.2e-5 | 1821 (6), 8539, 1.1e-6 |
| scipy BDF | 85, 243, 8.7e-5 | 176, 548, 2.6e-6 | 346, 1080, 3.9e-8 |
| scipy Radau | 37, 393, 1.3e-6 | 84, 734, 6.4e-9 | 243, 1883, 1.4e-10 |
| scipy LSODA | 114, 238, 1.3e-4 | 245, 498, 8.4e-7 | 467, 902, 3.2e-9 |

HIRES, n = 8, t = 321.8122:

| config | rtol 1e-4 | rtol 1e-6 | rtol 1e-8 |
|---|---|---|---|
| lib + Halving | 327 (95), 22 062, 9.5e-5 | 3212 (1059), 187 292, 1.0e-6 | 32 587 (10 854), 1 773 376, 9.9e-9 |
| lib + I(1/2) | 252 (16), 14 301, 1.7e-4 | 2470 (1), 110 851, 1.6e-6 | 24 770 (3), 1 019 113, 1.5e-8 |
| mn + I(1/2) | 298 (31), 4352, 6.0e-5, 66 J | 2470 (1), 18 939, 1.6e-6, 19 J | 24 770 (3), 151 103, 1.4e-8, 9 J |
| pc + I(1/3) | 92 (8), 922, 7.0e-4 | 370 (3), 1721, 3.7e-5 | 1672 (4), 4455, 1.6e-6 |
| pc + PI34 | 127 (6), 998, 6.9e-4 | 482 (0), 1924, 2.6e-5 | 2150 (3), 5474, 9.9e-7 |
| scipy BDF | 80, 242, 1.3e-3 | 164, 450, 1.0e-5 | 308, 824, 1.0e-7 |
| scipy Radau | 31, 399, 8.4e-6 | 81, 803, 6.5e-8 | 240, 2027, 2.5e-10 |
| scipy LSODA | 110, 374, 2.9e-4 | 267, 701, 5.3e-6 | 495, 1225, 9.3e-8 |

van der Pol, mu = 1000, t = 2000 (reference y = (1.706167732170524, -8.928097010247538e-4)):

| config | rtol 1e-4 | rtol 1e-6 | rtol 1e-8 |
|---|---|---|---|
| lib + Halving | 4145 (1349), 116 202, 2.6e-4 | 39 041 (12 984), 978 933, 2.8e-6 | 420 057 (139 997), 9 015 709, 3.1e-8 |
| lib + I(1/2) | 3148 (58), 68 325, 4.7e-4 | 31 572 (4), 592 084, 4.6e-6 | 316 243 (6), 5 094 962, 4.5e-8 |
| mn + I(1/2) | 3179 (66), 37 116, 4.6e-4, 161 J | 31 572 (4), 231 437, 4.6e-6, 42 J | 316 241 (6), 1 668 404, 3.8e-8, 22 J |
| pc + I(1/3) | 716 (112), 6299, 4.5e-3 | 3267 (14), 13 760, 2.1e-4 | 15 203 (20), 45 389, 8.9e-6 |
| pc + PI34 | 883 (63), 6306, 3.0e-3 | 4146 (5), 16 644, 1.4e-4 | 19 424 (13), 46 147, 4.2e-6 |
| scipy BDF | 428, 1394, 2.3e-3 | 847, 2620, 4.2e-5 | 1846, 5472, 7.1e-7 |
| scipy Radau | 218, 2084, 1.5e-5 | 616, 5167, 1.4e-7 | 1853, 14 443, 6.9e-10 |
| scipy LSODA | 502, 1149, 2.3e-3 | 880, 1726, 5.0e-5 | 1596, 2948, 8.4e-7 |

Brusselator MOL, 40 cells, n = 80, alpha = 1/50, t = 10 (reference scipy Radau rtol 1e-10,
LSODA agrees to 2e-9):

| config | rtol | steps (rej) | rhs | Jacobians | Newton it/step | CPU | error |
|---|---|---|---|---|---|---|---|
| lib + Halving | 1e-4 | 737 (230) | 453 298 | (every iteration) | — | 7.9 s | 6.0e-5 |
| lib + I(1/2) | 1e-4 | 601 (0) | 292 495 | (every iteration) | about 6 | 2.9-3.3 s | 7.5e-5 |
| mn + I(1/2) | 1e-4 | 601 (0) | 6932 | 3 | 11.1 | 7.6 s (98 % in solves) | 7.5e-5 |
| pc + I(1/3) | 1e-4 | 151 (7) | 1460 | 5 | 6.7 | 1.2 s (98 %) | 1.2e-3 |
| pc + I(1/3) | 1e-6 | 669 (1) | 3234 | 1 | 4.7 | 3.6 s (98 %) | 6.2e-5 |
| scipy BDF | 1e-4 / 1e-6 | 94 / 197 | 283 / 538 | 5 / 2 | — | 29 / 43 ms | 1.7e-4 / 4.4e-6 |
| scipy Radau | 1e-4 / 1e-6 | 42 / 123 | 359 / 938 | 16 / 22 | — | 51 / 91 ms | 3.7e-6 / 1.2e-8 |
| scipy LSODA | 1e-4 / 1e-6 | 148 / 324 | 1368 / 1618 | 14 / 13 | — | 26 / 35 ms | 6.5e-5 / 1.8e-6 |

The high Newton iteration counts of the prototypes come from the fixed 1e-10 convergence test
on the step; the tolerance-linked rate test of 2.3 is expected to bring them to 2-4 **[M, from
CVODE's and DASSL's limits]**, which is part of M3's acceptance. The CPU cells above were measured
earlier, probably under load (Appendix A).

Controllers at rtol = atol = 1e-6, steps (rejected), rhs calls, error (`lib`
estimate = backward Euler gap with the library Newton, `pc` = divided-difference estimate with
carried Jacobian; D(e) = DASSL's discrete rule with exponent e):

| config | Robertson | HIRES | van der Pol |
|---|---|---|---|
| lib + Halving | 2531 (831), 84 267, 6.5e-7 | 3212 (1059), 187 292, 1.0e-6 | 39 041 (12 984), 978 933, 2.8e-6 |
| lib + I(1/2) | 1957 (3), 49 055, 9.5e-7 | 2470 (1), 110 851, 1.6e-6 | 31 572 (4), 592 084, 4.6e-6 |
| lib + D(1/2) | 3580 (3), 89 630, 3.3e-7 | 4399 (0), 180 941, 8.8e-7 | 44 966 (4), 828 449, 2.4e-6 |
| pc + I(1/3) | 320 (4), 1890, 3.8e-5 | 370 (3), 1721, 3.7e-5 | 3267 (14), 13 760, 2.1e-4 |
| pc + D(1/3) | 503 (3), 2793, 1.6e-5 | 576 (0), 2012, 2.5e-5 | 4157 (4), 16 783, 1.4e-4 |

scipy rhs evaluations counted with a wrapper (reported `nfev` in brackets):

| config | Robertson 1e-4 / 1e-6 / 1e-8 | HIRES 1e-4 / 1e-6 / 1e-8 | van der Pol 1e-4 / 1e-6 / 1e-8 | Brusselator 1e-4 / 1e-6 |
|---|---|---|---|---|
| scipy BDF | 297 (243) / 589 (548) / 1129 (1080) | 387 (242) / 649 (450) / 1087 (824) | 1535 (1394) / 2791 (2620) / 5706 (5472) | 688 (283) / 700 (538) |
| scipy Radau | 473 (393) / 846 (734) / 2031 (1883) | 536 (399) / 1028 (803) / 2509 (2027) | 2210 (2084) / 5411 (5167) / 15 003 (14 443) | 1639 (359) / 2698 (938) |

## Appendix C. Sources

- Hairer, Wanner, *Solving Ordinary Differential Equations II*, 2nd ed., Springer 1996:
  IV.8 (Radau IIA, simplified Newton), V.1-V.2 (BDF stability), V.5 (variable-step BDF
  implementations), IV.10 and the test problems. **[M]**
- Hairer, Norsett, Wanner, *Solving ODEs I*, III.5 (variable step multistep). **[M]**
- Brenan, Campbell, Petzold, *Numerical Solution of Initial-Value Problems in
  Differential-Algebraic Equations*, SIAM Classics 1996, ch. 5 (DASSL implementation; chapter
  confirmed by search) **[M]**; DDASSL source at netlib.org/ode/ddassl.f **[D]**.
- Hindmarsh et al., *CVODE User Guide*, "Mathematical considerations",
  sundials.readthedocs.io, and the sources `src/cvode/cvode.c`, `cvode_impl.h`, `cvode_ls.c`
  (LLNL/sundials, main) **[D]**.
- Radhakrishnan, Hindmarsh, *Description and Use of LSODE*, NASA RP-1327 / LLNL 1993. **[M]**
- Jackson, Sacks-Davis, "An alternative implementation of variable step-size multistep formulas
  for stiff ODEs", ACM TOMS 6(3), 1980 (fixed leading coefficient). **[M]**
- Shampine, Reichelt, "The MATLAB ODE Suite", SIAM J. Sci. Comput. 18(1), 1997 (NDFs,
  quasi-constant steps; scipy's BDF follows it). **[M; the PDF is behind a login]**
- Soderlind, "Automatic control and adaptive time-stepping", Numer. Algorithms 31, 2002;
  "Digital filters in adaptive time-stepping", ACM TOMS 29(1), 2003; coefficients as
  implemented in PETSc `src/ts/adapt/impls/dsp/adaptdsp.c` **[D]**.
- Gustafsson, "Control-theoretic techniques for stepsize selection in implicit Runge–Kutta
  methods", ACM TOMS 20(4), 1994. **[M]**
- Akrivis, Katsoprinakis, "Maximum angles of A(theta)-stability of backward difference
  formulae", BIT 60 (2020) 93-99 (cs.uoi.gr/~akrivis/AK2a.pdf; angles 86.03, 73.35, 51.84,
  17.84). **[D]**
- Li, Liao, "Stability of variable-step BDF2 and BDF3 methods", arXiv 2201.00527 (BDF3: "almost
  all adjacent step ratios less than 2.553"; stability for dissipative problems, not
  zero-stability; the abstract does not mention Grigorieff). **[D, abstract]**
- Grigorieff, "Stability of multistep-methods on variable grids", Numer. Math. 42 (1983)
  359-377 (zero-stability of variable-step BDF2 for ratios < 1 + sqrt 2, sharp). **[D: confirmed
  by literature search; paper not read]**
- Mazzia, Magherini, Iavernaro, *Test Set for Initial Value Problem Solvers*, Bari; R package
  deTestSet (Soetaert, Mazzia) **[D, package index]**. Chemical Akzo Nobel (index-1 DAE, 6
  equations), Medical Akzo Nobel (ODE, 400) and the ring modulator (ODE of 15 equations for
  Cs = 2e-12, index-2 DAE for Cs = 0) were checked on search-result extracts of the Test Set
  pages and report; archimede.uniba.it was not reachable.
- Curtis, Powell, Reid, "On the estimation of sparse Jacobian matrices", J. Inst. Math. Appl.
  13, 1974. **[M]**
