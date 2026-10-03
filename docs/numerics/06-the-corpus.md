# 6. The test corpus

Chapter 6 of 6 in the numerical-methods track. Previous: [5. Error estimates and step-size control](05-step-control.md). Index and reading order: [docs/README.md](../README.md). Terms and symbols: [glossary](../glossary.md).

**Summary.** The corpus is a table of 47 cases whose printed lines are pinned in [`test/corpus.expected`](../../test/corpus.expected). This chapter says, line by line, what each one proves, what a failure points at and what it cannot see: nine lines that climb from Newton to a stiff chemical system, and thirty-eight regression pins that each guard one rule. It introduces the four ODEs behind them (the canary, the logistic equation, van der Pol and Robertson), the external reference value and why it is never edited, the four soak lines, a table of deliberate bugs with the line that notices each, and an honest list of what the corpus does not catch, each item a possible contribution.

Running the tests, reading a failure, adding a case and the policy for changing expectations are in [docs/testing.md](../testing.md), which this chapter does not repeat. The mathematics comes from chapters [1](01-odes-and-stiffness.md) to [5](05-step-control.md). Snippets run in the probe project (set up on day 1 from Setup in [exercises.md](../exercises.md), or "Running the snippets" in chapter 1); those that use `Refs` also need a link to `test/refs.ml`, made from the repository root next to the one for `problems.ml`:

```sh
ln -sfn "$(pwd)/test/refs.ml" ../vstiff-scratch/p/refs.ml
```

OCaml notes for them ([docs/ocaml.md](../ocaml.md) has the details): `(module Bdf2)` passes a module as an argument; `{ Ode.rhs = ...; t0 = ... }` builds a record, the `Ode.` prefix saying which module its fields belong to; `(s : T)` is a type annotation that tells OCaml which record type `s` has, so that `s.stats` can be read.

## 1. What the corpus is for

[`test/corpus.ml`](../../test/corpus.ml) is a table of cases, each a name and a function that computes the text printed after it. `Report.lines` prints `name: text` for each, in order, and `dune runtest` compares the output with `test/corpus.expected`. Several cases print a boolean (`true` means the criterion held), so the criterion lives in the code and the expected file pins its outcome; others print digits, which pin more. The first nine lines, from six development steps, climb the code layer by layer, the first two layers in the kernel (`Numerics`) and the rest in the solver (`Vstiff`):

```
Newton + Linalg  →  Jac  →  Stage + Bdf1  →  Bdf2  →  Adaptive + Halving  →  a stiff chemical system
```

The other 38 are regression pins, each added for a specific mistake (section 9 lists what catches what). The lines are numbered in the order of `corpus.expected`, which is the order of the named lists in `corpus.ml`: `newton` (1 to 3), `jacobian` (4 and 5), `backward_euler` (6), `bdf2_order` (7), `van_der_pol` (8), `robertson` (9), then the pins `stiff_canary` (10), `orientation` (11), `pivoting` (12 and 13), `give_up` (14 to 20), `step_control` (21 and 22), `robertson_accuracy` (23), `newton_overflow` (24), `clock` (25 to 27), `too_small` (28), `fixed_clock` (29 to 31), `zero_floor` (32), `arguments` (33 to 37), `newton_guards` (38 to 40), `jacobian_step` (41), `adaptive_time` (42), `adaptive_canary` (43) and `start_and_limits` (44 to 47). Two rules of thumb: when several lines differ, start with the first, which involves the fewest layers; and a case earns its place by pinning something the others do not.

The problems live in [`test/problems.ml`](../../test/problems.ml). Each module exposes `rhs`, `y0`, a `problem` record and, where the exact solution has a formula, `exact`. None supplies a Jacobian: the solver must work for any black-box right-hand side. All four right-hand sides ignore their time argument (`rhs _t y`); the cases that need one that depends on `t` define it inline (lines 15, 29, 30, 32 and 42, section 8). Van der Pol, the two Robertson lines and the canary through the adaptive driver (8, 9, 23 and 43), like the termination and argument cases of section 8 (lines 14 to 20, 25 to 28 and 32 to 37), run on the call budget of `Guard`, so a bug that makes one of them crawl prints a line within seconds instead of stalling the run.

| Problem | Components | Time scales | Corpus interval | Stiff? | Closed form |
|---------|-----------|-------------|-----------------|--------|-------------|
| Canary | 3 | rates 1, 100, 1e4 | `[0, 1]` | yes, ratio `1e4` | yes |
| Logistic | 1 | about 1 | `[0, 5]` | no | yes |
| Van der Pol, `μ = 1000` | 2 | slow drifts, jumps on `1/μ = 1e-3` | `[0, 2000]` | yes | no |
| Robertson | 3 | rate constants 0.04 to 3e7 | `[0, 1e4]` | yes | no |

## 2. Newton and Linalg: lines 1 to 3, 12, 13, 24, 38 to 40 and 45 to 47

```
newton linear 2d: Ok [2.000000000000; 3.000000000000]
newton quadratic x0=1: Ok [1.414213562373]
newton quadratic x0=0: Error Diverged
linalg zero leading pivot: [1.000000000000; 2.000000000000; 3.000000000000]
linalg tiny leading pivot: [1.000000000000; 1.000000000000]
newton step that converges into an overflow: Error Nan
newton atan x0=1.5 (needs damping): Ok [0.000000000000]
newton residual nan at the start: Error Nan
newton step that overflows in the linear solve: Error Diverged
newton x^2 = 2 from x0=1e-3 (deep line search): Ok [1.414213562373]
newton x^3 = 0 from x0=1 (needs more than 50 iterations): Error Diverged
newton log x = 0 from x0=3 (full step leaves the domain): Ok [1.000000000000]
```

Lines 1 to 3 test `Newton.solve` with a hand-written Jacobian, `Vec` and `Linalg`, and nothing above them. Line 1 is `G(x) = A x - b` with `A = [[3, 1], [1, 2]]`, `b = (9, 8)`, from `(0, 0)`: one Newton step ([chapter 2](02-newton.md)), printed with 12 decimals. Line 2 is `x² - 2` from 1, whose iterates are in chapter 2; 12 printed decimals make the line sensitive to Newton's tolerance. Line 3 is the same function from 0, where the 1×1 Jacobian is singular: `Linalg.solve` returns `None` and the failure comes back as the value `Error Diverged`. What they do not exercise: `A` is symmetric with its largest entry already in the corner, so neither orientation nor a row exchange is tested, and no step needs shortening (line 38 below does).

Lines 12 and 13 close the row-exchange gap, both worked by hand in chapter 2: a zero in the corner that only a row exchange gets past (without it the answer is `None`), and a tiny corner entry that elimination without the exchange turns into a first unknown of `0` instead of `1`.

Line 24 closes a gap in the stopping test of chapter 2, which is relative, `|dx|_inf <= 1e-10 (1 + |x|_inf)`: next to the largest float, `max_float` (about `1.8e308`), a step of up to about `1.8e298` counts as tiny, and adding it can overflow. The case solves `x - max_float = 0` from `max_float - 1e297` and gives Newton a Jacobian of `0.25` where the true one is `1`, so the step comes out four times too long, about `4e297`. It passes the test, and `x + dx` lies beyond `max_float`, which is infinity. A solver must not call that a root, so `Newton.solve` checks the sum and returns `Error Nan`; without the check the line prints `Ok [inf]`.

Lines 38 to 40 pin safeguards of [chapter 2](02-newton.md) that no integration step needs. Line 38 is `atan` from 1.5, just past the point where plain Newton fails: the full step goes from 1.5 to `-1.69`, where `|G|` is `1.04` against the `0.98` it started from, and the plain iterates (`-1.69, 2.32, -5.11, 32.3, ...`) run away until `x²` overflows, the Jacobian is exactly 0 and the solve gives up. With damping the half step passes and the root comes out as `0`; without damping, or with a line search that accepts any finite trial point, the line prints `Error Diverged`. A minimum damping of `1/2` still passes (the half step is all that is needed), so line 45 pins the depth of the search. Line 39 pins a name: a residual that is `nan` at the start is `Error Nan`, not the `Error Diverged` of a stall. Line 40 pins the other: `x - 1` with a Jacobian of `1e-320`, a float that is not 0, so `Linalg.solve` succeeds and returns the step `1 / 1e-320`, which overflows to infinity. A step that cannot be taken is `Error Diverged`. The line pins that name and not the check on the step, which is redundant: without it the line search meets the infinite trial point and gives up with the same `Diverged`.

Lines 45 to 47 pin three limits of the same chapter. Line 45 is `x² - 2` from `1e-3`: Newton's step is about `+1000`, the new point needs `|G| <= 2`, so only `λ <= 2e-3` passes, and the first factor that small is `1/512`, the tenth trial ([chapter 2](02-newton.md), section 5); from there Newton converges to `1.414213562373`. With a minimum damping of `1/2` the line prints `Error Diverged`. Line 46 is `x³` from 1, a triple root where Newton converges only linearly, the error shrinking by `2/3` per iteration: it takes 55 iterations to meet the stopping test and the limit is 50, so the line prints `Error Diverged`; with a limit of 100 it prints `Ok [0.000000000138]`, which is `(2/3)^56` (chapter 2, section 6). Line 47 is `log x` from 3: the full step, `-3 ln 3`, lands at `-0.30`, outside the domain of `log`, where `Float.log` is `nan`. The line search must refuse that trial point and take the half step, `1.35`, from which Newton converges to `1.000000000000`; a line search whose test lets a `nan` residual pass (`not (norm > bound)` instead of `norm <= bound`) takes the full step and the line prints `Error Nan`, as it does with no damping at all (which also changes line 38).

## 3. The Jacobian: lines 4, 5, 11 and 41

```
jac canary at origin, max entry error < 1e-6: true
jac canary at y0, max entry error relative to |J_ij| < 1e-6: true
jac non-symmetric at (1, 2, 3), max entry error relative to |J_ij| < 1e-6: true
jac forward difference of y^3 at y = -2 (exact 12): 11.9999998224
```

Lines 4 and 5 compare `Jac.forward (rhs 0.) y` with the canary's analytic Jacobian, `-Λ`. At the origin the entries come out exact and the absolute bound is a strict check of layout, step and divisor. At `y0 = (1, 1, 1)` cancellation limits the `1e4` entry to an absolute error of `3.0e-5`, so each error is divided by `1 + |J_ij|` ([chapter 3](03-jacobians-and-floating-point.md), sections 5 and 6, work this out). They cannot see orientation: the canary's Jacobian is diagonal, so a transposed matrix passes. Line 11 fixes that with `f(y) = (y_2, -y_1 - 3 y_2, 5 y_1 + y_3²)` at `(1, 2, 3)`, whose Jacobian `[[0, 1, 0], [-1, -3, 0], [5, 0, 6]]` is not symmetric. It also has a curved term, so a coarse step (`1e-2`) shows as truncation error there. The three lines also pass with a step of `1e-6`, so line 41 pins the step: a forward difference of `y³` at `y = -2` (exact derivative 12) is `12 - 6δ + δ²` plus round-off, and the line prints it with ten decimals, `11.9999998224` for `δ = 1e-8 (1 + |y|) = 3e-8`. A step of `1e-6 (1 + |y|)` prints `11.9999820002`, a step without the `abs` `12.0000000000` and a backward difference `12.0000001776`; dividing by the nominal step instead of the stored perturbation prints `11.9999998383`.

## 4. The canary: lines 6, 10 and 43

The canary is the decoupled linear system `y' = -Λ y`, `Λ = diag(1, 100, 1e4)`, `y(0) = (1, 1, 1)`, exact `y_i(t) = exp(-λ_i t)` ([chapter 1](01-odes-and-stiffness.md)). The three rates are far apart so that a mixed-up row or column puts a fast rate on a slow component. By `t = 1` the second component is `3.7e-44` and the third underflows to `0`, so only `y_1(1) = e^-1` matters: for the lines below the largest error is the first component's.

```
bdf1 canary t=1 dt=2e-6: max error 3.68e-07 < 1e-06: true
bdf2 canary t=1 dt=1e-3 (h lambda = 10): max error 1.53e-07 < 1e-06: true
adaptive canary t=1 tol=1e-6: max error 4.66e-07 < 1e-6: true, rhs calls 106445
```

Line 6 runs backward Euler with `Stepper.fixed`. The `3.68e-07` is predicted, not just observed: backward Euler gives `y_n = (1 + h)^(-n)` on the slow component, `(1 + h)^(-1/h) = e^-1 (1 + h/2 + ...)`, so the error at `t = 1` is about `0.184 h`, which is `3.68e-7` for `h = 2e-6` (chapter 1, section 6). The step is chosen by accuracy (`0.184 h < 1e-6` needs `h` below about `5.4e-6`) and costs `500,000` steps, which dominate the test time. It is also a hundred times below explicit Euler's stability limit `2 / 1e4 = 2e-4` (`h λ = 0.02` for the fast rate), so the line tests accuracy and the Jacobian path, not stiffness.

Line 10 is where stiffness bites. `dt = 1e-3` puts `h λ = 10` on the fast rate, five times beyond the explicit limit, and runs BDF2 with equal steps. Newton must converge on every step, and it does only with a Jacobian that is about right: the Jacobian sets how fast Newton converges, not where. For a stage equation whose multiplier of `f` is `gamma`, swapping rows 0 and 2 of the Jacobian multiplies Newton's error per iteration by about `gamma` times the fast rate. At line 6 that is `2e-6 · 1e4 = 0.02`: Newton still converges, slowly, to the same root. At line 10, `gamma` is of the order of `h = 1e-3`, the factor is several times 1, and Newton diverges. Hence swapping rows 0 and 2 in `Jac.forward` (section 9) leaves line 6 passing, with different digits (`3.68e-07` becomes `3.83e-07`, because the slower iteration stops at slightly different points), and turns line 10 into `Error Diverged`.

Line 43 runs the canary through `Adaptive.integrate` with BDF2 and `Halving` at `tol = 1e-6`, on the call budget of `Guard`, and prints the right-hand-side calls that `Instrument.count` saw next to the error. The step grows from its default start, `1e-6`, to the size the tolerance allows, so Newton meets every `gamma` on the way, and the count adds up all its iterations: each costs a Jacobian, four calls for three components, and one or more evaluations of the residual. A swap of two rows does not move the root, so the error stays under the bound; it slows Newton down or breaks it, and the count moves: `186690` with rows 0 and 1 swapped, `634048` with rows 0 and 2, `551696` with rows 1 and 2, against `106445`. That is what the line is for, because lines 6 and 10 do not fail for every swap. Rows 0 and 1 hold the rates `1` and `100`, so the factor of the paragraph above is about `gamma` times `100`, of the order of `0.1` at line 10: Newton still converges at line 10, which stays `true` with other digits. Rows 1 and 2 hold `100` and `1e4`, the factor is several times 1 and line 10 prints `Error Diverged`.

## 5. The logistic equation: lines 7, 21 and 22

```
bdf2 logistic [0,5]: error dt=0.01 3.116e-06, dt=0.005 7.810e-07, ratio 3.99 in [3.5, 4.5]: true
```

`y' = y (1 - y)`, `y(0) = 0.1`, exact `y(t) = 1 / (1 + 9 e^-t)`: an S-curve from `0.1` towards `1`, smooth and not stiff, so truncation error dominates everything else, and the right-hand side is nonlinear, so Newton genuinely iterates. Line 7 measures the order of BDF2: the errors at `t = 5` for `dt = 0.01` and `0.005`, and their ratio. A second-order method divides its error by `2² = 4` when `h` is halved, a first-order one by 2. The window `[3.5, 4.5]` accepts orders between about 1.8 and 2.2 and rejects 1 and 3. The first step is backward Euler, whose local error `O(h²)` adds only `O(h²)` at the end, so the ratio stays near 4 ([chapter 4](04-bdf.md)). `Stepper.fixed` takes equal steps (up to rounding: each is the difference of two grid times), so after the first step the step ratio `ω` is 1 up to rounding: the variable-step coefficients are exercised only by the adaptive lines.

```
adaptive logistic [0,5] dt0=0.5 tol=1e-6: accepted 1109, rejected 374, max error 8.96e-08
adaptive logistic [0,5] dt_max=1e-3 tol=1e-6: accepted 5021, rejected 0
```

Lines 21 and 22 pin the step-control rules by counting. From `dt0 = 0.5` (also the default `dt_max` here) the startup estimate (the first step has no history, so its error estimate is half the gap between backward Euler and explicit Euler) rejects the first attempts; after that the halve and double rules decide every step, so the counts and the 3-digit error move if any rule or any BDF2 coefficient at `ω ≠ 1` does ([chapter 5](05-step-control.md)). With `dt_max = 1e-3` accuracy alone would allow longer steps, so the cap sets every step once `dt` has grown to it, and the count can be derived. The default `dt0` is `1e-6 · 5 = 5e-6`; it doubles after every third accept (`3 · 5e-6 · (2^8 - 1) = 3.825e-3` covered by 24 steps of growing size) until the next doubling would pass `1e-3`; the remaining `4.996175` takes 4996 steps of `1e-3` and one shortened step: `24 + 4996 + 1 = 5021`.

## 6. Van der Pol: line 8

```
vdp mu=1000 [0,2000] tol=1e-4: Ok at t=2000, rejected >= 1: true, y finite: true
```

`y_1' = y_2`, `y_2' = μ (1 - y_1²) y_2 - y_1`, `μ = 1000`, `y(0) = (2, 0)`: a relaxation oscillation. `y_1` drifts slowly while `y_2` is small, then jumps on a time scale of about `1/μ`. In the leading-order theory for large `μ`, which keeps only the dominant terms, the first jump comes near `(3/2 - ln 2) μ ≈ 807` and the second, as long again after it, at the period `(3 - 2 ln 2) μ ≈ 1614`; during a jump `|y_2|` reaches about 1333. At the start the Jacobian is `[[0, 1], [-1, -3000]]`, so the problem is stiff from the first step ([chapter 1](01-odes-and-stiffness.md)). There is no closed form.

The case forces the step-size controller to work: long steps while `y_1` drifts, a drastic cut at each jump. It runs `Adaptive.integrate` on `[0, 2000]`, which covers both jumps, and records three facts: the run finished (`Ok` at `t = 2000`), at least one step was rejected (the rejection path really ran) and the final state is finite. It does not pin the number of steps, and nothing compares the answer with a reference. Like lines 9 and 23 it runs on the call budget of `Guard`, so a bug that makes it crawl prints `no answer within 5e6 rhs calls` instead of stalling the run. To see the numbers instead (this one needs the `refs.ml` link too), a probe prints the step counts of both adaptive runs and the errors that lines 9 and 23 bound:

```ocaml
open Vstiff

let () =
  let bdf2_halving = Adaptive.integrate (module Bdf2) (module Halving) in
  let counts name (s : Halving.stats Adaptive.solution) =
    Printf.printf "%s: accepted %d, rejected %d\n" name s.stats.accepted_steps s.stats.rejected_steps
  in
  (match bdf2_halving ~tol:1e-4 Problems.VanDerPol.problem with
   | Ok s -> counts "van der Pol" s
   | Error e -> print_endline (Fail.to_string e));
  match bdf2_halving ~tol:1e-6 Problems.Robertson.problem with
  | Ok s ->
      counts "Robertson" s;
      Printf.printf "|y1 - ref| = %.2e   |y1 + y2 + y3 - 1| = %.2e\n"
        (Float.abs (s.y.(0) -. Refs.robertson_y1_at_1e4))
        (Float.abs (Array.fold_left ( +. ) 0. s.y -. 1.))
  | Error e -> print_endline (Fail.to_string e)
```

## 7. Robertson: lines 9 and 23

```
robertson t=1e4 tol=1e-6: |y1 - ref| < 1e-3: true, |y1 + y2 + y3 - 1| < 1e-8: true
robertson t=1e4 tol=1e-6 accuracy: |y1 - ref| = 7.2e-07 < 1e-5: true
```

Robertson's chemical kinetics has three species and rate constants from `0.04` to `3e7` (a ratio of `7.5e8`). With `a = 0.04 y_1`, `b = 1e4 y_2 y_3` and `c = 3e7 y_2²` (the names in `test/problems.ml`):

```
y1' = b - a,    y2' = a - b - c,    y3' = c,    y(0) = (1, 0, 0)
```

`y_2` shoots up to a tiny maximum (about `3.6e-5`, near `t = 0.0045`) and stays tiny, `y_1` decays and `y_3` grows. Once `y_3` is near 1 the Jacobian entry `-(1e4 y_3 + 6e7 y_2)` is of order `1e4`, so a method limited by the fast rate would need steps below `2/1e4` all the way to `t = 1e4`: the problem is stiff in the sense of chapter 1.

**The reference.** There is no closed form, so the value to compare with comes from outside vstiff: `Refs.robertson_y1_at_1e4 = 0.10730042854` in [`test/refs.ml`](../../test/refs.ml), whose comment records the provenance: SciPy 1.18's `solve_ivp` with three different stiff methods, Radau (`rtol` 1e-13), BDF and LSODA (`rtol` 1e-12), all with `atol` 1e-20 and finite-difference Jacobians; the three agree to `4e-12`. Agreement of three different methods is what makes the number trustworthy. It is never edited, because a reference is only worth something while it does not depend on the code under test: changing it to make a failing run pass would turn a check against the truth into a check that vstiff agrees with vstiff. For the same reason a reference is never computed with vstiff.

**Two bounds.** Line 9 asks for `|y1 - ref| < 1e-3`, about one percent of `y_1`: loose, so a change that made the answer a hundred times worse would still pass. Line 23 asks for `1e-5` against a measured error of `7.2e-07`, which catches that (section 9). `tol` limits each step's estimated error, not the final error, which is why the two are not simply related. Only `y_1` is compared with a reference.

**The sum.** The right-hand sides add up to zero, so `y1 + y2 + y3` is constant. Summing the three components of the stage equation `x = psi + gamma f(t_{n+1}, x)` gives `Σx = Σpsi + gamma Σf = Σpsi`, and `Σpsi = 1` because `psi` combines old states with weights that sum to one ([chapter 4](04-bdf.md)). Newton keeps this property: the columns of `J` sum to zero too (differentiate `Σf = 0`), so the columns of `I - gamma J` sum to one and a full Newton step lands on `Σpsi` (exactly for the exact Jacobian; the finite-difference one is accurate to about eight digits, and Newton iterates until its step is tiny). The sum therefore stays at 1 up to round-off, and the bound `1e-8` leaves room for round-off to accumulate over the steps. It tests a property of the method that no accuracy bound would find, though it cannot see an error that only moves amount from one species to another.

## 8. Termination, the clock, the fixed-step grid and argument checks: lines 14 to 20, 25 to 37, 42 and 44

The right-hand sides of lines 14 to 20, 25 to 28 and 32 to 37 run under `Guard.budget`, a limit of 5,000,000 calls: a loop that never ends gives the text `no answer within 5e6 rhs calls` instead of hanging the run, and `Guard.run` turns an `Invalid_argument` into text too; `Report` prints it.

```
adaptive blow-up y' = y^2 from y(0)=1 to t=2: Error StepRejected 1
adaptive rhs NaN past t=0.5: Error StepRejected 46
adaptive dt0 = dt_max = 1e30 on [0, 1]: Ok at t=1
adaptive empty span: Ok at t=1
adaptive t_end < t0: Invalid_argument Adaptive.integrate: t_end is before t0
adaptive dt0 = 0: Invalid_argument Adaptive.integrate: dt0 and dt_max must be positive
bdf1 dt = 0: Invalid_argument Stepper.fixed: need dt > 0 and t_end >= t0
```

The first two lines are runs with no way forward. `y' = y²` from `y(0) = 1` has the exact solution `1 / (1 - t)` and blows up at `t = 1`; in the second, the right-hand side is `nan` past `t = 0.5`, so every step that reaches beyond it fails in Newton. In both the controller halves the failed step until it is below `16 eps |t|`, where it gives up with `StepRejected n`, `n` being the rejections in a row (in the first, the very first rejection of the final streak already meets the floor; [chapter 5](05-step-control.md) derives both counts). Without that floor both runs would still end, but late: the driver rejects every halved step below `16 eps |t|` as `Too_small` and `Halving` halves it again, until a step of 0 or `max_rejects` ends the run, and the lines print `Error StepRejected 7` (the blow-up) and `Error StepRejected 51` (the wall) (section 9). The third line pins which step is halved: the first attempt is cut to the span `1`, and a rejection must halve that cut step, not the requested `1e30`, or the same failure repeats until `max_rejects` runs out (`Error StepRejected 51`). Line 17 says that a zero span is valid. The last three lines are the messages of `Check`: a span that runs backwards (line 18), a first step of zero (line 19) and a zero step for `Stepper.fixed` (line 20). Without the checks they print `Ok at t=1`, `Error StepRejected 1` and `Ok [0.500000000000]` instead (section 9).

```
adaptive y' = 1 over one ulp from t=1: Ok, at t_end: true, y = 2.22045e-16
adaptive y' = 1 from t=1e10, dt_max = 5e-7 is below t's resolution: Error StepRejected 1
adaptive y' = 1 from t=1e15 over 100, dt0 = dt_max = 3.7: Ok, at t_end: true, y = 100
```

Lines 25 to 27 pin the rules that keep the state and the clock together ([chapter 5](05-step-control.md), section 9). The problem is `y' = 1` from `y = 0`, so `y` is the time elapsed and a correct run ends with `y = t_end - t0`; each line prints whether the run landed on `t_end`, and `y`. The clock is a float, so near `t` it moves in whole ulps: a step below half an ulp of `t` does not move it, and a step of a few ulps moves it by a different amount than its length.

Line 26 starts at `t = 1e10`, where the floats are `1.9e-6` apart ([chapter 3](03-jacobians-and-floating-point.md)) and the resolution is `3.6e-5`, with steps of at most `5e-7`. Each snaps to `h = (t + dt) - t = 0`, is rejected as `Too_small` without calling the method, and `Halving` halves 0, which is below the floor, so the first rejection ends the run with `StepRejected 1`. Unsnapped, a step of `5e-7` is still below the resolution and is rejected the same way; only with snapping and the rejection of short steps both gone are these steps accepted without moving `t`, so that the run never finishes and the budget text takes the line's place.

Line 27 starts at `t = 1e15`, where the floats are `0.125` apart and the resolution is `3.55`, with steps of `3.7`, just above it, so `t + 3.7` is `t + 3.75`. An unsnapped step advances the state by `3.7` and the clock by `3.75`: after 26 steps the state is `26 · 3.7 = 96.2` where the clock has moved `97.5`, and the last step, the remaining `2.5`, ends the state at `98.7`, an `Ok` with a wrong answer. Snapped, every step is `3.75`, the last one `2.5` again, and the state ends at `100`. Steps below the resolution, `dt0 = dt_max = 0.19` for instance, are rejected before they can be snapped: that run ends with `Error StepRejected 1` ([chapter 5](05-step-control.md), section 9). Line 25 pins the remainder rule. The span is one ulp of 1, `2.2e-16`, and the default first step, a millionth of the span, snaps to 0. A remainder of at most `Clock.resolution t`, that is `16 eps |t|` (`3.6e-15` at `t = 1`), is taken whole as the last step, so the run ends `Ok` on `t_end` with `y` equal to that one ulp; without the rule it ends with `StepRejected 1`. Section 9 lists what the three lines print when the rules are removed.

Line 28 pins the `Too_small` rejection itself, which the outcome of line 26 cannot tell apart: without it the method is called with `h = 0`, the second such step makes `ω = 0 / 0`, BDF2 fails as `Solver Nan` and the run ends with the same `StepRejected 1`. So line 28 runs the problem of line 26 with its right-hand side wrapped by `Instrument.count` and prints the number of calls. With the rejection there is one, the driver's check of the initial state, because the step that cannot move `t` never reaches the method.

```
bdf1 y' = sqrt(1 - t) on [0, 1], dt = 1/9: Ok [0.603925945409]
bdf2 y' = 2t on [0, 1]: error dt=0.1 1.500e-02, dt=0.05 3.750e-03, ratio 4.00 in [3.5, 4.5]: true
bdf1 dt = 1e-7 at t=1e10, below t's resolution: Invalid_argument Stepper.fixed: dt is below the resolution of t
adaptive rhs NaN past t=0, max_rejects = 5000: Error StepRejected 1055
adaptive tol = 0: Invalid_argument Adaptive.integrate: tol must be positive
adaptive rhs of the wrong length: Invalid_argument Adaptive.integrate: rhs returns a vector of the wrong length
bdf1 rhs of the wrong length: Invalid_argument Stepper.fixed: rhs returns a vector of the wrong length
bdf1 y0 = nan on an empty span: Error Nan
adaptive over a span of 1e-320: Ok at t=9.99989e-321
adaptive y' = 2t on [0, 1] tol=1e-6: max error 1.57e-12
adaptive y0 = nan: Error Nan
```

Lines 29 and 31 pin the grid of `Stepper.fixed` ([chapter 5](05-step-control.md), section 9): step `k` ends at `t0 + k h` with `h = span / n`, the last step at `t_end` itself, and the method is given the difference of the two end times. Line 29 is `y' = sqrt(1 - t)` on `[0, 1]` with `dt = 1/9`. The right-hand side does not depend on `y`, so backward Euler's state is `h` times the sum of `sqrt(1 - t_k)` over the nine end times `t_k = k/9`, the last term being `sqrt 0`: `(sqrt 0 + sqrt 1 + ... + sqrt 8) / 27 = 0.603926`, the printed value (the integral is `2/3`). A running sum of nine steps of `1/9` ends at `1.0000000000000002`, `sqrt` of a negative number is `nan`, and the line printed `Error Nan`. Line 31 starts at `t = 1e10`, where the floats are `1.9e-6` apart, with `dt = 1e-7`: grid times that close repeat or jump by a whole `1.9e-6`, so the method would get steps of 0 and `1.9e-6` and return the error of those long steps, about `0.184 · 1.9e-6 = 3.5e-7` above `e^-1` (section 4), without a word. `Check.fixed` refuses a `dt` below `Clock.resolution` of the larger of `|t0|` and `|t_end|`, which also keeps the number of steps at most `2^49`.

Lines 30 and 42 pin the time argument of the right-hand side. `y' = 2t` with `y(0) = 0` has the solution `t²`, a quadratic, which BDF2 reproduces from exact data, so its error is what the start leaves: the first step is backward Euler, whose local error here is `h²` in size, and BDF2 adds none. The printed errors are `1.5 h²` for `h = 0.1` and `0.05` to the digits shown, and their ratio is 4 to the three digits printed. Evaluate the stage at the old time `t_n` instead of `t_{n+1}` and BDF2 is no longer exact on quadratics: the ratio falls to 1.75 (to 1.92 if backward Euler's stage is evaluated at the old time too). Line 42 is the adaptive driver on the same problem. The default first step is `1e-6`, with a start error of `(1e-6)²`, and the run ends with `1.57e-12`, about `1.5 (1e-6)²`; with the stage at the old time it ends with `1.02e-03`.

Line 32 is the NaN wall of line 15 moved to `t = 0`: the right-hand side is `nan` for every `t > 0`, `max_rejects` is 5000 and the default first step `1e-6` fails and is halved. At `t = 0` the floor `16 eps |t|` is 0, so only the limit and the underflow of the step to 0 can end the run. The smallest positive float is `2^-1074`, about `5e-324`; halving `1e-6` reaches it at the 1054th halving and gives exactly 0 at the 1055th, where `Halving` gives up: `Error StepRejected 1055`. Without that test the 0 was rejected as `Too_small` again and again until the limit ran out, `Error StepRejected 5001`, and for a huge limit it never ended ([chapter 5](05-step-control.md), section 9).

Lines 33 to 37 pin what the drivers check before they step, and line 44 adds one. Line 33: `tol = 0` raises (without the check the run gave up with `Error StepRejected 4`). Lines 34 and 35: each driver evaluates `rhs t0 y0` once, and `Check.output` raises, naming the driver, when the vector has another length than `y0` (here a right-hand side that returns one component for a state of two; it printed `Invalid_argument index out of bounds`). Line 36: `Stepper.fixed` returns `Error Nan` for a state that is not finite, as `Adaptive.integrate` does, even when the span is empty (it printed `Ok [nan]`). Line 37: on a span of `1e-320`, a denormal float, the default first step, a millionth of the span, underflows to 0, and a run would be refused as `dt0 = 0` (`Invalid_argument`); the default now falls back to the span, and the run ends `Ok` at the float nearest `1e-320`, `9.99989e-321`. Line 44 is the check of `Adaptive.integrate`: `y0 = nan` is `Error Nan` (without it every attempt fails, and the rejections run up to the limit: `Error StepRejected 51`).

## 9. The soak, and what each pin catches

```
soak canary x10: passed 10/10, identical: true
soak logistic order x10: passed 10/10, identical: true
soak van der Pol x10: passed 10/10, identical: true
soak robertson x10: passed 10/10, identical: true
```

[`test/soak.ml`](../../test/soak.ml) reruns the canary, the logistic order check, van der Pol and Robertson ten times each, with the same criteria as lines 6, 7, 8 and 9, and checks that every round passes and that all rounds are equal under `=`. The library has no hidden state, so `identical` can turn `false` only if someone adds hidden state, randomness or parallelism (or a NaN appears, since NaN is never equal to itself). The van der Pol and Robertson rounds run on the call budget of `Guard`: `Guard.bounded` turns a round whose budget ran out into `None`, which counts as not passed, so a bug that makes them crawl ends a round instead of stalling the soak (the canary and the logistic order check take fixed steps and always end). It is a tripwire, not extra coverage. The mechanics are in [docs/testing.md](../testing.md).

Each row is one deliberate change, an edit of a few lines at most, that you can repeat in a scratch copy of the repository, running the built test programs under a time limit; [docs/testing.md](../testing.md) has the method and a longer table. "Silent" means no line of either expected file changes.

| Change | What the corpus does |
|--------|----------------------|
| `Jac.forward`: swap rows 0 and 2 of a three-row result | both `jac canary` lines and line 11 print `false`; line 10 prints `Error Diverged`; line 6 still passes but `3.68e-07` becomes `3.83e-07`; line 43 prints `rhs calls 634048` and stays `true` |
| `Jac.forward`: swap rows 0 and 1, or rows 1 and 2, of a three-row result | the three Jacobian lines print `false`; line 43 prints `rhs calls 186690`, respectively `551696`, and stays `true`; line 10 still passes for rows 0 and 1, with other digits, and prints `Error Diverged` for rows 1 and 2; the Robertson lines (9 and 23) print the budget text, and the soak's Robertson line `passed 0/10` |
| `Jac.forward`: return the transpose | line 11 prints `false`, and the Robertson lines (9 and 23) print the budget text, in about 2 s; van der Pol (line 8) still passes; the soak ends in about 14 s with `soak robertson x10: passed 0/10, identical: true` instead of stalling |
| `Jac` step `1e-2` or `1e-14` instead of `1e-8` | line 11 prints `false` for both; line 5 also for `1e-14`; line 41 moves for both; lines 9 and 23 print the budget text for `1e-2` (and the soak's Robertson line `passed 0/10`); line 43 moves for `1e-14`. A step `1e-6`, no `abs`, a backward difference or the nominal step as divisor moves only line 41 |
| `Linalg`: pivot always on row 0 | line 12 prints `None`, line 13 prints `[0.000000000000; 1.000000000000]` |
| Newton tolerance `1e-6` instead of `1e-10` | the 12th decimal of lines 2, 45 and 47; line 46 prints an `Ok` instead of `Error Diverged`; the count of line 43 moves |
| Newton: no stopping test on the step | line 2 and every line that takes an integration step (6 to 10, 14 to 16, 21 to 23, 29, 30, 37, 42 and 43), and line 45, print `Error Diverged` or `Error StepRejected n`; line 24 prints `Ok` with `max_float` in full (309 digits), line 27 `Error StepRejected 1`, and lines 25, 26 and 28 do not change |
| Newton: the step test applied only after the line search accepts a step | lines 6, 7, 10, 29 and 30 print `Error Diverged`; lines 8, 16, 21, 22 and 42 print `Error StepRejected`, lines 14 and 15 other counts; line 24 prints `max_float` in full and line 27 `Error StepRejected 1`; line 43 prints the budget text; the soak lines for the canary, the logistic order and van der Pol pass 0/10 |
| Newton: no finiteness check on the converged `x + dx` | line 24 prints `Ok [inf]` |
| Newton: no damping, or a line search that accepts any finite trial point | line 38 prints `Error Diverged`; without damping at all line 47 prints `Error Nan` too |
| Newton: a `nan` residual at the start named `Diverged`, or a non-finite step named `Nan` | line 39, respectively line 40, prints the other name |
| Newton: smallest damping factor `1/2`; iteration limit 100; a line search that lets a `nan` trial point pass | line 45 prints `Error Diverged`; line 46 prints `Ok [0.000000000138]`; line 47 prints `Error Nan` |
| `Halving`: double after two accepts, not three | lines 21, 22, 42 and 43, and the counts of the blow-up (line 14, `StepRejected 2`) and the NaN wall (line 15) |
| `Bdf2`: startup estimate factor `1e-3` instead of `0.5`; `Halving`: no cap by `dt_max` | lines 21 and 43; line 22 |
| `Halving`: error weight 1 instead of `1 + abs(y_i)` | the blow-up and NaN-wall lines (14 and 15), lines 21, 23, 42 and 43 |
| `Halving`: halve the requested step, not the failed one | the `1e30` line (16): `Error StepRejected 51` |
| `Halving`: no `16 eps abs(t)` floor | lines 14 and 15 print `Error StepRejected 7` and `Error StepRejected 51` instead of 1 and 46: the driver rejects the halved steps below the floor as `Too_small`, and they are halved down to a step of 0 or to `max_rejects`; lines 26 and 28 still end at once, because the halved step is 0 |
| `Halving`: no give-up on a halved step of 0 | line 32 prints `Error StepRejected 5001` instead of 1055 |
| `Adaptive`: do not shorten the last step (`h = dt` where it takes `t_end - t`) | the `1e30` line (16) and lines 23, 42 and 43; on the clock lines the state advances by `dt` instead of the remainder, `y = 2.22045e-22` on line 25 and `y = 101.2` on line 27 (26 steps of `3.75` and a last one of `3.7`) |
| `Adaptive`: no snapping (`h = dt` for a step that is not the last) | line 27 prints `y = 98.7`, and lines 42 and 43 change; lines 26 and 28 still end at once, because their steps are below the resolution and rejected anyway |
| `Adaptive`: no remainder rule (`last = dt >= remaining`) | line 25 prints `Error StepRejected 1` |
| `Adaptive`: no snapping, and a step below the resolution rejected only when `h <= 0` | lines 26 and 28 print the budget line, line 27 prints `y = 98.7`, and lines 42 and 43 change |
| `Adaptive`: neither snapping, nor the remainder rule, nor any `Too_small` rejection | lines 25, 26 and 28 print the budget line, line 27 prints `y = 98.7`, and lines 42 and 43 change |
| `Adaptive`: no `Too_small` rejection (call the method with `h = 0`) | line 28 prints `rhs calls: 5` |
| `Adaptive`: a non-final step rejected only when `h <= 0`, not below the resolution | silent |
| `Adaptive`: no fallback to the span when a default step underflows to 0 | line 37 prints `Invalid_argument Adaptive.integrate: dt0 and dt_max must be positive` |
| `Adaptive`: no check of a non-finite `y0` or `rhs t0 y0` | line 44 prints `Error StepRejected 51` instead of `Error Nan` |
| `Check`: remove the argument checks | lines 18 to 20 print `Ok at t=1`, `Error StepRejected 1` and `Ok [0.500000000000]`; lines 31 and 33 to 35 print `Ok [0.367879792008]`, `Error StepRejected 4` and `Invalid_argument index out of bounds` twice |
| `Stepper.fixed`: a running sum for `t` instead of the grid; no finiteness check on `y0` and `rhs t0 y0` | line 29 prints `Error Nan`; line 36 prints `Ok [nan]` |
| `Bdf2`: return backward Euler's result, not BDF2's | line 23 prints `1.6e-04 < 1e-5: false`, lines 21 and 42 change, line 43 prints `false`, and the blow-up (line 14) prints `StepRejected 2`, while line 9 still passes |
| `Bdf2.coeffs`: `beta` off by 1 part in `1e5` times `(ω - 1)` | the digits of lines 21, 42 and 43 and the count of the blow-up line |
| `Bdf2` or `Bdf1`: stage at `t_n` instead of `t_{n+1}` | `Bdf2`: line 30 prints the ratio `1.75` and line 42 `1.02e-03`; `Bdf1`: line 29 prints `Ok [0.715037056520]`, line 32 `StepRejected 48` and line 42 `1.60e-12` |
| `Halving`: a rejection does not reset the accept streak | line 21: the `rejected` count changes |
| Newton: iteration limit 1 | lines 8, 9, 14 to 16, 23 and 43 print the budget text, but lines 21, 22 and 42 have none and run on: the corpus does not finish in practice; the soak ends, with every line at `passed 0/10` |

## 10. What the corpus does not catch

These changes leave every line of both expected files unchanged. Each suggested case gives a different result on the broken version (the Armijo constant is the one exception, as noted).

- **Newton's safeguards.** Lines 38 and 45 to 47 pin damping, its depth, the limit of 50 iterations and the refusal of a `nan` trial point. What is left is the Armijo constant: `0` instead of `1e-4` changes no line, and none of the cases of [chapter 2](02-newton.md) tells the two apart. The finiteness tests on the step and on the trial residual are redundant in behaviour (`nan <= bound` is already false, [chapter 3](03-jacobians-and-floating-point.md), and the line search refuses the trial point of a non-finite step), and the one on `x` at the start never fires alone in the corpus, so no line pins them; line 47 pins that the trial point is refused, not which test refuses it. The test on `G(x)` has line 39 and the check on a converged `x + dx` line 24. The order of the step test and the line search is not a gap: reversing it moves many lines (section 9).
- **Cost.** Starting Newton from the previous state instead of the extrapolated guess changes the number of right-hand-side calls but no line. Line 43 counts calls, but on the canary, which is linear: Newton lands on the root from any guess. A case that counts the calls of a nonlinear run with `Instrument.count`, the logistic equation for instance, would tell the two apart.
- **Step-control defaults and limits.** The default `dt_max` (`span / 10`), the `max_rejects` give-up (the NaN wall ends by the floor, at 46 rejections, below the default 50, and line 32 by the underflow to 0, at 1055, below its limit of 5000), and accepting a `nan` error estimate are all silent. Cases: `y' = 0` on `[0, 1]` accepts every step, so the count follows from the rules alone, like line 22 (`dt0 = 1e-6` doubles every third accept up to `0.1`: `3 · 17 = 51` steps cover `0.39`, and 7 more finish, 58 in all; another `dt_max` or doubling rule gives another count); `~max_rejects:3` on the NaN wall gives `Error StepRejected 4`; `Halving.acceptable` with `err = [| nan |]` is `false`.
- **Short steps.** The rejection of a non-final step that is below `Clock.resolution t` but still moves `t` is silent: a driver that rejects only `h <= 0` leaves both files unchanged, because the steps of lines 26 and 28 snap to 0 and line 27 runs above the resolution. The rule keeps BDF2's step ratio at most `66/31` ([chapter 5](05-step-control.md), section 9), so it deserves a case. Case: `y' = 1` from `t = 1e15` over 100 with `dt0 = dt_max = 0.19`, where the resolution is `3.55`, ends `Error StepRejected 1` without a call to the method; with only `h <= 0` rejected it ends `Ok` with `y = 100`, the steps snapping to `0.25`. A probe that records `ω` (the wrapper of [chapter 4](04-bdf.md), section 7) on the same problem with `dt0 = 0.18` and a large `dt_max` prints the ratio `3` when only `h <= 0` is rejected, and nothing to record otherwise, because the run ends at once.
- **The fixed-step driver.** Three rules of `Stepper.fixed` are silent: at least one step, steps of `span / n` rather than of `dt`, and a last grid time that is `t_end` itself with each step the difference of its end times. In every corpus case `dt` is a whole fraction of the span and the grid reaches `t_end` by itself. Cases: `dt = 4` on the canary's `[0, 1]` must take one step of size 1 (without that rule the state comes back unchanged); `dt = 0.3` must take three steps of about `1/3`, not steps of `0.3`, `0.3` and `0.4`; `y' = sqrt(0.9 - t)` on `[0, 0.9]` in seven steps (`dt = 0.9 / 7`) is `Ok`, and `Error Nan` when the last rule is dropped, because the last stage then lies past `t_end` by an ulp (`7 h` is `0.9000000000000001`). The program after this list runs the last one.
- **Vectors.** `Vec` has no case of its own. `Vec.norm_inf [| nan; 1. |]` should be `nan`. A version built on the generic `max`, or on `if Float.abs v > m then Float.abs v else m`, returns `1.` (with the NaN last, `[| 1.; nan |]`, the generic `max` would still pass).
- **What the lines cannot see.** Van der Pol is not compared with a reference; only `y_1` of Robertson is; the soak is a tripwire, not coverage; `Instrument` counts right-hand-side calls only, not Newton iterations. The `500,000`-step canary dominates the run time.

The case for the last rule of the fixed-step driver, as a program (it prints `Ok [0.499365140919]`):

```ocaml
open Vstiff

let () =
  let problem = { Ode.rhs = (fun t _ -> [| Float.sqrt (0.9 -. t) |]); t0 = 0.; t_end = 0.9; y0 = [| 0. |] } in
  match Stepper.fixed (module Bdf1) ~dt:(0.9 /. 7.) problem with
  | Ok y -> Printf.printf "Ok [%.12f]\n" y.(0)
  | Error e -> print_endline ("Error " ^ Fail.to_string e)
```

## In the code

| Idea | Where |
|------|-------|
| The four problems: `rhs`, `y0`, `problem`, and `exact` for the first two | `Canary`, `Logistic`, `VanDerPol`, `Robertson` in [`test/problems.ml`](../../test/problems.ml) |
| The external Robertson reference | `Refs.robertson_y1_at_1e4` in [`test/refs.ml`](../../test/refs.ml), provenance in its doc comment |
| The cases | the named lists of [`test/corpus.ml`](../../test/corpus.ml), mapped to line numbers in section 1; `corpus` concatenates them in output order |
| Test effects | `Guard.budget`, `Guard.run` and `Guard.bounded` in [`test/guard.ml`](../../test/guard.ml); `Report.lines` in [`test/report.ml`](../../test/report.ml) |
| The expected output | [`test/corpus.expected`](../../test/corpus.expected) (47 lines), [`test/soak.expected`](../../test/soak.expected) (4 lines) |
| The soak | [`test/soak.ml`](../../test/soak.ml): the module type `Case`, `repeat`, `soak`, `on_budget` and four case modules |
| How both run | [`test/dune`](../../test/dune): `(tests (names corpus soak) (libraries vstiff numerics))`; `corpus.ml` and `soak.ml` start with `open Vstiff` and `open Numerics` |

## Check yourself

1. **With Jacobian rows 0 and 2 swapped, which corpus lines change, and which changed lines still say `true`?**
   Both `jac canary` lines and line 11 turn `false`, and line 10 becomes `Error Diverged`. Line 6 changes from `3.68e-07` to `3.83e-07` and stays `true`: the roots of the stage equation do not depend on the Jacobian, and the slower iteration only stops at slightly different points. Line 43 stays `true` too and prints `rhs calls 634048` instead of `106445`.

2. **A row swap in `Jac.forward` leaves line 6 passing but breaks line 10. Why?**
   The swapped matrix only slows Newton down or breaks it, depending on the step. The error is multiplied per iteration by about `gamma` times the fast rate: `2e-6 · 1e4 = 0.02` on line 6, so Newton still converges to the same root, and several times 1 on line 10, so it diverges.

3. **Why must the Robertson reference never be computed by vstiff?**
   A reference produced by the code under test shares that code's mistakes, so it cannot detect them. It comes from SciPy, from three different methods that agree to `4e-12`.

4. **What ratio would a first-order method print on line 7, and why does the window start at 3.5?**
   About 2: halving `h` halves the error. The window `[3.5, 4.5]` is centred on the second-order value 4 and leaves room for higher-order terms at `dt = 0.01`.

5. **Why is the conservation bound `1e-8` and not `1e-15`?**
   The sum is preserved up to round-off at each step, but round-off accumulates over many steps, so a bound at the level of one rounding error would be fragile. `1e-8` still fails a conservation bug of any real size.

6. **Derive the `5021` of line 22.**
   `dt0 = 5e-6` and `dt_max = 1e-3`: 8 groups of 3 steps (`5e-6 · 2^j`) cover `3.825e-3`, then 4996 steps of `1e-3` and one shortened step: `24 + 4996 + 1 = 5021`.

7. **Name one bug the corpus cannot see, and a case that would.**
   `Halving` accepting a step whose estimate is `nan` (`not (est > tol)` instead of `est <= tol`): `Halving.acceptable` with `err = [| nan |]` is `false`, and `true` with the bug.

8. **What does line 24 pin, and what does it print without the pin?**
   That `Newton.solve` does not return an overflowed `x + dx` as a root: it prints `Error Nan`, where the unchecked sum prints `Ok [inf]`. The step test is relative, so near the largest float a step of `4e297` counts as tiny.

9. **Line 26 ends in `StepRejected 1` although no step is ever computed. Why, and which rules keep it from running on?**
   At `t = 1e10` the floats are `1.9e-6` apart and the resolution is `3.6e-5`, so a step of `5e-7` is far too short: it snaps to `h = 0`, the driver rejects it as `Too_small` without calling the method, and `Halving` halves 0, below the floor, so the first rejection ends the run. Without snapping the step is still below the resolution and is rejected the same way; without both rules the steps are accepted without moving `t` and the run never ends: the budget stops it (`no answer within 5e6 rhs calls`).

10. **Line 29 prints `Error Nan` when `t` is a running sum. Why?**
   Nine additions of `1/9` give `1.0000000000000002`, so the last stage is evaluated past `t_end = 1`, where `sqrt (1 - t)` is `nan`. The grid `t0 + k h` reaches 1 itself and its last time is `t_end`.

11. **Line 32 ends with `StepRejected 1055` although `max_rejects` is 5000. Where does 1055 come from?**
   At `t = 0` the floor is 0, so the halving of the default first step `1e-6` goes on until the step is 0: 1054 halvings reach the smallest float and the 1055th gives 0, where `Halving` gives up.

12. **Line 43 prints a count of right-hand-side calls. What does that catch that lines 6 and 10 do not?**
   A row swap leaves the root alone, so the error stays under the bound; it changes how many iterations Newton needs along the run, so the count moves for every pair of rows (`186690`, `634048` and `551696` for rows 0 and 1, 0 and 2, 1 and 2, against `106445`), while line 10 passes for rows 0 and 1 and line 6 for every pair.

13. **Which of Newton's limits do lines 45 to 47 pin, and which is still free?**
   The depth of the line search (`x² - 2` from `1e-3` needs the factor `1/512`), the limit of 50 iterations (`x³` from 1 needs 55) and the refusal of a `nan` trial point (`log x` from 3). The Armijo constant `1e-4` is free: `0` changes no line.

That is the end of the numerical-methods track. Next: [docs/testing.md](../testing.md) for running the tests and changing expectations, and [docs/exercises.md](../exercises.md) for exercises and starter contributions.
