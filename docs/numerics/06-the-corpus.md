# 6. The test corpus

Chapter 6 of 6 in the numerical-methods track. Previous: [5. Error estimates and step-size control](05-step-control.md). Index and reading order: [docs/README.md](../README.md). Terms and symbols: [glossary](../glossary.md).

**Summary.** The corpus is a table of 23 cases whose printed lines are pinned in [`test/corpus.expected`](../../test/corpus.expected). This chapter says, line by line, what each one proves, what a failure points at and what it cannot see: nine lines that climb from Newton to a stiff chemical system, and fourteen regression pins that each guard one rule. It introduces the four ODEs behind them (the canary, the logistic equation, van der Pol and Robertson), the external reference value and why it is never edited, the four soak lines, a table of deliberate bugs with the line that notices each, and an honest list of what the corpus does not catch, each item a possible contribution.

Running the tests, reading a failure, adding a case and the policy for changing expectations are in [docs/testing.md](../testing.md), which this chapter does not repeat. The mathematics comes from chapters [1](01-odes-and-stiffness.md) to [5](05-step-control.md). Snippets run in the probe project (set up on day 1 from Setup in [exercises.md](../exercises.md), or "Running the snippets" in chapter 1); those that use `Refs` also need a link to `test/refs.ml`, made from the repository root next to the one for `problems.ml`:

```sh
ln -sfn "$(pwd)/test/refs.ml" ../vstiff-scratch/p/refs.ml
```

OCaml notes for them ([docs/ocaml.md](../ocaml.md) has the details): `(module Bdf2)` passes a module as an argument; `{ Ode.rhs = ...; t0 = ... }` builds a record, the `Ode.` prefix saying which module its fields belong to; `(s : T)` is a type annotation that tells OCaml which record type `s` has, so that `s.stats` can be read.

## 1. What the corpus is for

[`test/corpus.ml`](../../test/corpus.ml) is a table of cases, each a name and a function that computes the text printed after it. `Report.lines` prints `name: text` for each, in order, and `dune runtest` compares the output with `test/corpus.expected`. Several cases print a boolean (`true` means the criterion held), so the criterion lives in the code and the expected file pins its outcome; others print digits, which pin more. The first nine lines, from six development steps, climb the library layer by layer:

```
Newton + Linalg  →  Jac  →  Stage + Bdf1  →  Bdf2  →  Adaptive + Halving  →  a stiff chemical system
```

The other fourteen are regression pins, each added for a specific mistake (section 9 lists what catches what). The lines are numbered in the order of `corpus.expected`, which is the order of the named lists in `corpus.ml`: `newton` (1 to 3), `jacobian` (4 and 5), `backward_euler` (6), `bdf2_order` (7), `van_der_pol` (8), `robertson` (9), then the pins `stiff_canary` (10), `orientation` (11), `pivoting` (12 and 13), `give_up` (14 to 20), `step_control` (21 and 22) and `robertson_accuracy` (23). Two rules of thumb: when several lines differ, start with the first, which involves the fewest layers; and a case earns its place by pinning something the others do not.

The problems live in [`test/problems.ml`](../../test/problems.ml). Each module exposes `rhs`, `y0`, a `problem` record and, where the exact solution has a formula, `exact`. None supplies a Jacobian: the solver must work for any black-box right-hand side. All four right-hand sides ignore their time argument (`rhs _t y`); section 10 comes back to that.

| Problem | Components | Time scales | Corpus interval | Stiff? | Closed form |
|---------|-----------|-------------|-----------------|--------|-------------|
| Canary | 3 | rates 1, 100, 1e4 | `[0, 1]` | yes, ratio `1e4` | yes |
| Logistic | 1 | about 1 | `[0, 5]` | no | yes |
| Van der Pol, `μ = 1000` | 2 | slow drifts, jumps on `1/μ = 1e-3` | `[0, 2000]` | yes | no |
| Robertson | 3 | rate constants 0.04 to 3e7 | `[0, 1e4]` | yes | no |

## 2. Newton and Linalg: lines 1 to 3, 12 and 13

```
newton linear 2d: Ok [2.000000000000; 3.000000000000]
newton quadratic x0=1: Ok [1.414213562373]
newton quadratic x0=0: Error Diverged
linalg zero leading pivot: [1.000000000000; 2.000000000000; 3.000000000000]
linalg tiny leading pivot: [1.000000000000; 1.000000000000]
```

Lines 1 to 3 test `Newton.solve` with a hand-written Jacobian, `Vec` and `Linalg`, and nothing above them. Line 1 is `G(x) = A x - b` with `A = [[3, 1], [1, 2]]`, `b = (9, 8)`, from `(0, 0)`: one Newton step ([chapter 2](02-newton.md)), printed with 12 decimals. Line 2 is `x² - 2` from 1, whose iterates are in chapter 2; 12 printed decimals make the line sensitive to Newton's tolerance. Line 3 is the same function from 0, where the 1×1 Jacobian is singular: `Linalg.solve` returns `None` and the failure comes back as the value `Error Diverged`. What they do not exercise: `A` is symmetric with its largest entry already in the corner, so neither orientation nor a row exchange is tested, and no step needs shortening.

Lines 12 and 13 close the row-exchange gap, both worked by hand in chapter 2: a zero in the corner that only a row exchange gets past (without it the answer is `None`), and a tiny corner entry that elimination without the exchange turns into a first unknown of `0` instead of `1`.

## 3. The Jacobian: lines 4, 5 and 11

```
jac canary at origin, max entry error < 1e-6: true
jac canary at y0, max entry error relative to |J_ij| < 1e-6: true
jac non-symmetric at (1, 2, 3), max entry error relative to |J_ij| < 1e-6: true
```

Lines 4 and 5 compare `Jac.forward (rhs 0.) y` with the canary's analytic Jacobian, `-Λ`. At the origin the entries come out exact and the absolute bound is a strict check of layout, step and divisor. At `y0 = (1, 1, 1)` cancellation limits the `1e4` entry to an absolute error of `3.0e-5`, so each error is divided by `1 + |J_ij|` ([chapter 3](03-jacobians-and-floating-point.md), sections 5 and 6, work this out). They cannot see orientation: the canary's Jacobian is diagonal, so a transposed matrix passes. Line 11 fixes that with `f(y) = (y_2, -y_1 - 3 y_2, 5 y_1 + y_3²)` at `(1, 2, 3)`, whose Jacobian `[[0, 1, 0], [-1, -3, 0], [5, 0, 6]]` is not symmetric. It also has a curved term, so a coarse step (`1e-2`) shows as truncation error there.

## 4. The canary: lines 6 and 10

The canary is the decoupled linear system `y' = -Λ y`, `Λ = diag(1, 100, 1e4)`, `y(0) = (1, 1, 1)`, exact `y_i(t) = exp(-λ_i t)` ([chapter 1](01-odes-and-stiffness.md)). The three rates are far apart so that a mixed-up row or column puts a fast rate on a slow component. By `t = 1` the second component is `3.7e-44` and the third underflows to `0`, so only `y_1(1) = e^-1` matters: for both lines below the largest error is the first component's.

```
bdf1 canary t=1 dt=2e-6: max error 3.68e-07 < 1e-06: true
bdf2 canary t=1 dt=1e-3 (h lambda = 10): max error 1.53e-07 < 1e-06: true
```

Line 6 runs backward Euler with `Stepper.fixed`. The `3.68e-07` is predicted, not just observed: backward Euler gives `y_n = (1 + h)^(-n)` on the slow component, `(1 + h)^(-1/h) = e^-1 (1 + h/2 + ...)`, so the error at `t = 1` is about `0.184 h`, which is `3.68e-7` for `h = 2e-6` (chapter 1, section 6). The step is chosen by accuracy (`0.184 h < 1e-6` needs `h` below about `5.4e-6`) and costs `500,000` steps, which dominate the test time. It is also a hundred times below explicit Euler's stability limit `2 / 1e4 = 2e-4` (`h λ = 0.02` for the fast rate), so the line tests accuracy and the Jacobian path, not stiffness.

Line 10 is where stiffness bites. `dt = 1e-3` puts `h λ = 10` on the fast rate, five times beyond the explicit limit, and runs BDF2 with equal steps. Newton must converge on every step, and it does only with a Jacobian that is about right: the Jacobian sets how fast Newton converges, not where. For a stage equation whose multiplier of `f` is `gamma`, swapping rows 0 and 2 of the Jacobian multiplies Newton's error per iteration by about `gamma` times the fast rate. At line 6 that is `2e-6 · 1e4 = 0.02`: Newton still converges, slowly, to the same root. At line 10, `gamma` is of the order of `h = 1e-3`, the factor is several times 1, and Newton diverges. Hence a row swap in `Jac.forward` (section 9) leaves line 6 passing, with different digits (`3.68e-07` becomes `3.83e-07`, because the slower iteration stops at slightly different points), and turns line 10 into `Error Diverged`.

## 5. The logistic equation: lines 7, 21 and 22

```
bdf2 logistic [0,5]: error dt=0.01 3.116e-06, dt=0.005 7.810e-07, ratio 3.99 in [3.5, 4.5]: true
```

`y' = y (1 - y)`, `y(0) = 0.1`, exact `y(t) = 1 / (1 + 9 e^-t)`: an S-curve from `0.1` towards `1`, smooth and not stiff, so truncation error dominates everything else, and the right-hand side is nonlinear, so Newton genuinely iterates. Line 7 measures the order of BDF2: the errors at `t = 5` for `dt = 0.01` and `0.005`, and their ratio. A second-order method divides its error by `2² = 4` when `h` is halved, a first-order one by 2. The window `[3.5, 4.5]` accepts orders between about 1.8 and 2.2 and rejects 1 and 3. The first step is backward Euler, whose local error `O(h²)` adds only `O(h²)` at the end, so the ratio stays near 4 ([chapter 4](04-bdf.md)). `Stepper.fixed` takes equal steps, so after the first step the step ratio `ω` is 1: the variable-step coefficients are exercised only by the adaptive lines.

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

The case forces the step-size controller to work: long steps while `y_1` drifts, a drastic cut at each jump. It runs `Adaptive.integrate` on `[0, 2000]`, which covers both jumps, and records three facts: the run finished (`Ok` at `t = 2000`), at least one step was rejected (the rejection path really ran) and the final state is finite. It does not pin the number of steps, and nothing compares the answer with a reference. To see the numbers instead (this one needs the `refs.ml` link too), a probe prints the step counts of both adaptive runs and the errors that lines 9 and 23 bound:

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

## 8. Termination and argument checks: lines 14 to 20

Each right-hand side here runs under `Guard.budget`, a limit of 5,000,000 calls: a loop that never ends gives the text `no answer within 5e6 rhs calls` instead of hanging the run, and `Guard.run` turns an `Invalid_argument` into text too; `Report` prints it.

```
adaptive blow-up y' = y^2 from y(0)=1 to t=2: Error StepRejected 1
adaptive rhs NaN past t=0.5: Error StepRejected 46
adaptive dt0 = dt_max = 1e30 on [0, 1]: Ok at t=1
adaptive empty span: Ok at t=1
adaptive t_end < t0: Invalid_argument Adaptive.integrate: t_end is before t0
adaptive dt0 = 0: Invalid_argument Adaptive.integrate: dt0 and dt_max must be positive
bdf1 dt = 0: Invalid_argument Stepper.fixed: need dt > 0 and t_end >= t0
```

The first two lines are runs with no way forward. `y' = y²` from `y(0) = 1` has the exact solution `1 / (1 - t)` and blows up at `t = 1`; in the second, the right-hand side is `nan` past `t = 0.5`, so every step that reaches beyond it fails in Newton. In both the controller halves the failed step until it is below `16 eps |t|`, where it gives up with `StepRejected n`, `n` being the rejections in a row (in the first, the very first rejection of the final streak already meets the floor; [chapter 5](05-step-control.md) derives both counts). Without that floor both runs would never end: the budget stops them and the lines show the budget text instead. The third line pins which step is halved: the first attempt is cut to the span `1`, and a rejection must halve that cut step, not the requested `1e30`, or the same failure repeats until `max_rejects` runs out (`Error StepRejected 51`). Line 17 says that a zero span is valid. The last three lines are the messages of `Check`: a span that runs backwards (line 18), a first step of zero (line 19) and a zero step for `Stepper.fixed` (line 20). Without the checks they print `Ok ...` or `Error StepRejected 51` instead (section 9).

## 9. The soak, and what each pin catches

```
soak canary x10: passed 10/10, identical: true
soak logistic order x10: passed 10/10, identical: true
soak van der Pol x10: passed 10/10, identical: true
soak robertson x10: passed 10/10, identical: true
```

[`test/soak.ml`](../../test/soak.ml) reruns the canary, the logistic order check, van der Pol and Robertson ten times each, with the same criteria as lines 6, 7, 8 and 9, and checks that every round passes and that all rounds are equal under `=`. The library has no hidden state, so `identical` can turn `false` only if someone adds hidden state, randomness or parallelism (or a NaN appears, since NaN is never equal to itself). It is a tripwire, not extra coverage. The mechanics are in [docs/testing.md](../testing.md).

Each row is one deliberate change that you can repeat in a scratch copy of the repository, running the built test programs under a time limit; [docs/testing.md](../testing.md) has the method and a longer table. "Silent" means no line of either expected file changes.

| Change | What the corpus does |
|--------|----------------------|
| `Jac.forward`: swap rows 0 and 2 of the result | both `jac canary` lines and line 11 print `false`; line 10 prints `Error Diverged`; line 6 still passes but `3.68e-07` becomes `3.83e-07` |
| `Jac.forward`: return the transpose | line 11 fails if reached, but the Robertson run no longer finishes in practice (van der Pol still does, with many more steps), so a full run needs a timeout |
| `Jac` step `1e-2` or `1e-14` instead of `1e-8` | line 11 prints `false` for both; line 5 also for `1e-14`; line 23 moves its digits for `1e-2` |
| `Linalg`: pivot always on row 0 | line 12 prints `None`, line 13 prints `[0.000000000000; 1.000000000000]` |
| Newton tolerance `1e-6` instead of `1e-10` | only the 12th decimal of line 2 |
| Newton: no stopping test on the step | line 2 and every line that takes an integration step print `Error Diverged` or `Error StepRejected 51` |
| Newton: the step test applied only after the line search accepts a step | lines 6, 7 and 10 print `Error Diverged`; lines 8, 16, 21 and 22 print `Error StepRejected`, lines 14 and 15 other counts; the soak lines for the canary, the logistic order and van der Pol pass 0/10 |
| `Halving`: double after two accepts, not three | lines 21 and 22, and the count of the NaN wall (line 15) |
| `Bdf2`: startup estimate factor `1e-3` instead of `0.5` | line 21 |
| `Halving`: no cap by `dt_max` | line 22 |
| `Halving`: error weight 1 instead of `1 + abs(y_i)` | the blow-up and NaN-wall lines (14 and 15), lines 21 and 23 |
| `Halving`: halve the requested step, not the failed one | the `1e30` line (16): `Error StepRejected 51` |
| `Halving`: no `16 eps abs(t)` floor | the blow-up and NaN-wall lines (14 and 15) print the budget line |
| `Adaptive`: do not shorten the last step | the `1e30` line (16) and line 23 |
| `Check`: remove the argument checks | the three `Invalid_argument` lines (18 to 20) print `Ok ...` or `Error StepRejected 51` |
| `Bdf2`: return backward Euler's result, not BDF2's | line 23 prints `1.6e-04 < 1e-5: false` and line 21 changes, while line 9 still passes |
| `Bdf2.coeffs`: `beta` off by 1 part in `1e5` times `(ω - 1)` | the digits of line 21 and the count of the blow-up line |
| `Halving`: a rejection does not reset the accept streak | line 21: the `rejected` count changes |
| Newton: iteration limit 1 | the run does not finish in practice |

## 10. What the corpus does not catch

These changes leave every line of both expected files unchanged. Each suggested case gives a different result on the broken version (the Armijo constant is the one exception, as noted).

- **Newton's safeguards.** Without damping (`min_damping = 1.`), with a NaN trial residual accepted, or with an iteration limit of 100 instead of 50, nothing changes. Cases from [chapter 2](02-newton.md): `atan` from 2 and `sqrt x - 1` from 9 fail without damping; `sqrt x - 1` from 9 ends in `Error Nan` if a NaN trial point is accepted; `x³` from 1 is `Error Diverged`, and `Ok` with a limit of 100; a start at `nan` is `Error Nan`, and `Error Diverged` without the check. None of these cases tells the Armijo constant from `0`. The finiteness test on the trial residual is redundant (`nan <= bound` is already false, [chapter 3](03-jacobians-and-floating-point.md)), and removing the one on the step changed no result either. The order of the step test and the line search is not a gap: reversing it moves many lines (section 9).
- **Jacobian details.** Neither the relative factor in the step nor the division by the stored perturbation is pinned. Cases: `Jac.forward (fun y -> y) [| 1e10 |]` must be `[[1]]` up to rounding (an absolute step gives `nan`, the nominal divisor `1 - 1e-10`), and at `[| 123.456 |]` the entry is off by about `3e-9` with the nominal divisor.
- **The time argument.** Every problem ignores `t` (the NaN wall only uses it as a threshold), so passing `t_n` instead of `t_{n+1}` to the right-hand side in the stage equation changes no line. Case: `y' = -y + t`, `y(0) = 1`, exact `t - 1 + 2 e^-t`: halving `dt` divides BDF2's error by about 4, and by only about 2 with that mistake (first order); the program after this list measures it.
- **Cost.** Starting Newton from the previous state instead of the extrapolated guess changes the number of right-hand-side calls but no line. A case could count them with `Instrument.count`.
- **Step-control defaults and limits.** The default `dt_max` (`span / 10`), the `max_rejects` give-up (the NaN wall ends by the floor, at 46 rejections, below the default 50), accepting a NaN error estimate, and the `Error Nan` check on `y0` and `rhs t0 y0` are all silent. Cases: `y' = 0` on `[0, 1]` accepts every step, so the count follows from the rules alone, like line 22 (`dt0 = 1e-6` doubles every third accept up to `0.1`: `3 · 17 = 51` steps cover `0.39`, and 7 more finish, 58 in all; another `dt_max` or doubling rule gives another count); `~max_rejects:3` on the NaN wall gives `Error StepRejected 4`; `Halving.acceptable` with `err = [| nan |]` is `false`; a right-hand side that is `nan` at `t0` gives `Error Nan`.
- **The fixed-step driver.** That `Stepper.fixed` takes at least one step and that its equal steps cover the span is silent. Cases: `dt = 4` on the canary's `[0, 1]` must take one step of size 1 (without that rule the state comes back unchanged); `dt = 0.3` must take three steps of `1/3`, not steps of `0.3`.
- **Vectors.** `Vec` has no case of its own. `Vec.norm_inf [| nan; 1. |]` should be `nan`. A version built on the generic `max`, or on `if Float.abs v > m then Float.abs v else m`, returns `1.` (with the NaN last, `[| 1.; nan |]`, the generic `max` would still pass).
- **What the lines cannot see.** Van der Pol is not compared with a reference; only `y_1` of Robertson is; a transposed Jacobian stalls the run instead of failing one line; the soak is a tripwire; `Instrument` counts right-hand-side calls only, not Newton iterations. The `500,000`-step canary dominates the run time.

The case for the time argument, as a program:

```ocaml
open Vstiff

let () =
  let problem = { Ode.rhs = (fun t y -> [| -.y.(0) +. t |]); t0 = 0.; t_end = 1.; y0 = [| 1. |] } in
  let exact = 2. *. exp (-1.) in
  let error dt =
    match Stepper.fixed (module Bdf2) ~dt problem with
    | Ok y -> Float.abs (y.(0) -. exact)
    | Error _ -> Float.nan
  in
  Printf.printf "BDF2 error ratio on halving dt: %.1f\n" (error 0.01 /. error 0.005)
```

```text
BDF2 error ratio on halving dt: 4.0
```

## In the code

| Idea | Where |
|------|-------|
| The four problems: `rhs`, `y0`, `problem`, and `exact` for the first two | `Canary`, `Logistic`, `VanDerPol`, `Robertson` in [`test/problems.ml`](../../test/problems.ml) |
| The external Robertson reference | `Refs.robertson_y1_at_1e4` in [`test/refs.ml`](../../test/refs.ml), provenance in its doc comment |
| The cases | the named lists of [`test/corpus.ml`](../../test/corpus.ml), mapped to line numbers in section 1; `corpus` concatenates them in output order |
| Test effects | `Guard.budget` and `Guard.run` in [`test/guard.ml`](../../test/guard.ml); `Report.lines` in [`test/report.ml`](../../test/report.ml) |
| The expected output | [`test/corpus.expected`](../../test/corpus.expected) (23 lines), [`test/soak.expected`](../../test/soak.expected) (4 lines) |
| The soak | [`test/soak.ml`](../../test/soak.ml): the module type `Case`, `repeat`, `soak` and four case modules |
| How both run | [`test/dune`](../../test/dune): `(tests (names corpus soak) (libraries vstiff))` |

## Check yourself

1. **With two Jacobian rows swapped, which corpus lines change, and which changed lines still say `true`?**
   Both `jac canary` lines and line 11 turn `false`, and line 10 becomes `Error Diverged`. Line 6 changes from `3.68e-07` to `3.83e-07` and stays `true`: the roots of the stage equation do not depend on the Jacobian, and the slower iteration only stops at slightly different points.

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
   Passing `t_n` instead of `t_{n+1}` to the right-hand side: every problem ignores `t`. The case `y' = -y + t` shows the error ratio falling from about 4 to about 2.

That is the end of the numerical-methods track. Next: [docs/testing.md](../testing.md) for running the tests and changing expectations, and [docs/exercises.md](../exercises.md) for exercises and starter contributions.
