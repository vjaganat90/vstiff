# 5. Error estimates and step-size control

Chapter 5 of 6 in the numerical-methods track. Previous: [4. BDF methods](04-bdf.md). Next: [6. The test corpus](06-the-corpus.md). Index: [docs/README.md](../README.md). Terms and symbols: [glossary](../glossary.md).

**Summary.** An integrator has to estimate its own error to choose its steps. The idea is a pair of methods of different order whose gap estimates the error of the lower one; here the pair is backward Euler and BDF2, and `Bdf2.step_with_error` returns the gap. A controller turns each estimate into a decision: `Halving` accepts a step when the scaled estimate is at most `tol`, halves the step that failed, doubles after three accepts in a row, and gives up in two stated cases. This chapter derives the estimate (including the special one for the first step), spells out the rules exactly as the code implements them, shows how `Adaptive.integrate` drives any controller through the `Ode.Controller` contract, and ends with what to expect on van der Pol and how production controllers differ.

Prerequisites: [chapter 1](01-odes-and-stiffness.md) (local error, order) and [chapter 4](04-bdf.md) (BDF1, BDF2, the history, the step ratio `ω`). The probes use the probe project (set up on day 1 from Setup in [exercises.md](../exercises.md), or "Running the snippets" in chapter 1).

## 1. Why estimate the error?

A fixed step has to be chosen in advance, and choosing it well requires knowing the solution. The corpus does that by hand for smooth problems (`dt = 0.01` for the logistic order test, `dt = 2e-6` for the canary). Many problems have phases that need very different steps. Van der Pol with `μ = 1000` drifts slowly, where long steps are fine, then jumps abruptly on a time scale of about `1/μ`, where tiny steps are needed. A step small enough for the jumps wastes effort on the drift; a step suited to the drift misses the jumps. The cure is to adapt: measure how big the error of each step is, shrink the step when it is too big and grow it when it is comfortably small. We cannot measure the error, because we do not know the exact solution. So we **estimate** it from quantities we do have. Each try of a step is an **attempt**; it is either accepted or rejected.

## 2. What is controlled

The estimate is of the **local error** from chapter 1: the error one step commits. A step is kept when the estimate, scaled as in section 6, is at most `tol`, the tolerance the caller passes to `Adaptive.integrate`. That is standard practice in ODE codes, and it has a consequence: the final error is the accumulation of many local errors, each carried along by the dynamics of the problem (damped if it decays, amplified if it is unstable), so `tol` is not a bound on the final error.

How does the final error relate to `tol`? A heuristic: the estimate behaves like `h²` (section 4), so keeping it near `tol` makes `h` proportional to `sqrt(tol)`. One BDF2 step's true error is only `O(h³)`, but there are about `T / h` steps, so the global error is `O(h²)`, which is proportional to `tol` times a constant that depends on the problem and the interval. The work grows like `1 / h`, so for small `tol` every factor of 100 taken off `tol` multiplies the number of steps by about 10. Nothing guarantees this when errors grow along the solution. The probe prints the final error next to `tol`, with the numbers of accepted and rejected attempts (section 8), so you can check the prediction:

```ocaml
open Vstiff

let () =
  let open Problems.Logistic in
  List.iter
    (fun tol ->
      match Adaptive.integrate (module Bdf2) (module Halving) ~tol problem with
      | Ok s ->
          Printf.printf "tol = %-6g final error = %.2e   (%d accepted, %d rejected)\n%!" tol
            (Float.abs (s.y.(0) -. (exact 5.).(0)))
            s.stats.accepted_steps s.stats.rejected_steps
      | Error e -> print_endline (Fail.to_string e))
    [ 1e-3; 1e-4; 1e-5; 1e-6 ]
```

`Adaptive.integrate (module Bdf2) (module Halving)` picks the method and the controller; `s.y` and `s.stats.accepted_steps` read fields of the solution record. Tolerances near or below the accuracy of the arithmetic or of the nonlinear solves (`eps`, called `Float.epsilon` in the code, is about `2.2e-16`; Newton's tolerance is `1e-10`) can cost an enormous number of steps without delivering the accuracy asked for: add `1e-10`, `1e-12` and `1e-14` to the list and compare. The last run takes much longer, so run the probe under a time limit, for example `perl -e 'alarm 120; exec @ARGV' ./_build/default/p/probe.exe`; the `%!` at the end of the format flushes the output, so the finished lines appear even if the limit kills the program.

## 3. The pair idea

Suppose two methods L and H take one step from the same data, and L has order `p` while H has order `p + 1`. Let `y(t + h)` be the exact value that step is trying to reach. Then, to leading order,

```
y(t + h) - y_L = c h^(p+1) + O(h^(p+2))          (local error of L)
y(t + h) - y_H = O(h^(p+2))                       (local error of H, one power smaller)
```

Subtracting, `y_H - y_L = c h^(p+1) + O(h^(p+2))`: **the gap between the two results is the local error of the lower-order method**, to leading order, even though we never see the exact value. Keep the more accurate result H and use the gap as the estimate. Keeping H while the gap measures its weaker partner L is called local extrapolation, and the estimate then overstates the error of what we kept.

In Runge-Kutta codes (the other big family of ODE methods, not part of this library) such pairs are called *embedded*, because the two methods share most of their work. Here they share nothing, although the contract is still called `Ode.Embedded`: each attempt solves two nonlinear stage equations, one for backward Euler and one for BDF2. That is simple, because both methods are already in the code, but it is not cheap.

## 4. The BDF1/BDF2 gap

`Bdf2.step_with_error rhs h history at` is the one function `Ode.Embedded` adds to `Ode.Method`. It returns `(y, err, next)`: the new state, an estimate of the step's local error and the next history. Here it is with its comments removed ([`src/bdf2.ml`](../../src/bdf2.ml)):

```text
let step_with_error rhs h history at =
  let* be = backward_euler rhs h at in
  let next = After { h_prev = h; y_prev = at.y } in
  match history with
  | Start -> Ok (be, Vec.scale 0.5 (Vec.sub be (Vec.axpy h (rhs at.t at.y) at.y)), next)
  | After { h_prev; y_prev } ->
      let+ y = bdf2 rhs h ~h_prev ~y_prev at in
      (y, Vec.sub y be, next)
```

`let* x = e in body` is `Result.bind e (fun x -> body)` and `let+` is the `Result.map` form (both from `Fail.Syntax`, [`src/fail.mli`](../../src/fail.mli); [docs/ocaml.md](../ocaml.md) explains them): if `e` is an `Error`, the whole function returns that error at once. So if either Newton solve fails, the attempt is an `Error`. `backward_euler` is `Bdf1.step` from the same point. The `After` branch, discussed now, keeps the BDF2 result `y` and reports the gap `y - be` as the estimate. The `Start` branch is the first step (section 5).

Take the local errors of backward Euler (BE) and of BDF2 (chapter 4), both from the same exact data, with all derivatives at `t_n`. Section 5 derives the BE formula; chapter 4 checked it on `y = t²`:

```
exact - BE   = -(h²/2) y'' + O(h³)
exact - BDF2 =  C(ω) h³ y''' + O(h⁴)
```

so `y - be = (y - exact) + (exact - be) ≈ -(h²/2) y''`: the gap is backward Euler's local error, of size `(h²/2) |y''|`. The kept result, BDF2, has a local error `O(h³)`. For small `h` the real error of the step is therefore much smaller than the estimate. The estimate is **conservative**, and increasingly so as `h` shrinks: when `h` is halved (with `ω` fixed), the gap falls by about a factor of 4 and the true BDF2 error by about a factor of 8. These expansions assume a small `h`; for a step that is long compared with the fastest time scale of the problem, the estimate is only a rough guide.

**Probe: the gap, and the startup estimate.** The first block uses exact data on the logistic equation at `t = 1`. `Bdf2.step` from `Bdf2.start`, begun at the exact previous point, leaves a history whose `y_prev` is exact (the state it computes is discarded); `step_with_error` then runs from the exact current point, and the gap is compared with the true local error of BDF2. The second block checks the startup estimate of section 5.

```ocaml
open Vstiff

let () =
  let open Problems.Logistic in
  let t = 1. in
  List.iter
    (fun h ->
      let _, history = Result.get_ok (Bdf2.step rhs h Bdf2.start { Ode.t = t -. h; y = exact (t -. h) }) in
      let y, err, _ = Result.get_ok (Bdf2.step_with_error rhs h history { Ode.t; y = exact t }) in
      let true_error = Float.abs (y.(0) -. (exact (t +. h)).(0)) in
      Printf.printf "h = %-6g gap = %.3e   true BDF2 error = %.3e   gap / error = %.0f\n" h (Float.abs err.(0))
        true_error (Float.abs err.(0) /. true_error))
    [ 0.2; 0.1; 0.05; 0.025 ];
  List.iter
    (fun h ->
      let y, err, _ = Result.get_ok (Bdf2.step_with_error rhs h Bdf2.start { Ode.t = 0.; y = y0 }) in
      Printf.printf "h = %-6g estimate = %.4e   backward Euler error = %.4e\n" h (Float.abs err.(0))
        (Float.abs (y.(0) -. (exact h).(0))))
    [ 0.1; 0.05; 0.025 ]
```

OCaml notes: `{ Ode.t = t -. h; y = ... }` builds an `Ode.point` (the `Ode.` prefix says which module the field names belong to, and a bare `t` means `t = t`); `Result.get_ok` unwraps an `Ok` and raises on `Error`, which is fine in a probe. In the first block, each halving of `h` should shrink the gap by about 4 and the true error by about 8, so the last column grows roughly like `1/h`.

## 5. The first step: the startup estimate

With no history there is no BDF2. From `Start`, `step_with_error` returns backward Euler's result and, for its error, half the gap between backward Euler and the explicit Euler value `y + h f(t, y)` (the expression `Vec.axpy h (rhs at.t at.y) at.y`, one extra `rhs` call). Taylor series explain the `0.5`. Let `y`, `y'` and `y''` be evaluated at `t`:

```
exact:             y(t + h) = y + h y' + (h²/2) y'' + O(h³)
explicit Euler:    y + h y'                                  (error  +(h²/2) y'')
backward Euler:    y + h f(t + h, y_BE) = y + h y' + h² y'' + O(h³)
```

The last line uses `f(t + h, y_BE) = y'(t + h) + O(h²) = y' + h y'' + O(h²)`, which holds because backward Euler's own error is `O(h²)`. Subtracting,

```
y_BE - y_explicit = h² y'' + O(h³)          while the local error of backward Euler is -(h²/2) y''
```

so the gap is twice backward Euler's local error in magnitude, and half the gap estimates it. The sign comes out opposite, which does not matter because the controller takes absolute values: the two methods err by equal and opposite amounts at leading order, and [chapter 1](01-odes-and-stiffness.md) showed the opposite signs with numbers.

By hand on the logistic problem, `y_0 = 0.1`, `h = 0.1`. Explicit Euler gives `0.1 + 0.1 · 0.09 = 0.109`. Backward Euler solves `x = 0.1 + 0.1 x (1 - x)`, i.e. `x² + 9x - 1 = 0`, so `x = (-9 + √85) / 2 = 0.1097722`. The gap is `7.72e-4` and half of it is `3.86e-4`. The exact value is `y(0.1) = 0.1093669`, so backward Euler's true error is `4.05e-4`: the estimate is about 5% low. The leading-order formula gives `(h²/2) y''(0) = 0.005 · 0.072 = 3.6e-4` (with `y'' = y (1 - y) (1 - 2y) = 0.072` at `y = 0.1`). The neglected terms are one power of `h` smaller, so the relative error of the estimate shrinks like `h`; the second block of the probe in section 4 shows it. Unlike the BDF1/BDF2 gap, this estimate belongs to the very solution that is kept (backward Euler), so on the first step it is not conservative.

## 6. The scaled error and `tol`

The estimate is a vector, one entry per component. `Halving.acceptable` turns it into one number, using the candidate new state `y`, and compares it with `tol`:

```
max_i |err_i| / (1 + |y_i|)   <=   tol
```

Read the weight `1 / (1 + |y_i|)` as follows. Where `|y_i|` is much larger than 1 it is about `1 / |y_i|`, so `tol` acts as a *relative* tolerance. Where `|y_i|` is much smaller than 1 it is about 1, so `tol` acts as an *absolute* tolerance. The crossover is at `|y_i|` of order 1. The `1 +` also avoids dividing by zero when a component is exactly zero. In the vocabulary of production codes (`atol` is the absolute and `rtol` the relative tolerance) this is `|err_i| <= atol + rtol |y_i|` with `atol = rtol = tol`. The maximum means a single badly scaled component rejects the step.

By hand, for `y = (1000, 1e-5)` and `err = (0.5, 1e-6)`: the weights are `1/1001` and about 1, so the scaled entries are `4.995e-4` and about `1e-6`, and the norm is `4.995e-4`, from the first component. With `tol = 1e-3` that passes, although the second component is wrong by 10% of its own size. The weights treat a tiny component as nearly absolute, so it is controlled only loosely. That matters for Robertson's problem: `y2` never exceeds about `3.6e-5`, so with `tol = 1e-6` an estimated local error of `1e-6` in `y2`, about 3% of its peak and a much larger share of its value at other times, passes the test. There is one tolerance for all components and no separate `rtol` and `atol`; section 11 lists that as a limitation.

Last, NaN ("not a number") is what floating-point arithmetic gives for invalid operations such as `0. /. 0.` ([chapter 3](03-jacobians-and-floating-point.md)), and `nan <= x` is false whatever `x` is. An estimate that is NaN makes the maximum NaN, the comparison with `tol` fails, and the step is rejected, never accepted. `Halving.acceptable` pairs `err` with `y` entry by entry (`Array.map2`), so the two vectors must have the same length.

## 7. The controller contract and the driver

The method says how good a step was; a controller decides what to do about it. `Ode.Controller` ([`src/ode.mli`](../../src/ode.mli)) is the contract between the controller and the driver, shown here without its doc comments:

```text
module type Controller = sig
  type t
  type stats
  val init : tol:float -> dt0:float -> dt_max:float -> max_rejects:int -> t
  val proposal : t -> float
  val acceptable : t -> y:Vec.t -> err:Vec.t -> bool
  val accepted : t -> t
  val rejected : t -> rejection -> at:float -> h:float -> (t, Fail.t) result
  val stats : t -> stats
end
```

By the project's convention a controller is a pure state machine: each function returns a new state instead of changing the old one, and `rejected` can return an `Error` to give up. It is told why a step was rejected: `Too_large` (the estimate was above tolerance) or `Solver e` (the step itself failed with `e`). `Halving` ignores the reason; a controller could treat the two differently. Only `acceptable` is given the estimate `err`: `accepted` and `rejected` see just the outcome.

The driver, `Adaptive.integrate (module M) (module C) ?dt0 ?dt_max ?max_rejects ~tol problem` (`?` marks an optional argument, `~` a labelled one; [`src/adaptive.ml`](../../src/adaptive.ml)), starts from `(t0, y0)` with `M.start` and `C.init`, then repeats:

1. If `t >= t_end`, return `Ok { t; y; stats = C.stats c }`.
2. Ask `dt = C.proposal c`. If `dt >= t_end - t` this is the final step and `h = t_end - t`, so it lands exactly on `t_end`; otherwise `h = dt`.
3. Try `M.step_with_error problem.rhs h history at`. If it returns `Ok (y, err, next)` and `C.acceptable c ~y ~err`, **accept**: continue from `y` at `t + h` (after a final step the driver assigns `t_end` itself, because `t + h` could miss it by rounding), with `C.accepted c` and the history `next`.
4. Otherwise **reject**, with reason `Too_large` if the step was `Ok` but not acceptable and `Solver e` if it returned `Error e`: a failed step is handled exactly like a bad estimate. Call `C.rejected c reason ~at:t ~h`. On `Ok c'` retry from the *same* `t`, `y` and history; on `Error e` return `Error e`.

Before the loop, invalid arguments raise `Invalid_argument` (`t_end < t0`, or non-positive step bounds on a non-empty span; `Check.adaptive`), and a non-finite `y0` or `rhs t0 y0` gives `Error Nan`. `dt_max` defaults to `(t_end - t0) / 10`, the first proposal is `min dt_max dt0` with `dt0` defaulting to `1e-6 (t_end - t0)`, and `max_rejects` defaults to 50. The result type `(C.stats solution, Fail.t) result` depends on the controller passed in, which is why the controller is a modular explicit argument and not a plain value ([docs/ocaml.md](../ocaml.md)); the method is passed the same way.

## 8. Halving's rules

`Halving` ([`src/halving.ml`](../../src/halving.ml)) is the controller the library ships. Its state, a record hidden by the interface, holds the tolerance and the bounds, `dt` (the step to propose next), `streak` (accepts counted towards the next doubling: 0, 1 or 2), `failures` (rejections in a row) and the stats, `{ accepted_steps; rejected_steps }`. The rules:

- **`acceptable`**: the scaled error of section 6 is at most `tol`.
- **`accepted`**: `streak` increases by one. When it reaches 3, `dt` becomes `min (2 dt) dt_max` and `streak` goes back to 0; otherwise `dt` is unchanged. `failures` goes back to 0 and `accepted_steps` increases.
- **`rejected`**, for either reason: `failures` increases by one and the new `dt` is `h / 2`, half of the step that *failed*: if the driver had cut that step to land on `t_end`, this is less than half of the old `dt`. If `failures > max_rejects` or the new `dt` is below `16 eps |at|`, it returns `Error (StepRejected failures)` (section 9). Otherwise `streak` goes back to 0, `rejected_steps` increases and the driver retries.

Consequences worth noticing: steps change only by factors of 2 (apart from the cut that lands on `t_end`), which is what keeps the step ratio `ω` at most 2 (chapter 4); apart from the cut step and its retries, each step size is the starting step size or `dt_max` times a power of two; after any change of `dt`, three accepted steps are needed before the next doubling; and neither growth nor shrink depends on how far the estimate is from `tol`, which decides only accept or reject.

**A trace.** Drive `Halving` by hand: steps up to `5e-3` pass the test, longer ones fail, and `dt0 = 1e-3`. (In the probe, `let rec go ... in` defines a function that calls itself, which is how OCaml writes a loop; `begin ... end` groups statements like braces; `match` takes a `result` apart. [docs/ocaml.md](../ocaml.md) has the details.)

```ocaml
open Vstiff

let () =
  let rec go c attempts =
    let dt = Halving.proposal c in
    if attempts > 0 then
      if dt <= 5e-3 then begin
        Printf.printf "dt = %-6g accept\n" dt;
        go (Halving.accepted c) (attempts - 1)
      end
      else begin
        Printf.printf "dt = %-6g reject\n" dt;
        match Halving.rejected c Ode.Too_large ~at:1. ~h:dt with
        | Ok c -> go c (attempts - 1)
        | Error e -> print_endline (Fail.to_string e)
      end
  in
  go (Halving.init ~tol:1e-6 ~dt0:1e-3 ~dt_max:1. ~max_rejects:50) 14
```

The 14 lines it prints condense to this (`A` is an accept, `R` a reject):

```text
dt = 1e-3:  A A A    three accepts, so dt doubles to 2e-3
dt = 2e-3:  A A A    doubles to 4e-3
dt = 4e-3:  A A A    doubles to 8e-3
dt = 8e-3:  R        too long: halves to 4e-3 and the streak restarts
dt = 4e-3:  A A A    doubles to 8e-3 again
dt = 8e-3:  R        ... and so on: A A A R, A A A R
```

Whenever the best step lies between `dt` and `2 dt`, one attempt in four is a rejection in this steady state: the price of probing for longer steps. It is not a sign of trouble.

## 9. Giving up

The `rejected` rule has two exits to `Error (StepRejected n)`, where `n` is the length of the current run of consecutive rejections.

- **Too many in a row.** `failures > max_rejects` (default 50), so the 51st rejection in a row ends the run. Fifty halvings shrink a step by `2^50`, about `1e15`, so this only trips when no reasonable step size helps.
- **The floor.** The new `dt` is below `16 eps |at|`, with `at` the time we are stepping from. Floats near `at` are one **ulp** (unit in the last place) apart, which is between `eps |at| / 2` and `eps |at|` ([chapter 3](03-jacobians-and-floating-point.md) explains floating point), so the floor is 16 to 32 ulps. A step that short moves `t` by only 16 to 32 representable values, so `t + dt` is resolved to a few percent at best, and halving again cannot get past whatever caused the trouble: the solution is singular there, or the right-hand side fails just ahead.

At `t = 0` the floor is zero and `dt < 0` is never true, so only `max_rejects` can stop the run. A run that cannot move from `t0 = 0` ends with `StepRejected 51`. From `t0 = 1` with `t_end = 2`, the default first step is `1e-6` and the floor is `16 eps = 3.6e-15`. After `k` halvings the step is `1e-6 / 2^k`, which first drops below the floor at `k = 29` (`2^28 = 2.7e8` is below `1e-6 / 3.6e-15 = 2.8e8`, and `2^29` is above it), so a run that cannot move ends with `StepRejected 29`.

`StepRejected n` can have a small `n`, because an accepted step resets the count. Take the blow-up `y' = y²`, `y(0) = 1`, whose solution `1 / (1 - t)` is unbounded at `t = 1`, so `t_end = 2` cannot be reached. Here `y'' = 2 y³`, so once `y` is large the scaled estimate is about `(h y)²`: the accepted steps shrink like `1 / y` as the computed solution grows without bound, until a rejection whose halved step falls below the floor ends the run. The corpus blow-up ends with `StepRejected 1`: the attempt before that rejection was accepted. A wall is different. In the corpus the right-hand side turns NaN past `t = 0.5` (with `dt0 = 0.25`); the steps are powers of two, so the run reaches `t = 0.5` exactly, every attempt from there fails and all the halvings come in a row. The floor there is `16 eps · 0.5 = 2^-49`, and exactly 46 halvings take a step below it when the first failing step lies in `[2^-4, 2^-3)`, which is what the corpus line's `StepRejected 46` says about this run. Either way the error carries only `n`; the time reached is lost.

**Newton failures.** Treating a failed stage solve as a rejection makes sense because the stage equation `x = psi + gamma f(t_{n+1}, x)` gets easier as `h` shrinks: `gamma` shrinks, the Jacobian `I - gamma J` gets closer to the identity, and the initial guess gets closer to the answer. Sometimes there is no solution at all. For `y' = y²` from `y = 1`, backward Euler solves `x = 1 + h x²`, which has a real root only for `h <= 1/4`: at `h = 0.5` the equation `0.5 x² - x + 1 = 0` has discriminant `1 - 2 < 0`, so there is no root to find, `Newton.solve` returns an `Error`, and the controller has to shorten the step. Non-finite values inside a step are caught by Newton (`Error Nan` or `Diverged`) and turned into rejections too.

**Probe: a wall at the starting time.** It leaves no room to move, so the floor and the limit decide:

```ocaml
open Vstiff

let () =
  let wall w t y = if t > w then [| nan |] else [| -.y.(0) |] in
  List.iter
    (fun t0 ->
      let problem = { Ode.rhs = wall t0; t0; t_end = t0 +. 1.; y0 = [| 1. |] } in
      match Adaptive.integrate (module Bdf2) (module Halving) ~tol:1e-4 problem with
      | Ok s -> Printf.printf "wall at t0 = %g: Ok at t = %g\n" t0 s.t
      | Error e -> Printf.printf "wall at t0 = %g: Error %s\n" t0 (Fail.to_string e))
    [ 0.; 1. ]
```

Expect `StepRejected 51` and `StepRejected 29`, the numbers derived above. (`wall t0` is a partial application: a function still waiting for `t` and `y`; `-.` negates a float.) The blow-up and the NaN wall are in the corpus ([`test/corpus.ml`](../../test/corpus.ml)).

## 10. What to expect on van der Pol

With `μ = 1000` and `y(0) = (2, 0)`, van der Pol is a relaxation oscillation. `y1` drifts slowly while `y2` is small, then jumps on a time scale of about `1/μ`; leading-order theory puts the jumps near `t = 807` and `t = 1614`, and `|y2|` reaches about 1333 during a jump. At the initial state the Jacobian already has an eigenvalue near `-3000` (chapter 1), so the problem is stiff from the start. The corpus runs it on `[0, 2000]` with `tol = 1e-4`, and the expected output records only:

```text
vdp mu=1000 [0,2000] tol=1e-4: Ok at t=2000, rejected >= 1: true, y finite: true
```

In this run `dt_max` defaults to `200` and the first step to `2e-3`. During a drift the controller can take long steps. At a jump it has to shrink by orders of magnitude, one halving per rejection: cutting a step by a factor `F` takes about `log2 F` consecutive rejections, which stays within the limit of 50 for any `F` up to `2^50`. Afterwards it grows back by doubling every three accepted steps, so recovering a factor `F` takes about `3 log2 F` accepted steps. On top of that come the routine rejections of section 8: one attempt in four whenever the best step is between `dt` and `2 dt`. Many rejections are normal with this controller, because the doubling rule keeps probing larger steps; the corpus only demands at least one. To see the proportion yourself, adapt the probe of section 2: pass `Problems.VanDerPol.problem` with `~tol:1e-4`, drop the error column (there is no exact solution) and print the two counts.

## 11. How production controllers differ

This controller is deliberately simple: the factor 2 and the count 3 are the only numbers in its step-size rule. Production codes usually do more, and each difference is a possible contribution (the starter list is in [docs/exercises.md](../exercises.md)):

- **Error-proportional step choice.** If the estimate behaves like `C h^(p+1)`, where `p` is the order of the method whose error it measures, the step that would bring it exactly to `tol` is `h (tol / err)^(1/(p+1))`. Production codes use this to pick the next step directly, multiplied by a safety factor somewhat below 1 (0.8 to 0.9 is common) and clamped so the step cannot change by too large a factor at once. For the BDF1/BDF2 gap `p = 1`, so `dt_new = h · safety · sqrt(tol / err)`. Steps then settle near the largest acceptable value instead of cycling between `dt` and `2 dt`. This does not fit `Ode.Controller` as it stands, because `accepted` and `rejected` are not given `err`: the contract and the driver would have to change. Whatever the clamp on growth, it must keep `ω` below `1 + √2` (chapter 4).
- **Separate `rtol` and `atol`**, per component if needed, with weights `atol_i + rtol |y_i|`. The single weight `1 / (1 + |y_i|)` here cannot control tiny components such as Robertson's `y2` tightly without making everything else too strict.
- **Cheaper estimates.** The BDF1/BDF2 pair costs two nonlinear solves per attempt (one on the first step). Production BDF codes typically estimate the error from a single solve, for example from the difference between the converged value and a predictor of matching order (a predictor extrapolates the past values to the new time, as the Newton guess of chapter 4 does). The BDF1/BDF2 gap is also conservative for BDF2, as section 4 showed.
- **Telling the failures apart.** `Halving` counts a too-large estimate and a failed Newton solve together; a controller could count Newton failures separately.
- **Reuse of Jacobians.** Here Newton rebuilds the finite-difference Jacobian at every iteration, at one `rhs` call per component plus one each time (`Instrument.count` wraps an `rhs` to count its calls). Production codes keep one Jacobian across iterations and steps ([chapter 2](02-newton.md)).
- **Variable order and dense output.** Production codes typically choose the order automatically (up to 5 for BDF) and interpolate between steps to give values at requested times. Here the integrators return only the final state.

## In the code

| Idea | Where |
|------|-------|
| The step and its estimate; the controller contract | `Bdf2.step_with_error` in [`src/bdf2.ml`](../../src/bdf2.ml); the contracts `Ode.Embedded`, `Ode.Controller` and `Ode.rejection` (`Too_large`, `Solver of Fail.t`) in [`src/ode.mli`](../../src/ode.mli) |
| The accept and reject rules, the scaled error, the floor `16 eps abs(at)` | `Halving.acceptable`, `accepted`, `rejected` and the internal `dt_min` in [`src/halving.ml`](../../src/halving.ml) |
| The driver: loop, final step, defaults | `go` inside `Adaptive.integrate` in [`src/adaptive.ml`](../../src/adaptive.ml); argument checks in `Check.adaptive` ([`src/check.ml`](../../src/check.ml)) |
| The result and the failures | `Adaptive.solution` is `{ t; y; stats }`; for `Halving`, `stats` is `{ accepted_steps; rejected_steps }`. `Fail.StepRejected of int` is in [`src/fail.mli`](../../src/fail.mli); a failed stage solve is `Error Diverged` or `Error Nan` from `Newton.solve` ([`src/newton.ml`](../../src/newton.ml)) through `Stage.solve` |
| The corpus runs | The `vdp`, `robertson` and `adaptive ...` lines of [`test/corpus.ml`](../../test/corpus.ml); [`test/soak.ml`](../../test/soak.ml) reruns the van der Pol and Robertson runs ten times |

**What the tests pin.** Several rules are tied to corpus lines, and a change to them changes those lines ([`test/corpus.expected`](../../test/corpus.expected)). Doubling after two accepts instead of three changes both logistic step-count lines (`dt0=0.5` and `dt_max=1e-3`) and the NaN-wall count. Dropping the floor of section 9 makes the blow-up and NaN-wall runs never end (the corpus prints its call-budget line instead). Halving the current proposal instead of the step that failed makes the `dt0 = dt_max = 1e30` run fail with `StepRejected 51`. A startup-estimate factor of `1e-3` instead of `0.5` changes the `dt0=0.5` line, and removing the `dt_max` cap changes the `dt_max=1e-3` line. Keeping backward Euler's result instead of BDF2's moves the Robertson error from `7.2e-07` to `1.6e-04`, which fails the `1e-5` accuracy line although it would pass the `1e-3` acceptance bound. Removing the argument checks changes the three `Invalid_argument` lines. Two things stay loose: the van der Pol line checks only that the run ends at `t = 2000` with at least one rejection and a finite state, and van der Pol accuracy is not compared with a reference. The section Judging test strength by mutation of [docs/testing.md](../testing.md) lists more such changes and the lines that catch them.

## Check yourself

1. **What does `tol` bound, and is it a guarantee on the final error?**
   The scaled estimate of one step's local error. It is not a bound on the final error, which accumulates over steps.

2. **Why is `y - be` an estimate of backward Euler's local error, and why is it conservative for BDF2?**
   BDF2's local error is `O(h³)` and backward Euler's is `-(h²/2) y''`, so the gap is dominated by backward Euler's error. BDF2's real error is smaller by a factor that grows like `1/h`.

3. **Where does the factor `0.5` in the first-step estimate come from?**
   Explicit and backward Euler have equal and opposite local errors `±(h²/2) y''`, so their gap is `h² y''`, twice either error.

4. **After three accepts the step doubles and the next attempt fails. What do `proposal` and the streak look like afterwards?**
   `proposal` is back to the value it had before the doubling (unless `dt_max` cut the doubling short) and the streak is 0, so three more accepts are needed before the next doubling.

5. **A Newton solve fails during an attempt. What happens?**
   `step_with_error` returns an `Error`, the driver calls `rejected` with `Solver e`, and `Halving` handles it like a too-large estimate: the step halves and is retried from the same state and history.

6. **Why is there a floor of `16 eps |at|`, what happens at `t = 0`, and how can the result be `StepRejected 1`?**
   Steps of 16 to 32 ulps cannot get meaningfully smaller. At `t = 0` the floor is zero, so only `max_rejects` can end the run (`StepRejected 51` by default). `n` counts rejections in a row, and an accepted step resets it: if the steps shrank through accepted steps, the first rejection that falls below the floor ends the run with `n = 1`.

7. **Why are rejections normal on van der Pol, and what does the corpus check?**
   When the best step lies between `dt` and `2 dt`, the doubling rule produces the cycle accept, accept, accept, reject. The corpus checks only at least one rejection, `Ok` at `t = 2000` and a finite state.

Next: [6. The test corpus](06-the-corpus.md), the test problems and what each one proves.
