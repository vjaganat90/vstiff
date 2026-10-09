# Architecture

vstiff solves ordinary differential equations $y' = f(t, y)$: the state $y$ (a vector of numbers) changes with time $t$ at the rate $f(t, y)$, and the library computes $y$ at later times from $y$ at the start. It aims at *stiff* problems, where fast and slow changes are mixed ([numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md)). This page is the map of the code: the layers, the two libraries and the public API, the contracts the modules meet, what depends on what, what happens during one step, where effects live, which properties the code relies on, and where to edit. The mathematics is in the chapters under `numerics/`, the OCaml in [ocaml.md](ocaml.md), the tests in [testing.md](testing.md), terms and symbols in [glossary.md](glossary.md).

In one paragraph: an **implicit** step cannot compute the new state directly, because the new state appears on both sides of the *stage equation* $x = \psi + \gamma\thinspace f(t_{n+1}, x)$. `Stage` solves it with damped Newton's method ([numerics/02-newton.md](numerics/02-newton.md)), which needs a Jacobian (the matrix of partial derivatives of $f$). `Jac` builds it from differences of $f$, so the caller supplies only $f$, and the linear systems inside Newton are solved by Gaussian elimination. Newton, `Jac`, `Linalg` and `Vec` form a numerical kernel, a library of its own that knows nothing about ODEs. A *method* (`Bdf1`, `Bdf2`: backward differentiation formulas, [numerics/04-bdf.md](numerics/04-bdf.md)) turns a point into the next one, a *controller* (`Halving`) chooses step sizes, and two *drivers* (`Stepper.fixed`, `Adaptive.integrate`) repeat steps until the end time. Everything is pure except a few small modules (see "Effects").

## Layers

```text
 tests       Corpus  Soak  Props   shared: Problems  Refs  Gen  Prop  Dd   effects: Guard  Report
 ------------------------------------------------------------------------------
 vstiff      the solver library; Vstiff (src/vstiff.ml) lists the public modules
   drivers     Stepper.fixed (any Method)       Adaptive.integrate (Embedded + Controller)
   contracts   Ode: rhs, problem, point, rejection; Method, Embedded, Controller
   plug-ins    Bdf1 (Method)   Bdf2 (Embedded)   Halving (Controller)
   stage       Stage.solve (internal)     the stage equation of an implicit step
   support     Clock   Check (internal)   Instrument   Fail (the kernel's, re-exported)
 ------------------------------------------------------------------------------
 numerics    the kernel, which knows nothing about ODEs
               Newton.solve    Jac.forward    Linalg.solve    Vec    Fail (failures, let*, let+)
```

Read it upwards from the numerics row, the kernel: `Vec` and `Linalg` do arithmetic on arrays, `Newton` solves nonlinear equations with them and `Jac` supplies its Jacobians, and `Fail` names how that can go wrong. Above it the solver begins: `Stage` turns a BDF step into such an equation, `Clock` tells the controller, the adaptive driver and `Check` how short a step $t$ can still resolve, the methods set up the equation for their formula, and the drivers decide how steps are chained. The drivers never name a method or a controller: the caller passes them in as modules.

## The two libraries and the public API

The code is built as two dune libraries ([ocaml.md](ocaml.md) explains the mechanics with a miniature):

| Library | Defined in | Module | It holds |
|---|---|---|---|
| `numerics`, the kernel | [src/numerics/dune](../src/numerics/dune) | `Numerics` | `Fail`, `Vec`, `Linalg`, `Newton`, `Jac`: vectors, the dense linear solve, damped Newton, forward-difference Jacobians, and the failures they return. It knows nothing about ODEs and has no effect module (see "Effects"). |
| `vstiff`, the solver | [src/dune](../src/dune) | `Vstiff` | `Ode`, `Fail`, `Clock`, `Instrument`, `Stage`, `Check`, `Bdf1`, `Bdf2`, `Halving`, `Stepper`, `Adaptive`. It depends on `numerics`. |

The dependency runs one way. Dune rejects a cycle between libraries, so nothing in the kernel can name the solver, and `Newton` cannot know what an `Ode.problem` is. A solver module names what it takes from the kernel at its top, `module Vec = Numerics.Vec` (`Stage` also takes `Jac` and `Newton`).

**The main module.** [src/vstiff.ml](../src/vstiff.ml) is written by hand, one alias per public module (`module Ode = Ode`, `module Fail = Fail`, and so on), and [src/vstiff.mli](../src/vstiff.mli) documents each. They list `Ode`, `Fail`, `Clock`, `Instrument`, `Bdf1`, `Bdf2`, `Halving`, `Stepper` and `Adaptive`: that list is the public API. `Stage` and `Check` are not on it, so `Vstiff.Stage` and `Vstiff.Check` are unbound outside the library, and so is `Vstiff.Vec`, since the kernel's modules are reached as `Numerics.Vec`, not through the solver:

```text
Error: Unbound module Vstiff.Stage
```

**Private modules, and why they are not enough.** `src/dune` also says `(private_modules stage check)`, which keeps the compiled interfaces of the two modules out of the installed form of the library; the libraries have no `public_name`, so nothing is installed from this project today. On its own the field does not hide a module inside the workspace that builds the library, which is where the tests and every probe live: with `src/vstiff.ml` and `src/vstiff.mli` deleted, so that dune writes a main module of its own, `Vstiff.Stage.solve` compiles. What hides `Stage` and `Check` is the hand-written main module, which does not list them; `private_modules` adds the guarantee for an installed copy.

**`Fail` is re-exported.** The failures belong to the kernel, because `Newton.solve` returns them. [src/fail.ml](../src/fail.ml) is `include Numerics.Fail` and [src/fail.mli](../src/fail.mli) is `include module type of struct include Numerics.Fail end`, so `Vstiff.Fail` says `type t = Numerics.Fail.t = Diverged | StepRejected of int | Nan`: a failure from the kernel is a `Vstiff.Fail.t` without conversion, `Diverged` can be written `Vstiff.Fail.Diverged` or `Numerics.Fail.Diverged`, and a program that only uses the solver never names the kernel to read a failure ([ocaml.md](ocaml.md) explains the spelling).

**Vectors are `float array` in public signatures.** `Ode.rhs` is `float -> float array -> float array`, and `problem.y0`, `point.y`, `Method.step`, `Embedded.step_with_error`, `Controller.acceptable ~y ~err`, `Adaptive.solution.y` and the result of `Stepper.fixed` say `float array` too, not `Vec.t`. `Vec.t` is `float array`, so nothing converts, and a program that implements `Ode.Method` or `Ode.Controller` needs the solver library alone: the controller example under "Adding things" uses nothing from `Numerics`.

## The contracts in `Ode`

[src/ode.mli](../src/ode.mli) has no implementation: it declares the vocabulary and three module types ([ocaml.md](ocaml.md) explains the terms).

| Item | What it is |
|---|---|
| `Ode.rhs` | `float -> float array -> float array`: `rhs t y` is $f(t, y)$ as a fresh vector on every call, and does not modify `y` |
| `Ode.problem` | `{ rhs; t0; t_end; y0 }` |
| `Ode.point` | `{ t; y }`, a point on a solution |
| `Ode.rejection` | `Too_large` (the error estimate was over tolerance), `Solver of Fail.t` (the method could not take the step) or `Too_small` (the step is below the resolution of $t$, and might not move $t$ at all, so the method was not called) |

| Module type | It provides | Implemented by | Used by |
|---|---|---|---|
| `Ode.Method` | `type history`, `start`, `step rhs h history at` returning `(y, history)` | `Bdf1` (history is `unit`), `Bdf2` through `Embedded` | `Stepper.fixed` |
| `Ode.Embedded` | a `Method` plus `step_with_error`, returning `(y, err, history)`, where `err` estimates the local error of the step | `Bdf2` | `Adaptive.integrate` |
| `Ode.Controller` | `type t`, `type stats`, `init`, `proposal`, `acceptable`, `accepted`, `rejected`, `stats` | `Halving` | `Adaptive.integrate` |

A method keeps no state of its own: whatever it needs from the past it returns as `history` and the driver hands it back (`Bdf2`'s history is `Start`, or `After { h_prev; y_prev }`). A controller is the same: every transition returns a new `t`. `Controller.stats` is a type the controller chooses, which is why the result type of `Adaptive.integrate` depends on it (`Halving.stats` holds the accepted and rejected counts). `Fail.t` is `Diverged` (Newton found no root), `StepRejected n` (a controller gave up) or `Nan` (a value that had to be finite was not).

## Module dependency graph

The project modules each file needs. Within one library the source is `dune describe workspace --with-deps` (the `module_deps` of each module, interface and implementation together). What a solver file takes from the kernel is the `Numerics` that `ocamldep -modules src/*.ml src/*.mli src/numerics/*.ml src/numerics/*.mli` prints for it (it lists the standard library modules too), and shows in the `module Vec = Numerics.Vec` lines at the top of the file; the tables write it `Numerics.Vec`. A module can only use modules compiled before it, and cycles are errors. Each table is in a valid compile order, the kernel first.

The kernel, library `numerics`:

| Module | What it is | Uses | Explained in |
|---|---|---|---|
| `Fail` | named failures, `let*` and `let+` | | [ocaml.md](ocaml.md) |
| `Vec` | float-array vectors, never mutated | | [ocaml.md](ocaml.md) (arrays) |
| `Linalg` | Gaussian elimination with partial pivoting | `Vec` | [numerics/02-newton.md](numerics/02-newton.md) |
| `Newton` | damped Newton for $G(x) = 0$ | `Fail` `Linalg` `Vec` | [numerics/02-newton.md](numerics/02-newton.md) |
| `Jac` | forward-difference Jacobian | `Linalg` `Vec` | [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md) |

The solver, library `vstiff`:

| Module | What it is | Uses | Explained in |
|---|---|---|---|
| `Clock` | the resolution of time as floats: $16\thinspace\varepsilon\thinspace\lvert t \rvert$ | | [numerics/05-step-control.md](numerics/05-step-control.md) |
| `Fail` | the kernel's `Fail`, re-exported | `Numerics.Fail` | this page |
| `Ode` | the contracts above (interface only) | `Fail` | this page |
| `Check` (internal) | argument checks for the drivers | `Clock` `Ode` | this page |
| `Instrument` | counts calls of an `rhs` | `Ode` | this page |
| `Stage` (internal) | the stage equation, solved by Newton | `Fail` `Ode` `Numerics.Jac` `Numerics.Newton` `Numerics.Vec` | [numerics/04-bdf.md](numerics/04-bdf.md) |
| `Bdf1` | backward Euler, a `Method` | `Ode` `Stage` | [numerics/04-bdf.md](numerics/04-bdf.md) |
| `Bdf2` | variable-step BDF2, an `Embedded` | `Bdf1` `Fail` `Ode` `Stage` `Numerics.Vec` | [numerics/04-bdf.md](numerics/04-bdf.md) |
| `Halving` | the `Controller`: halve, double | `Clock` `Fail` `Ode` | [numerics/05-step-control.md](numerics/05-step-control.md) |
| `Stepper` | the fixed-step driver | `Check` `Fail` `Ode` `Numerics.Vec` | this page |
| `Adaptive` | the adaptive driver | `Check` `Clock` `Fail` `Ode` `Numerics.Vec` | [numerics/05-step-control.md](numerics/05-step-control.md) |
| `Vstiff` | the main module: the public API | `Adaptive` `Bdf1` `Bdf2` `Clock` `Fail` `Halving` `Instrument` `Ode` `Stepper` | this page |

The tests use `Problems` (`Ode`), `Refs` and `Report` (nothing), `Guard` (`Fail`, `Instrument`), and the three programs: `Corpus` uses `Adaptive`, `Bdf1`, `Bdf2`, `Fail`, `Guard`, `Halving`, `Instrument`, `Jac`, `Linalg`, `Newton`, `Ode`, `Problems`, `Refs`, `Report`, `Stepper`, `Vec`; `Soak` the same minus `Instrument`, `Jac`, `Linalg` and `Newton`; `Props` uses `Bdf2`, `Dd`, `Gen`, `Jac`, `Linalg`, `Newton`, `Prop`, `Report` and `Vec`. `Gen`, `Prop` and `Dd` use no library module. `Problems` and `Guard` name only the solver library; `corpus.ml`, `soak.ml` and `props.ml` link both and start with `open Vstiff` and `open Numerics`. Nothing in the library uses a test module.

What the graph shows: `Stepper` depends only on `Ode`, `Check`, `Fail` and the kernel's `Vec`, and `Adaptive` also on `Clock`, so they reach no method, no controller and no Newton code except through the modules they are given. `Clock` depends on nothing; `Halving`, `Adaptive` and `Check` use it, so the floor of the controller, the remainder rule and the rejection of too-short steps in the adaptive driver, and the shortest `dt` that `Stepper.fixed` accepts are one rule for the shortest step. `Bdf2` uses `Bdf1` for its first step and for its error estimate. Only `Stage` connects the methods to the kernel's `Newton` and `Jac` (`Bdf2`, `Stepper` and `Adaptive` use only its `Vec`), and the kernel uses nothing of the solver.

## The life of one adaptive step

`Adaptive.integrate (module M) (module C) ?dt0 ?dt_max ?max_rejects ~tol problem` first computes the span $t_{\mathrm{end}} - t_0$, takes `dt_max` (default $(t_{\mathrm{end}} - t_0) / 10$) and `dt0` (default $10^{-6} (t_{\mathrm{end}} - t_0)$, capped at `dt_max`; a default that underflows to 0 on a denormal span becomes the span), calls `Check.adaptive` (it raises `Invalid_argument` if `t_end < t0`, if the span is not empty and `dt0` is not positive, or if `tol` is not positive, `nan` included), evaluates `rhs t0 y0` once and passes it to `Check.output` (it raises if the vector is not as long as `y0`), returns `Error Nan` if `y0` or that vector is not finite, and starts the loop `go (C.init ~tol ~dt0 ~dt_max ~max_rejects) M.start { t = t0; y = y0 }`. One trip through `go c history at`:

1. **Done?** If `at.t >= t_end`, return `Ok { t; y; stats = C.stats c }`.
2. **Choose the step.** `dt = C.proposal c` and `remaining = t_end - at.t`. If `dt >= remaining`, or `remaining <= Clock.resolution at.t`, this is the last step: `h = remaining` and the new time is `t_end`. Otherwise the new time is `t_next = at.t + dt` and the step is snapped to the floats, `h = t_next - at.t`, so that the state advances by exactly the step the clock took.
3. **Too short?** If $h \le 0$, or the step is not the last and $h$ is below `Clock.resolution at.t`, the step is below the resolution of $t$ and might not move $t$ at all: the method is not called, and the attempt goes to item 7 with `Ode.Too_small`.
4. **Attempt it.** `M.step_with_error p.rhs h history at`. For `Bdf2`:
   - `backward_euler` calls `Bdf1.step`, which calls `Stage.solve rhs { Stage.t = at.t +. h; gamma = h; psi = at.y } at.y`. `Stage.solve` forms the residual $G(x) = x - \psi - \gamma\thinspace f(t_{n+1}, x)$ and the matrix $I - \gamma J$ with $J$ from `Jac.forward`, and `Newton.solve` iterates: `Linalg.solve` gives the step, a line search shortens it if the residual would grow.
   - With no history (`Start`, the first step and its retries) the result is the backward Euler state, and the error estimate is half its gap to the explicit Euler state $y_n + h\thinspace f(t_n, y_n)$.
   - With history `After { h_prev; y_prev }`, `bdf2` takes $\omega = h / h_{\mathrm{prev}}$ (`omega` in the code), the weights from `Bdf2.coeffs`, $\psi = a_1 y_n + a_0 y_{n-1}$, $\gamma = \beta h$, and a starting guess extrapolated along the last two points, and calls `Stage.solve` again. The result is the BDF2 state, and the estimate is its gap to the backward Euler state from the same point.
   - Either way the new history is `After { h_prev = h; y_prev = at.y }`.
5. **Judge it.** `C.acceptable c ~y ~err`; `Halving` accepts if $\max_i \lvert \mathrm{err}_i \rvert / (1 + \lvert y_i \rvert)$ is at most $\mathrm{tol}$.
6. **Accepted:** `go (C.accepted c) next { t = t_next; y }`, where `t_next` is `t_end` on the last step. `Halving` counts it, resets the failure count, and doubles `dt` (up to `dt_max`) at every third accept in a row.
7. **Rejected** (`Ode.Too_small` from item 3, `acceptable` said no, giving `Ode.Too_large`, or `step_with_error` returned `Error e`, giving `Ode.Solver e`): `C.rejected c reason ~at:at.t ~h`. `Halving` sets `dt = h / 2` and counts the failure, or gives up with `Error (StepRejected n)` when more than `max_rejects` rejections came in a row or the halved step is below `Clock.resolution at`, which is $16\thinspace\varepsilon\thinspace\lvert t \rvert$ ($\varepsilon$ is `Float.epsilon`, about $2.2 \times 10^{-16}$), or is 0 (the floor is 0 at $t = 0$). Otherwise `go` retries from the **same** `at` and `history`.

`go` calls itself last in every branch, so the number of steps is not limited by the stack. The calls, indented by who calls whom:

```text
Adaptive.integrate        checks the arguments, then loops with go c history at
  go                      ends with Ok when at.t >= t_end
    C.proposal c          the step to try; h is that snapped to the floats, or what is left
                          to land on t_end (a rest of at most Clock.resolution at.t)
                          h <= 0, or h below Clock.resolution unless last:
                            C.rejected c Too_small ~at ~h, and no call of the method
    M.step_with_error     Bdf2.step_with_error
      backward_euler      Bdf1.step -> Stage.solve -> Newton.solve (Jac.forward, Linalg.solve)
      bdf2                Stage.solve -> Newton.solve      (every step after the first)
    C.acceptable          yes: go (C.accepted c) next { t = t_next; y }
                          no, or Error e: C.rejected c reason ~at ~h
                            Ok c'   -> go c' history at    (the same point and history)
                            Error _ -> the run ends with that Error
```

To watch the loop, wrap `Halving` so that it prints each call of the contract (a probe, [testing.md](testing.md)). `include Halving` copies everything, and the four definitions after it replace the originals:

```ocaml
open Vstiff

module Traced = struct
  include Halving

  let proposal c =
    let dt = Halving.proposal c in
    Printf.printf "proposal %g\n" dt;
    dt

  let acceptable c ~y ~err =
    let ok = Halving.acceptable c ~y ~err in
    Printf.printf "acceptable: %b\n" ok;
    ok

  let accepted c =
    print_endline "accepted";
    Halving.accepted c

  let rejected c reason ~at ~h =
    Printf.printf "rejected: at t = %g, h = %g\n" at h;
    Halving.rejected c reason ~at ~h
end

let problem = { Ode.rhs = (fun _t y -> Array.map Float.neg y); t0 = 0.; t_end = 0.2; y0 = [| 1. |] }

let () =
  match Adaptive.integrate (module Bdf2) (module Traced) ~dt0:0.1 ~dt_max:0.1 ~tol:1e-3 problem with
  | Ok s -> Printf.printf "done at t = %g\n" s.t
  | Error e -> print_endline (Fail.to_string e)
```

For $y' = -y$ the first step has no history, and its scaled estimate works out to $h^2 / (2(2 + h))$: $2.4 \times 10^{-3}$ for $h = 0.1$, above $\mathrm{tol}$, so the first attempt ends in `rejected` and the retry has $h = 0.05$. After the third `accepted` in a row the proposal doubles, and the driver cuts the last step to land on $t_{\mathrm{end}}$, so the proposal printed for it is longer than the step taken.

`Stepper.fixed (module M) ~dt problem` is the simpler loop. `Check.fixed` raises unless `dt > 0` and `t_end >= t0`, and when the span is not empty it also raises if `dt` is below `Clock.resolution` of the larger of $\lvert t_0 \rvert$ and $\lvert t_{\mathrm{end}} \rvert$: a step that short cannot be told apart on the clock, and the same check keeps the number of steps at most $1 / (8 \varepsilon) = 2^{49}$, so the conversion to an integer is safe. The driver then evaluates `rhs t0 y0` once, passes it to `Check.output` (it raises unless the vector is as long as `y0`) and returns `Error Nan` if `y0` or that vector is not finite. An empty span returns `y0`; otherwise it takes $n = \max\bigl(1, \mathrm{round}((t_{\mathrm{end}} - t_0) / \mathtt{dt})\bigr)$ steps with `M.step`, stops at the first `Error`, and returns the final state. `dt` is only a target, and the method never sees it. Step $k$ ends at the time $t_0 + k h$ with $h = (t_{\mathrm{end}} - t_0) / n$, and the last one at $t_{\mathrm{end}}$ itself; the step handed to the method is the difference of its two end times, as in the snapping of the adaptive driver, so no running sum can drift off $t_{\mathrm{end}}$ and the state advances by exactly what the clock does.

## How errors flow

Numerical failures are `Fail.t` values in a `result`; they enter at the bottom and travel up through `Result.bind`, `Result.map`, `let*` and `let+`.

| Where | What happened | What comes out |
|---|---|---|
| `Linalg.solve` | an exactly zero pivot | `None` |
| `Newton.solve` | `None`, a non-finite step, no acceptable line-search factor, or `max_iter` reached | `Error Diverged` |
| `Newton.solve` | the iterate or its residual is not finite at the start of an iteration, or a converged step $x + \Delta x$ overflows | `Error Nan` |
| `Stage.solve`, `Bdf1.step`, `Bdf2.step_with_error` | pass it through | the same `Error` |
| `Stepper.fixed` | any step fails | the first `Error`; no further step is taken |
| `Adaptive.integrate` | a step fails, its estimate is not acceptable, or it is below the resolution of $t$ ($h \le 0$, or a non-final $h$ below `Clock.resolution`; the method is not called) | a rejection handed to `C.rejected` |
| `Halving.rejected` | too many in a row, the step floor, or a halved step of 0 | `Error (StepRejected n)`, which `Adaptive.integrate` returns |
| `Stepper.fixed`, `Adaptive.integrate` | `y0` or `rhs t0 y0` is not finite | `Error Nan` |
| `Check` | arguments that make no sense: `dt`, `dt0`, `dt_max` or `tol` not positive, `dt` below the resolution of $t$, `t_end < t0`, an `rhs t0 y0` of the wrong length | raises `Invalid_argument` |

The same Newton failure means different things at different levels. In a fixed-step run it ends the run; in `Adaptive` a smaller step usually fixes it, so it is a routine rejection. The `n` in `StepRejected n` is the length of the final run of rejections: `max_rejects + 1` when that limit stopped the run, less when the step floor or a step of 0 did. Five corpus lines end with `StepRejected`, none of them on the limit: `StepRejected 1` for a blow-up (the attempt before that rejection was accepted), `StepRejected 46` for a right-hand side that turns into `nan`, `StepRejected 1` twice for steps too short to move $t$, where the first `Too_small` rejection halves a step of 0, and `StepRejected 1055` for a right-hand side that is `nan` after a start at $t = 0$, where the floor is 0 and only the underflow of the halved step to 0 ends the run ([test/corpus.expected](../test/corpus.expected)). [numerics/05-step-control.md](numerics/05-step-control.md) (section 9) derives the counts and has a probe that ends on the limit.

## Effects

OCaml does not record effects in types (a function's type does not say whether it raises, mutates or prints), so the discipline is by module ([AGENTS.md](../AGENTS.md), H6): a function outside these modules neither raises on purpose, nor mutates, nor prints, and a new effect goes into one of them or into a new module that exists for it.

| Module | Effect | Why there |
|---|---|---|
| `Check` (library, internal) | raises `Invalid_argument` | a bad argument is a programming error, not a numerical failure; every deliberate raise is in one place, called at the start of a driver, before it takes a step |
| `Instrument` (library) | one `ref`: a call counter | counting needs state; wrapping the `rhs` a method receives keeps the methods pure. Nothing in the library calls it: `Guard.budget`, the `too_small` and `adaptive_canary` cases of the corpus and probes do |
| `Guard` (tests) | raises `Exhausted`, catches it and `Invalid_argument`, and catches whatever a property case raises | turns a runaway loop or a bad argument into a printed line (`run`), an exhausted budget into `None` (`bounded`), so that the soak test can count the round as not passed, and a raising property case into a failing one (`verdict`) |
| `Report` (tests) | prints | the only module that writes to standard output |
| `Gen` (tests) | advances the `Random.State.t` it is given | drawing a case moves the generator; `Prop.check` makes a fresh state for each property, so the cases do not depend on the order the properties run in |
| `props.ml` main (tests) | reads `VSTIFF_PROP_SEED` and `VSTIFF_PROP_SCALE` | once, before the suite runs: a nightly run sets another seed or more cases |

Everything else is pure: methods, controllers, the whole kernel (`Newton`, `Jac`, `Linalg`, `Vec`, `Fail`: `Fail.to_string` formats a string and prints nothing), and the drivers, which raise on purpose only by calling `Check`. Why:

- **Retrying is safe.** A rejected step leaves `at` and `history` untouched, and a failed attempt leaves nothing to clean up.
- **Sharing is safe.** The history keeps `y_prev` while `at.y` is the current state, and `Newton.solve` returns its starting vector itself if the residual there is exactly zero. This works only because nothing ever writes into an array: every `Vec` operation allocates its result.
- **Results are reproducible**, which the soak test relies on, and a test case is a pure function from `()` to text.
- **Effect types stay easy to annotate**: if effects are ever tracked in types, only these modules carry them.

## Invariants

A change that breaks one is a bug even if every corpus line still passes.

1. **Arrays are never mutated after creation** (above), and `rhs` returns a fresh one on every call. `Jac.forward` keeps `f y` while it evaluates `f` at the perturbed points, so a right-hand side that wrote every result into one buffer would make every difference zero.
2. **Jacobian orientation.** `J.(i).(j)` is $\partial f_i / \partial y_j$: the row is the output, the column the input. `Jac.forward` builds the columns first and then reads rows from them, `Stage.solve` forms $I - \gamma J$ entry by entry with the same indices, and `Linalg.solve a b` treats `a.(i)` as equation $i$. A transposed matrix still type-checks, and the canary's diagonal Jacobian is its own transpose. The root Newton converges to does not depend on the Jacobian, only the speed does, so a wrong Jacobian costs speed, not correctness: Newton converges slowly or not at all, `Adaptive` counts the failure as a rejection, and the run crawls ([testing.md](testing.md)).
3. **Both drivers land exactly on $t_{\mathrm{end}}$.** `Adaptive` cuts the last step to `t_end - at.t` and assigns the new time `t_end` rather than computing `at.t + h`, which could miss it by rounding. `Stepper.fixed` ends step $k$ at $t_0 + k h$ and the last one at $t_{\mathrm{end}}$ itself, with no running sum of steps: nine additions of $1/9$ give $1.0000000000000002$, and a right-hand side that is undefined past $t_{\mathrm{end}}$ (corpus line 29) returned `nan` at the last stage.
4. **The state advances by exactly the clock's step.** A step that is not the last is snapped: the new time is $t_{n+1} = t_n + \mathtt{dt}$ and $h = t_{n+1} - t_n$, the difference of the two clock readings (computed without rounding when they are within a factor of 2 of each other, as they are for any step much shorter than $\lvert t_n \rvert$), and the method moves the state by $h$. An unsnapped step $h = \mathtt{dt}$ would move the state by $\mathtt{dt}$ and the clock by $t_{n+1} - t_n$, which differs from $\mathtt{dt}$ by up to half an ulp of $t_{n+1}$ (an ulp of $t_n$ at most): at $t = 10^{15}$, where the floats are 0.125 apart, a step of 3.7 moves the clock by 3.75, and corpus line 27, which runs steps of that length, ends at $y = 98.7$ unsnapped and at $y = 100$ snapped. The price is that a snapped step may exceed `dt_max` by that much. A non-final step below `Clock.resolution` of $t$ never reaches a method (it is rejected as `Too_small`, as is every step with $h \le 0$), and above the resolution snapping changes a step by at most $1/32$ of its length (half an ulp against at least 16 ulps). `Stepper.fixed` keeps the invariant with its grid: the step handed to the method is the difference of the two grid times, not $h$.
5. **Determinism within one build** ([AGENTS.md](../AGENTS.md), H1). No randomness, hidden state or parallelism: the same arguments give identical results in the same build ([test/soak.ml](../test/soak.ml) checks it). Results may differ in the last bits between machines (see "Performance"), so the expected files print few digits.
6. **Vector lengths agree.** `Vec.add`, `sub`, `axpy` and `dot` index by the length of their first vector argument: a shorter second argument raises `Invalid_argument` (index out of bounds), a longer one is ignored past that length. `Halving.acceptable` uses `Array.map2`, which raises if `y` and `err` differ in length. These are caller bugs, not numerical failures. The drivers check one length, that of `rhs t0 y0` against `y0`, with `Invalid_argument` naming the driver; a right-hand side whose length changes later is not checked.
7. **History convention.** `After { h_prev; y_prev }` pairs the step just taken with the state it started from. A rejection leaves the history alone, so $\omega = h / h_{\mathrm{prev}}$ shrinks when $h$ does.
8. **Step ratios stay below the stability limit.** Variable-step BDF2 is zero-stable (earlier errors stay bounded) for $\omega \lt 1 + \sqrt{2}$, about 2.414 ([numerics/04-bdf.md](numerics/04-bdf.md)). `Bdf2` does not check it; the controller must. `Halving` proposes at most twice the last accepted step, so $\omega$ stays at most 2 for the proposals. The driver snaps them to the floats, which changes a step by at most $1/32$, so the $\omega$ it takes is at most $2(1 + 1/32) / (1 - 1/32) = 66/31$, about 2.13, inside the limit ([numerics/05-step-control.md](numerics/05-step-control.md), section 9). A new controller must keep it below the limit.
9. **The solver sees $f$ as a black box (H7).** Jacobians are always forward differences inside `Stage.solve`, and corpus problems never supply analytic ones.

## Conventions

The conventions are in [AGENTS.md](../AGENTS.md): the hard rules H5 (failures are values, invalid arguments raise), H6 (effects live in named modules, above), H9 (standard library only) and H10 (two libraries, one public API), and its "Design defaults", which cover the first two items below. This is how they show in the code.

- **Contracts and implementations.** Contracts are module types in `Ode`; methods and controllers are modules that implement them, passed to the drivers as modular explicits. That is the default because the result type can name the module's own types, and nothing has to be packed or applied ([ocaml.md](ocaml.md)). A functor, a first-class module in a data structure (when the method comes from a list at run time, say), a GADT or an effect handler is the better tool where it gives the clearer design. Types are shared between signatures (`with type`) only where a caller must see them: `Halving` does so for `stats`, so that the tests can read the counts.
- **Names.** Module types are CamelCase (`Method`, `Controller`), values snake_case. `dt` is a requested step size (`~dt`, `dt0`, `dt_max`, `Controller.proposal`), `h` a step taken or attempted. Labels go where two arguments of one type could be swapped (`~at ~h`, `~y ~err`), or name a bare literal at a call site (`~dt:2e-6`); other arguments are positional.
- **Failures and mistakes (H5).** A numerical failure is a `Fail.t` in a `result`; an invalid argument raises `Invalid_argument` from `Check`.
- **Two libraries, one public API (H10).** The kernel in `src/numerics/` knows nothing about ODEs, the solver offers only what `src/vstiff.ml` and `src/vstiff.mli` list, and a vector is a `float array` in public signatures (above).
- **Standard library only (H9).** `src/numerics/dune` lists no libraries, `src/dune` only `numerics` and `test/dune` `vstiff` and `numerics`.

## Adding things

The checklists are in [AGENTS.md](../AGENTS.md), under "Tests": adding a corpus case, a soak case, a method or a controller. The table and the notes after it say what to edit and why.

| You want to | Edit | Also |
|---|---|---|
| add a method | a module implementing `Ode.Method` (see below) | `include Ode.Method` in its `.mli`; list it in `src/vstiff.ml` and `src/vstiff.mli`; a corpus case; [numerics/04-bdf.md](numerics/04-bdf.md) if it is a BDF-type method |
| add a controller | a module implementing `Ode.Controller` (see below) | list it in `src/vstiff.ml` and `src/vstiff.mli`; pass it to `Adaptive.integrate`; a corpus case |
| make an internal module public (H10) | [src/vstiff.ml](../src/vstiff.ml) and [src/vstiff.mli](../src/vstiff.mli) | one `module X = X` line in each, with a doc comment in the `.mli`; take the module out of `private_modules` in `src/dune` if it is listed there |
| add a numerical building block | a module in `src/numerics/`, which must not name `Ode` | an `.mli`; a corpus case (the `newton` and `jacobian` groups call the kernel directly) |
| add a corpus problem | a module in [test/problems.ml](../test/problems.ml) with `rhs`, `y0`, `problem` and, if a closed form exists, `exact`; no Jacobian (H7) | a case in [test/corpus.ml](../test/corpus.ml) by the steps in [testing.md](testing.md) |
| add a kind of failure | `Fail.t` and `Fail.to_string` in [src/numerics/fail.ml](../src/numerics/fail.ml) and its `.mli` | every `match` on `Fail.t` without a catch-all stops compiling until it handles the new case, which is the point; `Vstiff.Fail` follows by itself, being a re-export |
| change Newton's tuning constants | the constants at the top of [src/numerics/newton.ml](../src/numerics/newton.ml) | [testing.md](testing.md) lists what notices each |
| change how the Jacobian is perturbed | `Jac.step` in [src/numerics/jac.ml](../src/numerics/jac.ml) | [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md) |

**A method.** Write `src/your_method.ml` and an `.mli` that says `include Ode.Method`, as [src/bdf1.mli](../src/bdf1.mli) does; dune picks up a new module in `src/` without a change to `src/dune`, and until the module is listed in `src/vstiff.ml` and `src/vstiff.mli` (`module Your_method = Your_method` and a documented line) it is internal like `Stage`: the tests and other programs cannot name it. An implicit method calls `Stage.solve rhs { Stage.t = at.t +. h; gamma; psi } guess` with its own $\psi$ and $\gamma$, as `Bdf1` and `Bdf2` do; an explicit one needs no Newton. [numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md) ("Running the snippets") implements explicit Euler as an `Ode.Method` in a few lines and runs it with `Stepper.fixed`. To use a method adaptively it must also provide `step_with_error`, which returns the estimated local error vector in the middle: `Ode.Embedded`. A method written outside the library, in a probe or another project, implements `Ode.Method` in the same way but cannot see `Stage`: it builds its stage equation with `Newton.solve` and `Jac.forward` from the kernel (capstone E1 in [exercises.md](exercises.md) does it for the trapezoidal rule).

**A controller.** Write `src/your_controller.ml` and an `.mli` like [src/halving.mli](../src/halving.mli): declare `stats` as a concrete type, then `include Ode.Controller with type stats := stats`, so that callers can read it, and list it in `src/vstiff.ml` and `src/vstiff.mli` to make it public. The example below is written outside the library and needs the solver library alone; it keeps everything in one module, whose signature `Ode.Controller with type stats = int` shows callers the type of `stats`. It accepts a step if the largest error entry is at most `tol`, never changes the step, and gives up at the first rejection, which a controller must do in the end: `Adaptive.integrate` terminates only because every run of rejections ends ([numerics/05-step-control.md](numerics/05-step-control.md), section 9). Its `stats` is an `int`, so `Adaptive.integrate (module Bdf2) (module Give_up_at_once) ~tol problem` has type `(int Adaptive.solution, Fail.t) result`:

```ocaml
open Vstiff

module Give_up_at_once : Ode.Controller with type stats = int = struct
  type t = { tol : float; dt : float; accepted : int }
  type stats = int

  let init ~tol ~dt0 ~dt_max:_ ~max_rejects:_ = { tol; dt = dt0; accepted = 0 }
  let proposal c = c.dt
  let acceptable c ~y:_ ~err = Array.for_all (fun e -> Float.abs e <= c.tol) err
  let accepted c = { c with accepted = c.accepted + 1 }
  let rejected _ _ ~at:_ ~h:_ = Error (Fail.StepRejected 1)
  let stats c = c.accepted
end

let () =
  let problem = { Ode.rhs = (fun _t y -> Array.map Float.neg y); t0 = 0.; t_end = 1.; y0 = [| 1. |] } in
  match Adaptive.integrate (module Bdf2) (module Give_up_at_once) ~dt0:0.01 ~tol:1e-3 problem with
  | Ok s -> Printf.printf "accepted %d steps, t = %g\n" s.stats s.t
  | Error e -> print_endline (Fail.to_string e)
```

Run it as a probe ([testing.md](testing.md)). With `~dt0:0.1` the first estimate, $h^2 / (2(1 + h))$ for this problem, already exceeds $\mathrm{tol}$, and `rejected` ends the run with `Error (StepRejected 1)`.

**A corpus problem.** Add a module to [test/problems.ml](../test/problems.ml) with `rhs`, `y0`, `problem = { Ode.rhs; t0; t_end; y0 }` and, if a closed form exists, `exact`, and a doc comment saying what it is and why the corpus needs it. Give it no Jacobian (H7): the solver must work for any black-box `rhs`. With no closed form, the reference answer is computed outside vstiff and goes into a new module, since [test/refs.ml](../test/refs.ml) is never edited (H4; [testing.md](testing.md), "Reference values"). A problem changes no output until a case in [test/corpus.ml](../test/corpus.ml) uses it ([testing.md](testing.md), "Adding a case, step by step").

## Performance

Where the time goes: a step costs its Newton iterations, each with $n + 1$ calls of `rhs` for the Jacobian and one dense linear solve, and an adaptive attempt after the first makes two such stage solves. The bullets take these parts in turn.

- **The canary dominates the test time.** The backward Euler run on the canary (the corpus problem of three independent decays, [numerics/06-the-corpus.md](numerics/06-the-corpus.md)) takes $(t_{\mathrm{end}} - t_0) / \mathtt{dt} = 1 / (2 \times 10^{-6}) = 500{,}000$ fixed steps, each with a Newton solve; the corpus runs it once and the soak test ten times. Measure with `time ./_build/default/test/soak.exe`.
- **The Jacobian is rebuilt in every Newton iteration**, and `Jac.forward` calls `rhs` once at $y$ and once per perturbed component: $n + 1$ calls for $n$ components, plus one call for each trial point of the line search. Production codes reuse a Jacobian across iterations and steps. The probe counts the calls of `Jac.forward` and of one backward Euler step (`Stage` is internal, so the step goes through `Bdf1.step`, which calls `Stage.solve`); the first line is derived, 4 for the canary's three components, and the second depends on how many iterations Newton needs:

```ocaml
open Vstiff
open Numerics

let () =
  let rhs, calls = Instrument.count Problems.Canary.rhs in
  ignore (Jac.forward (rhs 0.) Problems.Canary.y0);
  Printf.printf "Jac.forward, 3 components: %d calls of rhs\n" (calls ());
  let rhs, calls = Instrument.count Problems.Canary.rhs in
  ignore (Bdf1.step rhs 1e-3 Bdf1.start { Ode.t = 0.; y = Problems.Canary.y0 });
  Printf.printf "one backward Euler step: %d calls of rhs\n" (calls ())
```

- **An adaptive attempt after the first step costs two stage solves**, one for BDF2 and one for backward Euler; the first step needs only the backward Euler solve and one `rhs` call for the explicit Euler state. An attempt rejected as too large costs as much as an accepted one; one rejected as `Too_small` costs nothing, because the method is not called.
- **Dense, allocating linear algebra.** `Linalg.solve` recurses on the trailing submatrix: $O(n^3)$ arithmetic and, since every level builds fresh arrays, $O(n^3)$ allocation. It is meant for small systems.
- **Fused multiply-add.** On arm64, `ocamlopt` fuses `a +. b *. c` and `a -. b *. c` into one instruction with a single rounding when the product is a direct operand; binding the product with `let` first rounds twice, and can give a different result. So moving a product into or out of a sum can change the last bits of a result, results may differ between machines, and the expected files print few digits. That is allowed: no check may depend on the last bits ([AGENTS.md](../AGENTS.md), H1; [testing.md](testing.md) says what to do when a line moves), and `Float.fma` is welcome wherever it makes a result more accurate. Run `dune runtest` after any arithmetic edit all the same. [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md) (section 8) has a program that shows it.

## Limits of the design

Stated as limits and as starting points for contributions ([exercises.md](exercises.md); terms in [glossary.md](glossary.md)):

- Newton rebuilds the Jacobian at every iteration.
- `Halving` is halve/double. Production controllers scale $h$ by a safety factor times $(\mathrm{tol} / \mathrm{err})^{1/(p+1)}$, $p$ being the order of the method whose error is estimated; that needs the error estimate when a step is accepted, which `Ode.Controller.accepted` does not receive, so it needs a contract change ([numerics/05-step-control.md](numerics/05-step-control.md), section 11). The BDF1/BDF2 gap used as the estimate is conservative for BDF2.
- One mixed absolute/relative weight, $1 / (1 + \lvert y_i \rvert)$, controls the error: there is no separate $\mathrm{rtol}$ and $\mathrm{atol}$ (relative and absolute tolerance), so tiny components such as Robertson's $y_2$ are controlled only loosely.
- `Halving` counts every rejection together, whatever its reason; a controller could count Newton failures separately, since `Ode.rejection` says which it was.
- `Instrument` counts `rhs` calls only; Newton iterations per step are not visible from outside.
- No dense output (values between the steps): the drivers return only the final state.
