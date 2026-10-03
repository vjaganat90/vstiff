# 1. ODEs and stiffness

Chapter 1 of 6 in the numerical-methods track. Next: [2. Newton's method](02-newton.md). Index and reading order: [docs/README.md](../README.md). Terms and symbols: [glossary](../glossary.md).

**Summary.** An initial value problem is a rule for how a state changes plus a starting point. A computer follows its solution one small step at a time, and explicit Euler is the simplest way to do that. This chapter measures the error of such a method (local versus global error, the order), shows what a too-long step does (instability) and introduces stiff problems, where explicit methods are forced to take absurdly short steps. Backward Euler escapes that limit but charges a price: a nonlinear equation to solve at every step. That price is why the library needs Newton's method, a Jacobian and a linear solver.

You need first-year calculus (derivatives, Taylor series). The little linear algebra that appears, eigenvalues and the Jacobian, is explained where it first shows up. OCaml is needed only for the optional snippets: "Running the snippets" near the end says how to run them, and [docs/ocaml.md](../ocaml.md) teaches the language.

Notation (the glossary has the full list): `t` is time, `y` the state, `f` the right-hand side, so `y' = f(t, y)`; the code calls `f` by the name `rhs`. `h` is a step taken, `t_{n+1} = t_n + h`, and `y_n` is our approximation of the exact `y(t_n)`. `dt` is a step *requested* by the caller, as in `Stepper.fixed ~dt`.

## 1. The problem: an initial value problem

An ordinary differential equation (ODE) says how a quantity changes in time, as a function of the time and of its current value. "Ordinary" means there is one independent variable (here time); partial differential equations have several. An initial value problem (IVP) adds a starting point.

Take a population `y(t)` measured as a fraction of what its environment can sustain, so `y = 1` means "full". It grows in proportion to its size `y`, but growth stalls as it fills up, in proportion to the room left, `1 - y`:

```
y'(t) = y(t) (1 - y(t)),      y(0) = 0.1
```

This is the logistic equation, and it has the three ingredients every IVP has:

- **State**: `y(t)`. Here it is one number; in general it is a vector (a `float array` in the code).
- **Right-hand side**: `f(t, y)`, which says how fast the state changes given the time and the current state. Here `f(t, y) = y (1 - y)`; it happens not to use `t`.
- **Initial condition**: the state at the start, `y(0) = y_0 = 0.1`.

The solution is the function `y(t)` that starts at `y_0` and whose slope at every time equals `f(t, y(t))`. For this equation, and only because it is unusually friendly, we can write it down:

```
y(t) = 1 / (1 + 9 e^(-t))
```

At `t = 0` it gives `1 / (1 + 9) = 0.1`, and differentiating gives `y' = 9 e^(-t) / (1 + 9 e^(-t))²`, which equals `y (1 - y)`. The curve is an S: a slow start, fastest growth around `y = 1/2` (at `t = ln 9`, about 2.2), then flattening towards 1.

Most ODEs have no formula like this, which is why numerical methods exist: they compute approximate values of `y(t)` on a grid of times. The project's tests, called the **corpus** ([`test/corpus.ml`](../../test/corpus.ml)), run the code on four problems defined in [`test/problems.ml`](../../test/problems.ml): `Canary`, `Logistic`, `VanDerPol` and `Robertson`. `Canary` and `Logistic` have an `exact` function that gives the true solution, so the tests can measure the error of the code. For the other two the tests check different properties; [chapter 6](06-the-corpus.md) covers all four.

**Systems.** In general `y = (y_1, ..., y_m)` is a vector and `f` returns a vector of the same length. In the code `f` has the type `Ode.rhs = float -> Vec.t -> Vec.t`, where `Vec.t` is `float array` ([`src/ode.mli`](../../src/ode.mli), [`src/vec.mli`](../../src/vec.mli)). Read the arrows as: a function that takes a float (the time) and a vector (the state) and returns a vector (the derivative); [docs/ocaml.md](../ocaml.md) explains the notation. A whole problem is one record, `Ode.problem = { rhs; t0; t_end; y0 }`. The *canary* is a system of three independent equations packed into one vector:

```
y_1' = -1 y_1,     y_2' = -100 y_2,     y_3' = -10000 y_3,     y(0) = (1, 1, 1)
```

Its solution is `y_i(t) = e^(-λ_i t)` with the rates `λ = (1, 100, 1e4)` (in the code, `Canary.lambda`). All three components decay to zero, at wildly different speeds. We will come back to it.

## 2. Time stepping and explicit Euler

We cannot compute `y(t)` for every `t`. Instead we pick times `t_0 < t_1 < t_2 < ...` a step `h` apart and compute numbers `y_n` that approximate `y(t_n)`, starting from the known `y_0`.

At any point the ODE tells us the slope: `y'(t) = f(t, y(t))`. A straight line with that slope is a decent approximation of the curve for a short while. Follow it for time `h`, read the slope again at the new point, and repeat. That is **explicit Euler** (also called forward Euler):

```
y_{n+1} = y_n + h f(t_n, y_n)
```

"Explicit" means the right-hand side uses only values we already have, so `y_{n+1}` is computed directly, with no equation to solve.

By hand on the logistic problem with `h = 0.5`, using `f(y) = y (1 - y)`:

| n | t_n | y_n | f(y_n) | y_{n+1} = y_n + 0.5 f(y_n) | exact y(t_{n+1}) |
|---|-----|-----|--------|----------------------------|------------------|
| 0 | 0.0 | 0.1 | 0.09 | 0.145 | 0.1548281 |
| 1 | 0.5 | 0.145 | 0.123975 | 0.2069875 | 0.2319693 |
| 2 | 1.0 | 0.2069875 | 0.1641437 | 0.2890593 | 0.3324279 |

The error (exact minus computed) is 0.0098 after one step, 0.0250 after two and 0.0434 after three. Smaller steps shrink it. The next section says by how much, and why.

## 3. How wrong: local error, global error, order

**Local error.** Taylor's theorem says

```
y(t + h) = y(t) + h y'(t) + (h²/2) y''(t) + (h³/6) y'''(t) + ...
```

One Euler step started from the exact value `y(t)` produces `y(t) + h y'(t)`. So the **local error**, defined here as the exact value minus what one step computes when it starts on the exact solution, is

```
local error of Euler = (h²/2) y''(t) + O(h³)
```

`O(h³)` stands for all terms that shrink at least as fast as a constant times `h³` when `h` is small. (The local error is also called the local truncation error, LTE: the part of the error that comes from truncating the Taylor series.)

Check on `y' = -y`, `y(0) = 1`, `h = 0.1`. One step gives `1 + 0.1 · (-1) = 0.9`. The exact value is `e^(-0.1) = 0.9048374`, so the local error is 0.0048374. The formula predicts `(0.1² / 2) · y''(0) = 0.005 · 1 = 0.005`.

**Global error.** What we care about is the **global error**, the difference between `y_n` and `y(t_n)` at some fixed final time `T`. Reaching `T` takes `T / h` steps. Each one commits a local error of about `c h²`, and a stable method carries earlier errors forward without blowing them up, so the total is roughly `(T / h) · c h² = c T h`. The global error of Euler is proportional to `h`.

**Order.** A method has **order p** if its global error is `O(h^p)` (its local error is then `O(h^(p+1))`). Euler has order 1; [chapter 4](04-bdf.md) builds a method of order 2. The practical meaning: halving `h` divides the error by about `2^p`.

By hand on `y' = -y` up to `T = 1`: Euler gives `y_n = (1 - h)^n`, and `e^(-1) = 0.3678794`.

| h | steps | Euler value | error |
|---|-------|-------------|-------|
| 0.1 | 10 | 0.9^10 = 0.3486784 | 0.0192010 |
| 0.05 | 20 | 0.95^20 = 0.3584859 | 0.0093935 |

The ratio is 0.0192010 / 0.0093935 = 2.04, close to `2^1`. A ratio near 4 would indicate order 2; chapter 4 uses exactly this kind of ratio to test BDF2. Order describes the limit `h → 0` (at a given `h` the constant in front matters too), and "local errors add up" assumes errors are not amplified from step to step, which is a stability question, the next topic.

## 4. Stability: when a step is too long

The simplest decaying problem is the **test equation** `y' = -λ y` with `λ > 0`. Its solution `y(0) e^(-λ t)` decays to zero. Explicit Euler turns it into

```
y_{n+1} = y_n - h λ y_n = (1 - hλ) y_n          so          y_n = (1 - hλ)^n y_0
```

The numerical solution is multiplied by the same factor `1 - hλ` at every step. If the factor is in `[0, 1)` (`hλ <= 1`) the values decay smoothly, like the true solution. In `(-1, 0)` they decay but flip sign each step, which the true solution never does. At `-1` they bounce between `±y_0` forever. Below `-1` (`hλ > 2`) they grow without bound even though the true solution decays to zero. The numerical solution stays bounded exactly when `|1 - hλ| <= 1`, that is, when `h <= 2 / λ`. By hand with `λ = 100`, so the limit is `h = 0.02`:

| h | hλ | factor 1 - hλ | y_0, y_1, y_2, y_3, ... |
|---|----|---------------|--------------------------|
| 0.005 | 0.5 | 0.5 | 1, 0.5, 0.25, 0.125 |
| 0.01 | 1 | 0 | 1, 0, 0, 0 |
| 0.02 | 2 | -1 | 1, -1, 1, -1 |
| 0.05 | 5 | -4 | 1, -4, 16, -64 |

The last row is wrong in the worst way: the exact solution drops by `e^(-5) ≈ 0.0067` over that step, while Euler multiplies by `-4`. A step that is too long does not merely lose accuracy; it produces garbage that grows exponentially. This is a stability limit, and it has nothing to do with how accurate you want the answer to be.

**Systems, eigenvalues and the Jacobian.** Take a linear system `y' = A y`, where `A` is a matrix. An **eigenvector** of `A` is a nonzero vector `v` with `A v = λ v` for some number `λ`, its **eigenvalue**: `A` merely stretches `v` by the factor `λ`. A state that is a multiple of `v`, `y = c v`, then obeys `c' = λ c`. For a decaying problem `λ` is negative, so this is the test equation with rate `-λ`, and explicit Euler multiplies `c` by `1 + hλ = 1 - h|λ|` at every step. For a diagonal matrix the eigenvalues are the diagonal entries and the eigenvectors are the coordinate axes. In general eigenvalues can be complex numbers (section 6). For most matrices every state is a sum of multiples of eigenvectors, so the pieces evolve independently and the step must satisfy the limit for every eigenvalue.

For a nonlinear `f` nothing is exactly linear, but small disturbances of a state `y` evolve approximately like a linear system whose matrix is the **Jacobian** of `f` at `y`: the matrix of partial derivatives `J.(i).(j) = ∂f_i / ∂y_j`, the derivative of component `i` of `f` with respect to component `j` of `y` (row `i` is the output, column `j` the input). This is a rule of thumb rather than a theorem, but it is the one that guides practice. The code computes `J` numerically; see [chapter 3](03-jacobians-and-floating-point.md).

## 5. Stiffness

Return to the canary, `y' = -Λ y`, where `Λ = diag(1, 100, 1e4)` is the diagonal matrix with those entries. Its Jacobian is `-Λ`, so its eigenvalues are minus the three rates, and explicit Euler is stable only if all three limits hold:

```
component 1:  h <= 2 / 1     = 2
component 2:  h <= 2 / 100   = 0.02
component 3:  h <= 2 / 1e4   = 2e-4
```

The fastest component sets the limit for the whole system: `h <= 2e-4`.

Now look at what that component actually does. The **time scale** of a decay `e^(-λt)` is `1/λ`, the time in which it shrinks by a factor `e`, so the canary's time scales are `1`, `0.01` and `1e-4`. `y_3(t) = e^(-1e4 t)` is already down to `e^(-10) ≈ 4.5e-5` at `t = 0.001`, and `y_2` reaches the same size at `t = 0.1`. After that, what remains, `(e^(-t), ≈0, ≈0)`, is smooth and changes on a time scale of 1. Accuracy alone would allow steps far longer than `2e-4`. Yet explicit Euler must keep taking steps `1 / 2e-4 = 5000` times shorter than the time scale of what is left, because any tiny disturbance in `y_3` (a rounding error, or coupling to other components in a problem that is not diagonal) is multiplied by `1 - h λ_3` at every step. With `h = 0.01` that factor is `-99`: after 100 steps a disturbance has grown by `99^100 ≈ 3.7e199`.

That is **stiffness**. There is no single sharp definition; the working one used in this project is:

> A problem is stiff, on a given interval and for a given accuracy, when the step size an explicit method needs for stability is much smaller than the step size accuracy alone would require.

Typical symptoms are Jacobian eigenvalues with negative real parts and very different magnitudes: a fast decaying mode that dies quickly, next to the slow motion we actually want to follow. The canary has a ratio of `1e4` between its fastest and slowest rate. While a fast transient is still alive the step must be small anyway, so stiffness is about the stretch of the solution *after* the fast modes have died out. Because it depends on the interval and on the accuracy asked for, it is a property of the problem together with how you use it, not of the equation alone.

The corpus has two more stiff problems. Robertson's chemical kinetics has rate constants from 0.04 to 3e7, a ratio of 7.5e8. For van der Pol with `μ = 1000`, the right-hand side is `f(t, y) = (y_2, μ (1 - y_1²) y_2 - y_1)`; its Jacobian at the initial state `y = (2, 0)` is `[[0, 1], [-1, -3000]]`, with eigenvalues about `-3000` and about `-0.00033`. An explicit method would be limited to `h < 2 / 3000 ≈ 6.7e-4` there. See [chapter 6](06-the-corpus.md).

## 6. Backward Euler

Explicit Euler reads the slope at the *start* of the step. **Backward Euler** reads it at the *end*:

```
y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})
```

The unknown `y_{n+1}` now appears on both sides, so the method is **implicit**: each step requires solving an equation. On the test equation `y' = -λ y` the equation can be solved by hand:

```
y_{n+1} = y_n - h λ y_{n+1}          so          y_{n+1} = y_n / (1 + hλ)
```

For every `h > 0` the factor `1 / (1 + hλ)` lies in `(0, 1)`. The numerical solution decays whatever the step size; there is no stability limit. Compare the two methods by hand on `y' = -1e4 y` (the canary's fastest component) with `h = 0.01`, so `hλ = 100`:

| n | t_n | exact | explicit Euler (factor -99) | backward Euler (factor 1/101) |
|---|-----|-------|-----------------------------|-------------------------------|
| 0 | 0 | 1 | 1 | 1 |
| 1 | 0.01 | 3.7e-44 | -99 | 0.0099 |
| 2 | 0.02 | 1.4e-87 | 9801 | 9.8e-5 |
| 3 | 0.03 | 5.1e-131 | -970299 | 9.7e-7 |

Explicit Euler explodes. Backward Euler is not accurate on this component (0.0099 against `3.7e-44`), but it does what matters: the component is damped by a factor 101 per step and is irrelevant anyway once it has decayed. The step can therefore be chosen to suit the slow component alone. The reason it works is that the slope is read at the new point, where it already points back towards zero, the value the exact solution decays to; the larger `h`, the stronger the damping.

**A-stability (informally).** Real systems often have components that oscillate while they decay, such as a damped spring. Their eigenvalues are complex numbers `λ = a + ib` with `a < 0`, and the exact solution `e^(λt)` decays. Apply a one-step method such as Euler to `y' = λ y` with complex `λ` and it multiplies `y_n` by a factor `R` that depends only on the product `hλ`. The method is **A-stable** if `|R| <= 1` for every `h > 0` and every `λ` with negative real part. Backward Euler is A-stable: `R = 1 / (1 - hλ)`, and `|1 - hλ| >= 1` whenever the real part of `hλ` is negative. Explicit Euler is not: `R = 1 + hλ`, which has magnitude at most 1 only when `hλ` lies in the disc `|1 + hλ| <= 1`. (A multistep method such as BDF2 has a recurrence instead of a single factor, so its definition is phrased differently, but the meaning is the same: every numerical solution decays. [Chapter 4](04-bdf.md) returns to this.) On a stiff problem an A-stable method lets accuracy, not stability, decide the step size.

**Stable is not accurate.** Backward Euler is only first order. On the canary's slow component it computes `y_n = (1 + h)^(-n)`. With `n = 1/h` steps to reach `t = 1` that is about `e^(-1) (1 + h/2)`, so the error at `t = 1` is about `(h/2) e^(-1) ≈ 0.184 h`. The corpus line `bdf1 canary t=1 dt=2e-6: max error 3.68e-07 < 1e-06: true` is this formula at `h = 2e-6`. It also takes `1 / 2e-6 = 500,000` steps, which is why it dominates the run time of the test suite. Stability is not what forces the tiny step there. Accuracy is: a first-order method needs `0.184 h < 1e-6`, so `h` below about `5.4e-6`. Getting more accuracy per step is the motivation for BDF2 in [chapter 4](04-bdf.md).

## 7. The price: an equation to solve at every step

Backward Euler defines `y_{n+1}` through the equation `x = y_n + h f(t_{n+1}, x)`. Moving everything to one side gives the **residual** `G(x)`, which is zero exactly when `x` solves the equation:

```
G(x) = x - y_n - h f(t_{n+1}, x) = 0
```

How hard that is depends on `f`:

- **Linear `f`**, like the canary: `(I + hΛ) x = y_n`, where `I` is the identity matrix. It is a linear system (here diagonal), solved by Gaussian elimination.
- **Logistic**: `x - y_n - h x (1 - x) = 0` is the quadratic `h x² + (1 - h) x - y_n = 0`. For `h = 0.5` and `y_n = 0.1` it reads `x² + x - 0.2 = 0`, so `x = (-1 + √1.8) / 2 = 0.1708204`. The exact value is 0.1548281 and explicit Euler gave 0.145, so the two errors have opposite signs: backward Euler overshoots by 0.0160 where explicit Euler undershoots by 0.0098. [Chapter 5](05-step-control.md) uses that to estimate the error.
- **Anything else** (Robertson, van der Pol): no closed form. The equation has to be solved numerically, by **Newton's method**, the subject of [chapter 2](02-newton.md). Newton needs the derivative of `G`, which is `I - h J(t_{n+1}, x)` with `J` the Jacobian of `f`.

An implicit step therefore costs several Newton iterations, each needing evaluations of `f`, a Jacobian and a linear solve: far more than one explicit step. It pays off when the step can be thousands of times longer. For a non-stiff problem an explicit method is usually the cheaper choice.

## Running the snippets

The snippets in these chapters call the library directly. Build a small probe project outside the repository that links to it. Run this from the repository root, with your opam switch active ([docs/ocaml.md](../ocaml.md) explains opam and dune):

```sh
REPO=$(pwd)
mkdir -p ../vstiff-scratch/p && cd ../vstiff-scratch
printf '(lang dune 3.0)\n' > dune-project
ln -sfn "$REPO/src" src
ln -sfn "$REPO/test/problems.ml" p/problems.ml
printf '(executable (name probe) (libraries vstiff))\n' > p/dune
# then, for each snippet: put it in p/probe.ml, build and run
dune build --root . ./p/probe.exe && ./_build/default/p/probe.exe
```

Nothing is written inside the repository. The `problems.ml` link makes the corpus problems (`Problems.Canary` and friends) available to your snippet. Dune's default profile turns warnings into errors, and an unused top-level definition is one of them: delete it or use it. Later chapters reuse this setup; the section Probes of [docs/testing.md](../testing.md) describes the same setup and more ways to use it.

Explicit Euler is not part of the library, but it fits the contract `Ode.Method` in a few lines, and then `Stepper.fixed` can drive it like any other method:

```ocaml
open Vstiff

(* Explicit Euler: y_{n+1} = y_n + h f(t_n, y_n), with nothing to remember. *)
module Explicit_euler : Ode.Method = struct
  type history = unit

  let start = ()
  let step rhs h () (at : Ode.point) = Ok (Vec.axpy h (rhs at.t at.y) at.y, ())
end

let () =
  let open Problems.Canary in
  List.iter
    (fun dt ->
      match Stepper.fixed (module Explicit_euler) ~dt problem with
      | Ok y -> Printf.printf "explicit Euler, dt = %-8g y_3(1) = %g   (exact %g)\n" dt y.(2) (exact 1.).(2)
      | Error e -> print_endline (Fail.to_string e))
    [ 1.9e-4; 2.1e-4; 0.01 ]
```

A few OCaml notes (the full story is in [docs/ocaml.md](../ocaml.md)). `module Explicit_euler : Ode.Method = struct ... end` defines a module and promises that it provides what the contract in [`src/ode.mli`](../../src/ode.mli) lists: a `history` type, a value `start`, and a function `step`. The history is whatever a method remembers from one step to the next (BDF2 in chapter 4 remembers the previous step); Euler remembers nothing, so its history is `unit`, the type whose only value is `()`. `step rhs h () at` receives the right-hand side, the step, the history and the current point `at = { t; y }`, and returns `Ok (new state, next history)`; it is `Ok` because a step can fail in general, but this one never does. `Vec.axpy h x y` computes `h x + y` for vectors. `(module Explicit_euler)` hands the module to `Stepper.fixed` as an argument. `let open Problems.Canary in` brings `problem`, `y0` and `exact` into scope inside that block; `y.(2)` is the third element of an array (counting from 0); `%g` prints a float and `%-8g` pads it to 8 characters.

Per step, `y_3` is multiplied by about `-0.9`, `-1.1` and `-99` (work it out from `1 - h λ_3`; `Stepper.fixed` adjusts `dt` slightly so that a whole number of steps fits). So the first run should decay while flipping sign, and the other two should explode, even though the true `y_3(1) = e^(-10000)` is so small that it prints as exactly 0: it underflows, because the smallest positive float is about `5e-324`. To see backward Euler stay tame at a step that destroys explicit Euler, append this to the same file; both `let () =` blocks run, in order:

```ocaml
let () =
  let open Problems.Canary in
  match Stepper.fixed (module Bdf1) ~dt:0.01 problem with
  | Ok y ->
      Array.iteri (fun i v -> Printf.printf "backward Euler, dt = 0.01: y_%d(1) = %g   (exact %g)\n" (i + 1) v (exact 1.).(i)) y
  | Error e -> print_endline (Fail.to_string e)
```

All three components stay finite. The slow one is off by about `0.184 h ≈ 1.8e-3`, as section 6 predicts; the other two are tiny, as they should be.

## In the code

| Idea | Where |
|------|-------|
| State, right-hand side, problem | `Ode.rhs`, `Ode.problem` and `Ode.point` in [`src/ode.mli`](../../src/ode.mli); `Vec.t` is `float array` ([`src/vec.mli`](../../src/vec.mli)) |
| Corpus ODEs with their exact solutions | `Canary` and `Logistic` in [`test/problems.ml`](../../test/problems.ml): `rhs`, `y0`, `exact` and a `problem` record |
| What a method must provide | `Ode.Method` in [`src/ode.mli`](../../src/ode.mli): `history`, `start` and `step rhs h history at`, which returns the new state and the next history, or a `Fail.t` |
| Time stepping with a constant step | `Stepper.fixed (module M) ~dt problem` in [`src/stepper.ml`](../../src/stepper.ml): `n = round((t_end - t0) / dt)` steps (at least one for a non-empty span) of `h = (t_end - t0) / n`, which equals `dt` only when `dt` divides the span: step `k` ends at `t0 + k h`, the last at `t_end`, and the method is given the difference of the end times; it returns the final state or the first failure |
| One backward Euler step | `Bdf1.step` in [`src/bdf1.ml`](../../src/bdf1.ml), a `Method` whose history is `unit` inside the module (the interface keeps it abstract) |
| The implicit equation `G(x) = 0` and its Jacobian `I - h J` | `Stage.solve` in [`src/stage.ml`](../../src/stage.ml), which hands both to `Newton.solve` ([`src/newton.ml`](../../src/newton.ml)) |
| Failure of the implicit solve | `Fail.t` in [`src/fail.mli`](../../src/fail.mli), as `Error Diverged` or `Error Nan` |
| Explicit Euler | Not a library method. The snippet above writes it as one; the expression `Vec.axpy h (rhs t y) y` also appears in `Bdf2.step_with_error`, on the first step, as part of the error estimate ([chapter 5](05-step-control.md)) |
| The first-order accuracy check | The `bdf1 canary` line of [`test/corpus.ml`](../../test/corpus.ml); expected output in [`test/corpus.expected`](../../test/corpus.expected) |

`Bdf1.step` is a single expression:

```text
let step rhs h () at =
  Result.map (fun y -> (y, ())) (Stage.solve rhs { Stage.t = at.t +. h; gamma = h; psi = at.y } at.y)
```

Read it as: solve the stage equation at the new time `t + h`, with `psi = y_n` (the known part) and `gamma = h` (the multiplier of `f`), starting Newton's iteration from `y_n` itself; `Result.map f r` applies `f` inside an `Ok` and passes an `Error` through, here pairing the answer with the empty history. [Chapter 4](04-bdf.md) explains why every BDF step has that shape.

## Check yourself

1. **What are the three ingredients of an initial value problem, and which OCaml type holds the second?**
   The state `y`, the right-hand side `f(t, y)` and the initial condition `y(t_0) = y_0`. `f` is an `Ode.rhs`.

2. **One explicit Euler step from the exact solution of `y' = -y` with `h = 0.1`: how big is the local error?**
   About `(h²/2) y'' = 0.005`; the exact local error is 0.0048374.

3. **Euler with `h = 0.025` instead of `0.05` on a smooth problem: what error ratio do you expect?**
   About 2, since Euler has order 1. Order 2 would give about 4.

4. **What is the largest stable step for explicit Euler on `y' = -500 y`? What happens at `h = 0.01`?**
   `h <= 2 / 500 = 0.004`. At `h = 0.01`, `hλ = 5` and the factor is `-4`: the values grow by 4 in magnitude each step while alternating in sign.

5. **In the canary, which component limits explicit Euler, and why does the limit remain after that component has decayed?**
   The `λ = 1e4` component, giving `h <= 2e-4`. Any disturbance in it, such as a rounding error, is multiplied by `1 - hλ` per step and grows whenever `|1 - hλ| > 1`.

6. **Why is "stiff" not a property of the equation alone?**
   It compares the stability limit with the step accuracy would allow, and the latter depends on the interval and the accuracy required. While a fast transient is alive the step must be small anyway.

7. **Backward Euler gives 0.0099 where the exact value is `3.7e-44` (the `h = 0.01` table). Why is that acceptable?**
   The component is irrelevant once decayed, and backward Euler damps it by a factor 101 per step instead of letting it grow. It would be a problem only if we needed that component accurately, which would force a smaller step.

8. **Backward Euler is stable for any `h`. Why does the corpus canary test still need 500,000 steps?**
   Backward Euler is first order with error about `0.184 h` on the slow component, and the test demands an error below `1e-6`. Accuracy, not stability, forces `h` near `2e-6`.

9. **What is the residual `G` for a backward Euler step, and what is its Jacobian?**
   `G(x) = x - y_n - h f(t_{n+1}, x)` with Jacobian `I - h J`. See `Stage.solve` in [`src/stage.ml`](../../src/stage.ml).

Next: [2. Newton's method](02-newton.md), how the implicit equation of each step is solved. (On the onboarding plan you have read chapters 2 and 3 already: continue with [4. BDF methods](04-bdf.md).)
