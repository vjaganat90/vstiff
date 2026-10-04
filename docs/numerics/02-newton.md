# 2. Newton's method

Chapter 2 of 6 in the numerical-methods track. Previous: [1. ODEs and stiffness](01-odes-and-stiffness.md). Next: [3. Jacobians and floating point](03-jacobians-and-floating-point.md). Index and reading order: [docs/README.md](../README.md). Terms and symbols: [glossary](../glossary.md).

**Summary.** Newton's method solves `G(x) = 0` by replacing `G` with its tangent line again and again. This chapter shows why that finds a square root in five steps and when it fails (a flat tangent), how damping and the Armijo test rescue a step that overshoots, and how the idea carries over to many unknowns, where every step solves a linear system by Gaussian elimination with partial pivoting. It then covers singular matrices, the stopping test (why it looks at the step, not the residual, and the overflow it still lets through) and failure as `Fail.t` values, never exceptions. Finally it applies all this to the stage equation of every BDF step and explains why a wrong Jacobian makes Newton slower without changing the root it finds.

The snippets run in the probe project ("Setup" in [exercises.md](../exercises.md), or "Running the snippets" in [chapter 1](01-odes-and-stiffness.md)); most of them call the numerical kernel, the library `Numerics`, and start with `open Numerics`.

**Names.** Here `G` is the function whose root we want, and `G(x)` is the **residual**: zero exactly when `x` solves the equation, otherwise a measure of how far `x` is from doing so. `Newton.solve f jac x0` calls `G` by the name `f` and its derivative by `jac`; `Stage.solve` calls them `residual` and `jacobian`, because there `f` is the ODE's right-hand side. `G'(x)` is the derivative; for a vector `G` it is the matrix of entries `∂G_i/∂x_j` (row `i` an equation, column `j` an unknown). `|v|_inf` is the largest absolute component of `v`. Newton's own tolerance, `1e-10`, is a private constant, unrelated to the `tol` that decides whether the integrator accepts a step ([chapter 5](05-step-control.md)).

## 1. One unknown

### The tangent-line idea

Draw the graph of `G` and pick a starting guess `x_0`. Draw the tangent line there, `L(x) = G(x_0) + G'(x_0) (x - x_0)`. A curve is not a line, but near `x_0` it is close to one, and the line is easy: it crosses zero at

```
x_1 = x_0 - G(x_0) / G'(x_0)
```

Take `x_1` as the new guess and repeat: `x_{k+1} = x_k - G(x_k) / G'(x_k)`. That is Newton's method. The correction `dx = -G(x_k) / G'(x_k)` is the Newton step, the successive guesses are the **iterates**, and the method **converges** when they settle on a root.

### Worked example: the square root of 2

Take `G(x) = x² - 2`, so `G'(x) = 2x`, and start at `x_0 = 1`. Then `G(1) = -1` and `G'(1) = 2`, so `x_1 = 1 - (-1)/2 = 1.5`; at `x_1`, `G = 0.25` and `G' = 3`, so `x_2 = 1.5 - 0.25/3 = 1.41666...`. The library does the same arithmetic. `jac` is called once per iteration, so this program prints the point where each iteration starts, then the result; the `x_k` column of the table below is exactly its output:

```ocaml
open Numerics

let () =
  let f x = [| (x.(0) *. x.(0)) -. 2. |] in
  let jac x =
    Printf.printf "iterate %.17g\n" x.(0);
    [| [| 2. *. x.(0) |] |]
  in
  match Newton.solve f jac [| 1. |] with
  | Ok x -> Printf.printf "result  %.17g\n" x.(0)
  | Error e -> print_endline (Fail.to_string e)
```

OCaml notes ([docs/ocaml.md](../ocaml.md) has the details): `[| 1.; 2. |]` is an array of floats and `x.(0)` its first element; a matrix is an array of rows, so `a.(i).(j)` is row `i`, column `j`; float arithmetic has dotted operators (`+.`, `*.`); `Ok v` and `Error e` are the two shapes of a `result`, taken apart with `match`; `Printf.printf` takes a C-style format (`%g` short, `%.17g` with 17 significant digits, enough to tell any two floats apart).

With `e_k = |x_k - √2|` (the result, `x_5`, is the float closest to `√2`):

| k | x_k | e_k | e_{k+1} / e_k² |
|---|-----|-----|----------------|
| 0 | 1 | 4.1e-1 | 0.50 |
| 1 | 1.5 | 8.6e-2 | 0.33 |
| 2 | 1.4166666666666667 | 2.5e-3 | 0.35 |
| 3 | 1.4142156862745099 | 2.1e-6 | 0.35 |
| 4 | 1.4142135623746899 | 1.6e-12 | |
| 5 | 1.4142135623730951 | about 1e-16 (rounding) | |

The number of correct digits, `-log10(e_k)`, goes 0.4, 1.1, 2.6, 5.7, 11.8: it roughly doubles per step. This is **quadratic convergence**, and the last column shows why: each error is a constant times the square of the previous one. For this `G` the identity is exact, `x_{k+1} - √2 = (x_k - √2)² / (2 x_k)`, so the ratio is `1 / (2 x_k)`, settling at `1 / (2√2) ≈ 0.354`. For any `G` and a root `r` with `G'(r) ≠ 0`, Taylor's theorem gives the same shape: expand `0 = G(r) = G(x_k) + G'(x_k)(r - x_k) + ½ G''(ξ)(r - x_k)²` for some `ξ` between `x_k` and `r`, divide by `G'(x_k)` and rearrange to `x_{k+1} - r = (G''(ξ) / (2 G'(x_k))) (x_k - r)²`. The constant is small when `G` is nearly straight (`G''` small) and the tangent is not flat (`G'` not small). The catch is hidden in "near": the formula promises nothing for a bad `x_0`.

### When the tangent is flat

At `x_0 = 0`, `G'(0) = 0`: the tangent of `x² - 2` is the horizontal line `y = -2`, which never crosses zero. The formula divides by zero; in the code that division is a linear solve with the 1×1 matrix `[[0]]`, which `Linalg.solve` reports as `None` (section 5) and `Newton.solve` turns into `Error Diverged`. That is the corpus line `newton quadratic x0=0: Error Diverged`.

## 2. Overshoot and damping

A tangent is a good model only near the point where it is drawn. When `G'` is small the Newton step is huge and the new point may be worse than the old one. Take `G(x) = atan x`, whose root is 0. At `x_0 = 2`, `G = 1.107` and `G' = 1 / (1 + 4) = 0.2`, so the step is `-1.107 / 0.2 = -5.54` and the full step lands at `x = -3.54`, where `|G| = 1.295`: farther from zero than before. The plain formula gives `2, -3.54, 13.95, -279, 122017, ...`, growing without bound.

The safeguard is **damping**, a simple **line search**: if the full step is not good, try half of it, then a quarter, and so on. "Good" is judged by the one thing Newton can measure, the residual (it cannot measure the distance to a root it does not know). The code accepts the step `λ dx` when

```
|G(x + λ dx)|_inf  <=  (1 - 1e-4 λ) |G(x)|_inf
```

This is the **Armijo test**. If `G` were exactly linear, `G(x + λ dx)` would equal `(1 - λ) G(x)`, so the residual would drop by the fraction `λ`. The test demands only `1e-4` of that promised drop: almost any real decrease passes, an increase or a stall fails. The factors tried are `λ = 1, 1/2, ..., 1/1024` (eleven; `1/1024` is `min_damping`), and if none passes the result is `Error Diverged`. A trial point where `G` is not finite counts as a failure too. In [`src/numerics/newton.ml`](../../src/numerics/newton.ml) the test is one line of the local function `damp`:

```text
if Vec.finite fx' && Vec.norm_inf fx' <= (1. -. (armijo *. lambda)) *. r then iterate (k + 1) x' fx'
else damp (lambda /. 2.)
```

`fx'` is `G(x + λ dx)` (a prime is an ordinary part of an OCaml name), `r` is `|G(x)|_inf` and `armijo` is `1e-4`. An accepted point goes straight into the next iteration together with its residual, which is not computed again. This program prints every point where `G` is evaluated for `atan` from 2:

```ocaml
open Numerics

let () =
  let f x =
    Printf.printf "G(% .6f) = % .6f\n" x.(0) (atan x.(0));
    [| atan x.(0) |]
  in
  let jac x = [| [| 1. /. (1. +. (x.(0) *. x.(0))) |] |] in
  match Newton.solve f jac [| 2. |] with
  | Ok x -> Printf.printf "root %g\n" x.(0)
  | Error e -> print_endline (Fail.to_string e)
```

```text
G( 2.000000) =  1.107149
G(-3.535744) = -1.295169
G(-0.767872) = -0.654841
G( 0.273082) =  0.266582
G(-0.013380) = -0.013379
G( 0.000002) =  0.000002
G(-0.000000) = -0.000000
root 0
```

Line 1 is the start. Line 2 is the full step: `1.295 > (1 - 1e-4) · 1.107`, rejected. Line 3 is the half step: `0.655 <= (1 - 5e-5) · 1.107`, accepted. From there every full step is accepted on the first try. A second example, `G(x) = sqrt x - 1` from 9: the Newton step is `-12`, the full step lands at `x = -3`, and `sqrt` returns `nan` there instead of raising; a `nan` residual fails the test, so the point is refused and the half step, `x = 3`, is taken. The corpus pins this with `log x` from 3 (line 47 of [chapter 6](06-the-corpus.md)): the full step, `-3 ln 3`, lands at `-0.30`, `Float.log` returns `nan` there, the point is refused and the half step, `1.35`, is taken.

The corpus pins damping with a start just past the point where plain Newton fails on `atan`. A full Newton step takes `x` to `x - (1 + x²) atan x`, which is `-x` at `x ≈ 1.39` (solve `2x = (1 + x²) atan x`) and lands farther out beyond it; from 1.5 the plain iterates are `-1.69, 2.32, -5.11, 32.3, -1575, ...`, growing until `x²` overflows, the Jacobian is exactly 0 and the solve gives up with `Error Diverged`. With damping the half step passes (`|G|` falls from `0.983` to `0.097`) and the line `newton atan x0=1.5 (needs damping)` prints `Ok` with the root 0.

Why can a small enough `λ` pass? The Newton direction satisfies `G'(x) dx = -G(x)`, so near `x` the residual along it really does shrink like `(1 - λ)`. If none of the eleven factors passes, the Jacobian does not match `G`, or `G` is not smooth, or the step is so long that the linear model holds only for a `λ` below `1/1024` (section 5 has an example). Damping then gives up with `Diverged` instead of looping.

## 3. Several unknowns

With `n` unknowns `G` maps a vector to a vector, and the derivative becomes the Jacobian matrix `G'(x)`. Near `x` the best linear model is `G(x + dx) ≈ G(x) + G'(x) dx`; setting it to zero gives the Newton step as a linear system:

```
G'(x) dx = -G(x)          then   x_new = x + λ dx
```

For `n = 1` this is the formula of section 1. (Textbooks usually call the matrix `J`; this track keeps `J` for the Jacobian of the ODE's right-hand side `f`, so the matrix of `G` is written `G'`.) Each iteration costs a Jacobian, a linear solve and at least one evaluation of `G`. The corpus' first Newton case is `G(x) = A x - b` with `A = [[3, 1], [1, 2]]` and `b = (9, 8)`, the system `3u + v = 9`, `u + 2v = 8` for the unknown vector `(u, v)`, started at `(0, 0)` where `G = (-9, -8)`. Since `G` is linear, its tangent plane is `G` itself, and one step lands on the answer `(2, 3)`. The corpus line is `newton linear 2d: Ok [2.000000000000; 3.000000000000]`.

## 4. The linear solve: Gaussian elimination with partial pivoting

To solve `A dx = r` (here `r = -G(x)`), eliminate one unknown at a time, then substitute back. `Linalg.solve` first glues each right-hand-side entry to the end of its row, giving **augmented rows** `[coefficients | rhs]`. For the corpus system, whose first step from `(0, 0)` has `r = (9, 8)`:

```
row 0:   3  1 | 9
row 1:   1  2 | 8
```

1. **Pivot.** Look down the first column for the entry of largest absolute value: `3`, in row 0 (the row nearer the top wins a tie). That row is the **pivot row**, `3` the **pivot**. No exchange is needed.
2. **Reduce.** For every other row take the **multiplier** `m` = (its first entry) / pivot, subtract `m` times the pivot row, and drop the first column, which is now zero. Row 1: `m = 1/3` and `[2 - (1/3)·1 | 8 - (1/3)·9] = [5/3 | 5]`.
3. **Recurse.** What is left is the same problem with one unknown fewer, here the single row `[5/3 | 5]`: `v = 5 / (5/3) = 3`.
4. **Back-substitute.** The pivot row says `3u + 1·v = 9`, so `u = (9 - 1·3) / 3 = 2`.

Why the largest entry? Dividing by a small pivot magnifies rounding errors. Take `1e-20 u + v = 1`, `u + v = 2`, whose solution is very close to `(1, 1)`. With the tiny coefficient as pivot the multiplier is `1e20`; row 1 becomes `(1 - 1e20) v = 2 - 1e20`, and in floating point the `1` and the `2` are lost next to `1e20`, so both sides are `-1e20` and `v = 1`. Back-substitution then gives `u = (1 - 1) / 1e-20 = 0`: completely wrong. Choosing the larger entry (`1`, in row 1) as pivot keeps every multiplier at most 1 in size, and the library returns `[1; 1]` (the corpus line `linalg tiny leading pivot`).

A row exchange is needed when the corner is zero. The corpus line `linalg zero leading pivot` solves the augmented rows `[0 1 1 | 5]`, `[2 1 0 | 4]`, `[1 0 3 | 10]`. The first column is `(0, 2, 1)`, so the pivot is `2` in row 1 and rows 0 and 1 trade places. Row 0 has `m = 0` and becomes `[1 1 | 5]`; row 2 has `m = 1/2` and becomes `[-1/2 3 | 8]`. At the next level the pivot is `1`, `m = -1/2`, and `[3 + 1/2 | 8 + 5/2] = [7/2 | 21/2]` gives `z = 3`; then `y = (5 - 1·3) / 1 = 2` and `x = (4 - (1·2 + 0·3)) / 2 = 1`. The unknowns are never reordered, so there is nothing to undo at the end.

The code is this recursion, `solve_augmented` in [`src/numerics/linalg.ml`](../../src/numerics/linalg.ml). `pivot` returns the index of the pivot row and `swapped p i` says which original row sits at position `i` once rows 0 and `p` trade places, so the remaining rows keep their order except that row 0 takes the place of row `p`. The multiplier is `m` inside `reduce`, and the back-substitution is the `Option.map` at the end. Nothing is written into an array (a project rule, H6 in [AGENTS.md](../../AGENTS.md)), so each level allocates fresh ones, which does not matter for systems this small. This program runs the examples above (the second is the first with its equations in the other order) and two singular ones that section 5 explains:

```ocaml
open Numerics

let show = function
  | Some x -> print_endline (String.concat "; " (Array.to_list (Array.map (Printf.sprintf "%g") x)))
  | None -> print_endline "None"

let () =
  show (Linalg.solve [| [| 3.; 1. |]; [| 1.; 2. |] |] [| 9.; 8. |]);
  show (Linalg.solve [| [| 1.; 2. |]; [| 3.; 1. |] |] [| 8.; 9. |]);
  show (Linalg.solve [| [| 0.; 1.; 1. |]; [| 2.; 1.; 0. |]; [| 1.; 0.; 3. |] |] [| 5.; 4.; 10. |]);
  show (Linalg.solve [| [| 1.; 2. |]; [| 2.; 4. |] |] [| 3.; 6. |]);
  show (Linalg.solve [| [| 0.3; 0.1 |]; [| 3.; 1. |] |] [| 1.; 0. |])
```

```text
2; 3
2; 3
1; 2; 3
None
-2.40192e+16; 7.20576e+16
```

## 5. Singular and nearly singular matrices

A matrix is **singular** when the system has no unique solution: one equation is a combination of the others (`[[1, 2], [2, 4]]`: row 1 is twice row 0). `Linalg.solve` detects it one way only: when the best pivot available is exactly `0.` it returns `None`, and `Newton.solve` turns `None` into `Error Diverged`. In the fourth line above the elimination leaves the row `[0 | 0]`.

Rounding makes "exactly" fragile. The last line comes from `[[0.3, 0.1], [3, 1]]`, in which row 1 is exactly ten times row 0 in decimal, so the matrix is singular. But `0.3` and `0.1` are not exactly representable in binary: the multiplier is the rounded `0.3 / 3`, the entry that should cancel to `0` leaves about `1.4e-17`, `solve` divides by that and returns numbers of size `1e16`.

A matrix that is merely close to singular is just as risky. For `[[1, 1], [1, 1 + 1e-12]]` the right-hand side `(2, 2 + 1e-12)` gives `(1, 1)` and `(2, 2)` gives `(2, 0)`. A change of `1e-12` in the data moves the answer by 1, so rounding errors of relative size `1e-16` can move it by about `1e-4`: about twelve of the sixteen digits are gone. The factor by which a matrix can magnify relative errors in the data is its **condition number** (about `4e12` here, infinite for a singular matrix).

The library does not estimate it, so Newton's protection is only partial. A step containing NaN or infinity ends the solve with `Error Diverged`; any other step, however large, is judged by the Armijo test like every step, and an accepted step has reduced the residual. That does not rule out nonsense. Take the first Newton case of the corpus (section 3) with the matrix above and `b = (1, 0)`. On paper it has no solution, because row 1 of the matrix is ten times row 0 but `b_1` is not ten times `b_0`:

```ocaml
open Numerics

let () =
  let a = [| [| 0.3; 0.1 |]; [| 3.; 1. |] |] and b = [| 1.; 0. |] in
  let f x = Vec.sub (Array.map (fun row -> Vec.dot row x) a) b in
  match Newton.solve f (fun _ -> a) [| 0.; 0. |] with
  | Ok x -> Printf.printf "Ok [%g; %g]\n" x.(0) x.(1)
  | Error e -> print_endline (Fail.to_string e)
```

```text
Ok [-2.40192e+16; 7.20576e+16]
```

Newton reports success. The terms of `A x` are about `7e15` in row 0 and `7e16` in row 1, where neighbouring floats are 1 and 16 apart, so rounding the terms swallows the true residual (for the matrix as stored), about `(0.07, 4)`. The computed residual is exactly `0`, and a zero residual ends the solve with `Ok`. A Jacobian that is singular on paper is a bug to find upstream; the solver neither measures nor repairs it.

When `G` is strongly nonlinear a huge step usually lands where `G` is huge and is refused, and damping can reach its limit. For `G(x) = x² - 2` the step from `x_0` is about `1 / x_0`. From `x_0 = 1e-3` it is about `+1000`; the new point needs `|G| <= 2`, that is `x <= 2`, so `λ <= 2e-3`, and the first factor that small is `1/512`, the tenth trial, landing at `1.954`; Newton converges from there. This is corpus line 45, which prints `Error Diverged` if the smallest factor is `1/2`. From `x_0 = 1e-4` the step is about `+10000`, `λ` would have to be at most `2e-4`, below the smallest factor `1/1024`, and the solver gives up with `Error Diverged` although a perfectly good root sits at `1.414`.

## 6. When to stop

In every iteration `Newton.solve` checks these conditions, in this order (`k` counts iterations from 0, `norm_inf` is `|v|_inf`):

| Condition | Result |
|-----------|--------|
| `x` or `G(x)` contains NaN or infinity | `Error Nan` |
| `norm_inf(G(x))` is exactly `0` | `Ok x` |
| `k` has reached 50 | `Error Diverged` |
| the linear solve returns `None`, or the step `dx` contains NaN or infinity | `Error Diverged` |
| `norm_inf(dx) <= 1e-10 (1 + norm_inf(x))` | `Ok (x + dx)`, or `Error Nan` if `x + dx` is not finite |
| otherwise | damped step (section 2), then iteration `k + 1`; `Error Diverged` if damping gives up |

The converged result is `x + dx`: the tiny step is taken first. The factor `1 + |x|_inf` makes the test absolute when `x` is small and relative when `x` is large, and relative has a limit: next to the largest float, `max_float` (about `1.8e308`), a step of up to about `1.8e298` is tiny, and adding it can overflow. The corpus line `newton step that converges into an overflow` solves `x - max_float = 0` from `max_float - 1e297` with a Jacobian of `0.25` instead of `1`: the step is about `4e297`, the test passes and `x + dx` is infinity. So the sum is checked, and a non-finite one is `Error Nan` (without the check the line would print `Ok [inf]`). In the code the last two rows are two `match` cases of `iterate`; the first has a guard (`when`), which matches only if its condition holds. The corpus pins two of the names: a `nan` residual at the start is `Error Nan` (line 39), and a step that the linear solve cannot make finite is `Error Diverged` (line 40: `x - 1` with a Jacobian of `1e-320`, a float that is not 0, whose step `1 / 1e-320` overflows).

Why test the step and not the residual? First, **units**: `dx` is measured in the units of `x`, which is what we want accurate, while the size of `G(x)` depends on how the equation happens to be written (`x² - 2 = 0` and `1e6 (x² - 2) = 0` have the same root and the same Newton steps, with residuals a million times apart). Second, near the root the step estimates the error: since `G(x) ≈ G'(x)(x - r)`, the step `dx = -G'(x)^{-1} G(x)` has about the length of the current error, and taking it leaves an error much smaller. Third, rounding errors in `G` are about `eps` times the size of its terms, so a fixed residual threshold may sit below that noise (never reached) or far above it. The step test also has to come before the line search: close to the root the computed residual stops shrinking (it is rounding noise), so the Armijo test can refuse every damped step and an answer that is already accurate would be reported as `Diverged`.

The reasoning relies on fast convergence. `G(x) = x³` has a triple root at 0, `G'(0)` is `0` too, and the error shrinks by only a factor of two thirds per iteration (linear convergence): from `x_0 = 1` the iterates are `(2/3)^k`, each full step passing the Armijo test. The step at `x` is `x / 3`, and the stopping test needs it below `1e-10`, that is `x` below `3e-10`. Solving `(2/3)^k = 3e-10` gives `k ≈ 54`, so the first iterate below that is `(2/3)^55`, after 55 iterations, more than the limit of 50: with `f x = [| x.(0) *. x.(0) *. x.(0) |]` and `jac x = [| [| 3. *. x.(0) *. x.(0) |] |]`, `Newton.solve f jac [| 1. |]` returns `Error Diverged`. With a limit of 100 it would return `(2/3)^56 ≈ 1.4e-10`, twice the last step from the root: with linear convergence the step understates the remaining error. This is corpus line 46, which prints `Ok [0.000000000138]` with the limit of 100.

## 7. Newton inside the integrator: the stage equation

Every step of the implicit methods (backward Euler and BDF2, members of the BDF family of [chapter 4](04-bdf.md)) reduces to `x = psi + gamma f(t_{n+1}, x)`: backward Euler has `psi = y_n` and `gamma = h`; BDF2 has `gamma = beta h` and a `psi` built from two old states. `Stage.solve` hands Newton the residual and its Jacobian:

```
G(x)  = x - psi - gamma f(t_{n+1}, x)
G'(x) = I - gamma J          J = Jacobian of f, by forward differences (chapter 3)
```

**Why `I - gamma J` and not just `I`?** With `I` in place of `G'`, Newton's update `x - G(x)` would be `x ← psi + gamma f(t_{n+1}, x)`, plain fixed-point iteration. For `y' = -λ y` (this `λ` is the decay rate of chapter 1, not the damping factor of section 2) it multiplies the error by `-gamma λ` each time, so it converges only when `gamma λ < 1`: for backward Euler on the canary's fastest component (`λ = 1e4`, `gamma = h`) only for `h < 1e-4`, the same order as explicit Euler's stability limit `2e-4`. The `gamma J` term is what lets Newton take the long steps that stiff problems allow.

**Few iterations.** Only the term `gamma f` can be nonlinear, and `gamma` is small when the step is, so `G` is almost the linear function `x - psi`. In the error formula of section 1 the constant is `G'' / (2 G') = -gamma f'' / (2 (1 - gamma f'))` in the scalar case, which shrinks with `gamma`. The start is close too: `y_n` for backward Euler, a line through the last two points for BDF2. For `y' = -λ y` the stage equation is exactly linear, `G(x) = (1 + gamma λ) x - psi`, so one Newton step lands on `x = psi / (1 + gamma λ)`, for backward Euler `y_n / (1 + h λ)`, the closed form of [chapter 1](01-odes-and-stiffness.md). The Jacobian is accurate only to about eight digits, so the first step lands very close to the root rather than on it, and one or two more iterations confirm that the step is below the tolerance.

**Large steps are harder.** As `gamma` grows the constant above grows, Newton needs more iterations, and from a poor guess it can fail. `Adaptive.integrate` treats a failed solve as a rejection and halves the step ([chapter 5](05-step-control.md)); `Stepper.fixed` returns the error. Here is one backward Euler step of the logistic equation from `y_n = 0.5` (needs the `problems.ml` link of the probe project). `Stage` is internal to the solver and out of a probe's reach, so the probe calls the public `Bdf1.step`, which sets `psi = y_n` and `gamma = h` and hands the equation to `Stage.solve`, started at `y_n`:

```ocaml
open Vstiff

let () =
  List.iter
    (fun h ->
      match Bdf1.step Problems.Logistic.rhs h Bdf1.start { Ode.t = 0.; y = [| 0.5 |] } with
      | Ok (x, _) -> Printf.printf "h = %-6g Ok %.6f\n" h x.(0)
      | Error e -> Printf.printf "h = %-6g Error %s\n" h (Fail.to_string e))
    [ 1e-3; 1e-1; 10.; 1e4 ]
```

```text
h = 0.001  Ok 0.500250
h = 0.1    Ok 0.524938
h = 10     Ok 0.952494
h = 10000  Error Diverged
```

At `h = 1e4` the stage equation has real roots, near `1` and near `-5e-5`, yet Newton fails: at `x = 0.5` the logistic slope `1 - 2x` is zero, so `G' = 1`, the step is `+2500`, and no factor down to `1/1024` brings the residual below its starting size of `2500`.

**A wrong Jacobian changes the speed, not the equation.** Newton with a matrix `M` in place of the exact Jacobian iterates `x ← x - M⁻¹ G(x)`. For any invertible `M` this stands still only where `G(x) = 0`, so a wrong `M` cannot make the iteration settle on a non-root (apart from the stopping test's blind spot: an absurdly large `M` makes every step tiny). The speed changes. With the exact Jacobian the error is squared at each iteration; with a slightly wrong `M` it shrinks by a constant factor (linear convergence); with a plainly wrong one it grows or the step points uphill and damping gives up. Try a constant slope `2` for `x² - 2` (the true slope is `2x`). `calls` is a mutable counter: `incr` adds one to it and `!calls` reads it. The library's numerical code never does this ([docs/architecture.md](../architecture.md) says why), but a throwaway probe may:

```ocaml
open Numerics

let () =
  let f x = [| (x.(0) *. x.(0)) -. 2. |] in
  let calls = ref 0 in
  let jac _ =
    incr calls;
    [| [| 2. |] |]
  in
  match Newton.solve f jac [| 1. |] with
  | Ok x -> Printf.printf "x = %.12f after %d Jacobian calls\n" x.(0) !calls
  | Error e -> print_endline (Fail.to_string e)
```

```text
x = 1.414213562433 after 25 Jacobian calls
```

The iteration is `x ← x - (x² - 2)/2`, whose slope at the root is `1 - √2 ≈ -0.414`, so each iteration multiplies the error by about `0.414`; since `0.414^25 ≈ 3e-10`, about 25 iterations are needed against 5 with the true slope. The answer is within the tolerance of `√2` (`1.414213562433` against `1.414213562373`).

Two consequences. First, production codes build the Jacobian once and reuse it for several iterations or steps (an old Jacobian is a slightly wrong one); vstiff does not: `Stage.solve` calls `Jac.forward` in every iteration, at `n + 1` evaluations of the right-hand side each time. That is a known limitation. Second, a Jacobian bug can hide: the integrator still returns right answers, only slower or after more rejections. The corpus shows it: with two rows of the Jacobian swapped, the backward Euler canary still finishes within its error bound and only its printed digits move, while a stiffer step makes Newton fail ([chapter 6](06-the-corpus.md)).

## In the code

| Idea | Where |
|------|-------|
| The Newton loop; the exits of section 6, in that order | `Newton.solve`, local function `iterate`, in [`src/numerics/newton.ml`](../../src/numerics/newton.ml) |
| Damping and the Armijo test | local function `damp` (it halves `lambda`); the private constants `tol`, `max_iter`, `min_damping` and `armijo` at the top of the file |
| Gaussian elimination with partial pivoting | `Linalg.solve` and the recursion `solve_augmented`, with helpers `pivot` and `swapped`, in [`src/numerics/linalg.ml`](../../src/numerics/linalg.ml) |
| The residual `G` and the Jacobian `I - gamma J` | `Stage.solve` in [`src/stage.ml`](../../src/stage.ml): local `residual` and `jacobian`, passed to `Newton.solve` |
| Failure values | `Fail.t` in [`src/numerics/fail.ml`](../../src/numerics/fail.ml): `Diverged` and `Nan` (`StepRejected` belongs to chapter 5) |
| Vector helpers | `Vec.axpy a x y` is `a x + y`, `Vec.norm_inf` the largest absolute component, `Vec.finite` true when no entry is NaN or infinite ([`src/numerics/vec.ml`](../../src/numerics/vec.ml)) |

Changing this code needs care because the corpus pins most of it, not all: damping (line 38), its depth (line 45: `x² - 2` from `1e-3` needs the factor `1/512`), the limit of 50 iterations (line 46: `x³` from 1), the refusal of a trial point outside the domain (line 47: `log x` from 3), the names of two failures (lines 39 and 40) and the check on a converged `x + dx` (line 24). The Armijo constant is not pinned, and the finiteness checks on the step and on the trial residual change nothing when removed ([chapter 6](06-the-corpus.md), "What the corpus does not catch"). A behaviour change in a part the corpus does not pin comes with a case that tells the old code from the new ([AGENTS.md](../../AGENTS.md), H3).

## Check yourself

1. **Do one Newton step by hand on `G(x) = x² - 2` from `x_0 = 3`. How does the new error compare with the old one?**
   `x_1 = 3 - 7/6 = 11/6 ≈ 1.833`. The errors are `e_0 = 3 - √2 ≈ 1.586` and `e_1 ≈ 0.419`, and `e_1 / e_0² ≈ 0.167 = 1 / (2 x_0)`, as the identity `e_{k+1} = e_k² / (2 x_k)` says.

2. **Solve `2u + v = 5`, `4u + 3v = 11` by hand, as `Linalg.solve` would. Where is the row exchange?**
   The first column is `(2, 4)`, so the pivot is `4` in row 1: swap. The multiplier is `2/4 = 1/2`. Subtracting half of `4u + 3v = 11` from `2u + v = 5` leaves `-0.5 v = -0.5`, so `v = 1`, and then `4u + 3 = 11` gives `u = 2`.

3. **Why does the stopping test look at the step and not at the residual? When can it mislead?**
   The step is in the units of `x`, does not depend on how the equation is scaled and, under quadratic convergence, is much larger than the error left after taking it. Under linear convergence (a wrong Jacobian, a multiple root) it underestimates the remaining error.

4. **Why can a wrong Jacobian not make Newton settle on a point that is not a root? What can it do instead?**
   The iteration `x ← x - M⁻¹ G(x)` stands still only where `G(x) = 0`. A wrong `M` makes convergence slow or makes the solve fail (`Error Diverged`, through damping or the iteration limit); the adaptive integrator then rejects steps and takes more of them. The one caveat is the stopping test, which trusts small steps (question 3).

5. **The stopping test calls a step of `4e297` tiny when `x` is near `max_float`. What can still go wrong, and what does `Newton.solve` return?**
   The test is relative, `1e-10 (1 + |x|_inf)`, about `1.8e298` there, so the step passes; but `x + dx` can exceed the largest float and become infinity. `Newton.solve` checks the sum and returns `Error Nan`, not `Ok [inf]`.

Next: [3. Jacobians and floating point](03-jacobians-and-floating-point.md), where the Jacobian that Newton needs comes from, and why a computer cannot compute it exactly.
