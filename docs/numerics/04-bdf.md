# 4. BDF methods

Chapter 4 of 6 in the numerical-methods track. Previous: [3. Jacobians and floating point](03-jacobians-and-floating-point.md). Next: [5. Error estimates and step-size control](05-step-control.md). Index and reading order: [docs/README.md](../README.md). Terms and symbols: [glossary](../glossary.md).

**Summary.** A backward differentiation formula (BDF) turns the ODE into one implicit equation per step by differentiating a polynomial through the newest points. BDF1 is backward Euler. BDF2 uses three points on a grid whose steps may differ in size; its three weights depend only on the step ratio $\omega = h / h_{\mathrm{prev}}$, and the code computes them in `Bdf2.coeffs`. This chapter derives the conditions that fix the weights (without writing the weights down: the code is their single source of truth), explains why the order is two, why the first step is a backward Euler step and what the method remembers between steps, how every step reduces to the same stage equation for Newton's method, and the stability facts that limit $\omega$.

## 1. The idea

Notation: $h = t_{n+1} - t_n$ is the step being taken and $h_{\mathrm{prev}} = t_n - t_{n-1}$ the one before it. We know $y_n$, $y_{n-1}$, $\ldots$ and want $y_{n+1}$.

$$
\begin{array}{ccccc}
t_{n-1} & \xrightarrow{h_{\mathrm{prev}}} & t_n & \xrightarrow{\quad h \quad} & t_{n+1} \\
y_{n-1} & & y_n & & y_{n+1} \\
\text{(known)} & & \text{(known)} & & \text{(unknown)}
\end{array}
$$

To approximate a derivative from data, fit a polynomial through the points and differentiate the polynomial. BDF does exactly that, with a twist: the new point $(t_{n+1}, y_{n+1})$ is one of the points, even though $y_{n+1}$ is unknown.

1. Let $q$ be the polynomial of degree $k$ through the $k + 1$ points $(t_{n+1-k}, y_{n+1-k}), \ldots, (t_n, y_n), (t_{n+1}, y_{n+1})$.
2. Demand that $q$ satisfy the ODE at the new time: $q'(t_{n+1}) = f(t_{n+1}, y_{n+1})$.

The polynomial depends linearly on the unknown $y_{n+1}$, so step 2 is one equation for $y_{n+1}$. It is generally nonlinear, because $f$ is, and $y_{n+1}$ appears inside $f$. The method is implicit, and implicit methods can have the stability that stiff problems need (chapter 1). The name: it is a differentiation formula, since it approximates a derivative, and it looks only backwards in time. The member that uses degree $k$ is **$\mathrm{BDF}k$**. A multistep method such as BDF gets a higher order by remembering earlier steps, and each step still costs one implicit solve of the same shape whatever the order.

One consequence of this construction deserves a name: if the true solution $y(t)$ is a polynomial of degree at most $k$, then $q = y$ and the $\mathrm{BDF}k$ formula is exact. That "exact on polynomials up to degree $k$" property is how we derive the formula and why its order is $k$.

## 2. BDF1 is backward Euler

With $k = 1$, $q$ is the line through $(t_n, y_n)$ and $(t_{n+1}, y_{n+1})$. Its slope is $(y_{n+1} - y_n) / h$, so $q'(t_{n+1}) = f(t_{n+1}, y_{n+1})$ reads

$$
y_{n+1} = y_n + h\thinspace f(t_{n+1}, y_{n+1})
$$

which is backward Euler, as in chapter 1. It is exact for linear functions and not for quadratics: for $y = t^2$ (so $f = 2t$) one step from exact data overshoots by exactly $h^2$. That is $\frac{h^2}{2} y''$ with $y'' = 2$, the size of the local error in chapter 1 with the opposite sign: the local error (exact minus computed) of backward Euler is $-\frac{h^2}{2} y''$. Exact for degree 1 and wrong at degree 2 means order 1.

## 3. BDF2 on a variable grid

With $k = 2$, $q$ is the parabola through $(t_{n-1}, y_{n-1})$, $(t_n, y_n)$ and $(t_{n+1}, y_{n+1})$. Rearranging $q'(t_{n+1}) = f(t_{n+1}, y_{n+1})$ to put $y_{n+1}$ on the left gives a formula of this shape (with $y_{n+1}$ still inside $f$), the one `Bdf2` implements ([`src/bdf2.mli`](../../src/bdf2.mli) documents it):

$$
y_{n+1} = a_1 y_n + a_0 y_{n-1} + \beta h\thinspace f(t_{n+1}, y_{n+1})
$$

The three numbers $a_1$, $a_0$ and $\beta$ depend only on the step ratio $\omega = h / h_{\mathrm{prev}}$: the conditions below involve the time differences only as ratios, so scaling all steps together changes nothing. The code computes them in `Bdf2.coeffs omega`.

**Deriving them.** A formula of this shape that is exact for every quadratic must reproduce $y = 1$, $y = s$ and $y = s^2$, because any quadratic is a combination of those three and the formula is linear in $y$. Use shifted time $s = t - t_{n+1}$, so the new point is at $s = 0$, the previous point at $s = -h$, and the one before it at $s = -(h + h_{\mathrm{prev}}) = -r h$, where

$$
r = \frac{h + h_{\mathrm{prev}}}{h} = 1 + \frac{1}{\omega}
$$

is the distance from $t_{n+1}$ back to $t_{n-1}$ measured in steps $h$. Write the formula as $y(0) = a_1 y(-h) + a_0 y(-r h) + \beta h\thinspace y'(0)$, insert each test function and divide by $h$ or $h^2$:

$$
\begin{aligned}
y = 1: &\quad 1 = a_1 + a_0 \\
y = s: &\quad 0 = -a_1 - a_0 r + \beta \qquad (\text{here } y' = 1) \\
y = s^2: &\quad 0 = a_1 + a_0 r^2 \qquad (\text{here } y'(0) = 0 \text{, so } \beta \text{ drops out})
\end{aligned}
$$

Three linear equations in three unknowns; $\omega$ enters only through $r$. This page stops here, on purpose: `Bdf2.coeffs` is the single source of truth for the numbers, and a copy of its formulas or output in a document could go stale. Finish the algebra yourself (exercise 1) or ask the code.

**Probe: check the code against the conditions.** The residual of the formula on exact data is the exact $y(t_{n+1})$ minus the formula's right-hand side. It must vanish for $1$, $t$ and $t^2$ (so for every quadratic) and need not for $t^3$. `Bdf2.coeffs` supplies the weights. The probes of this chapter run in the probe project of [chapter 1](01-odes-and-stiffness.md) ("Running the snippets"):

```ocaml
open Vstiff

(* Residual of the BDF2 formula on exact data: y at the new time t, minus the formula's
   right-hand side. The previous step is h / omega. *)
let residual ~y ~dy ~t ~h ~omega =
  let { Bdf2.a1; a0; beta } = Bdf2.coeffs omega in
  let h_prev = h /. omega in
  y t -. ((a1 *. y (t -. h)) +. (a0 *. y (t -. h -. h_prev)) +. (beta *. h *. dy t))

let one _ = 1. and one' _ = 0.
let line t = t and line' _ = 1.
let square t = t *. t and square' t = 2. *. t
let cube t = t *. t *. t and cube' t = 3. *. t *. t

let () =
  let tests = [ ("1", one, one'); ("t", line, line'); ("t^2", square, square'); ("t^3", cube, cube') ] in
  List.iter
    (fun omega ->
      List.iter
        (fun (name, y, dy) ->
          Printf.printf "omega = %-4g y = %-3s residual = %9.2e\n" omega name (residual ~y ~dy ~t:1.5 ~h:0.5 ~omega))
        tests)
    [ 0.25; 1.; 2. ];
  (* The cubic's residual is a fixed multiple of h^3. *)
  List.iter
    (fun h -> Printf.printf "t^3, h = %-6g residual / h^3 = %.4f\n" h (residual ~y:cube ~dy:cube' ~t:1.5 ~h ~omega:1. /. (h ** 3.)))
    [ 0.1; 0.05; 0.025 ]
```

OCaml notes ([docs/ocaml.md](../ocaml.md) has the details): `let { Bdf2.a1; a0; beta } = Bdf2.coeffs omega in` unpacks the record that `coeffs` returns into three variables named after its fields (the `Bdf2.` prefix says which module the field names belong to); `~y`, `~dy`, `~t`, `~h` and `~omega` are labelled arguments, like Python keyword arguments (here they stop you from swapping two floats); the `'` in `square'` is part of the name, a habit for derivatives; `List.iter f [ ... ]` runs `f` once for each element of the list; `**` is the power operator on floats; `%-4g` and `%9.2e` are `Printf` formats for a float.

The residuals for $1$, $t$ and $t^2$ are zero or at rounding level, a small multiple of $10^{-16}$; the one for $t^3$ is not, and divided by $h^3$ it is the same number in every row (the last three lines), as an error of size $h^3$ should be. At $\omega = 1$, equal steps, the formula reduces to the constant-step BDF2 formula found in any numerical analysis text. The first condition, $a_1 + a_0 = 1$, says the weights on the two old values add up to one. That has a practical use: the three right-hand sides of Robertson's problem sum to zero, so a BDF step turns a total of 1 into $a_1 + a_0$ times 1, and $y_1 + y_2 + y_3$ stays at 1 ([chapter 6](06-the-corpus.md) checks it to $10^{-8}$).

**Exercises.**

1. Solve the three equations for $a_1$, $a_0$ and $\beta$ (hint: the first gives $a_1$ in terms of $a_0$). Compare your results with the fields of `Bdf2.coeffs omega` at $\omega = 1$ and $\omega = 2$ (print them with `Printf.printf "%g %g %g\n" a1 a0 beta`), and with the constant-step formula in a textbook.
2. Derive the same numbers the other way. Write $q(t) = y_{n+1} \ell_0(t) + y_n \ell_1(t) + y_{n-1} \ell_2(t)$, where $\ell_0, \ell_1, \ell_2$ are the quadratic Lagrange basis polynomials for the nodes $t_{n+1}, t_n, t_{n-1}$ (each is 1 at its own node and 0 at the other two). Differentiate at $t_{n+1}$, set the result equal to $f(t_{n+1}, y_{n+1})$ and solve for $y_{n+1}$. Dividing by the coefficient of $y_{n+1}$ must give your $a_1$, $a_0$ and $\beta h$.
3. What does `Bdf2.coeffs` return as $\omega$ tends to 0? Evaluate it at $\omega = 10^{-3}$ and $\omega = 10^{-9}$. Which familiar method does the formula approach, and why does that make sense when $h$ is tiny compared with $h_{\mathrm{prev}}$?

## 4. Order two

The formula reproduces everything up to $s^2$ exactly and misses the $s^3$ term of the Taylor expansion of $y$. Applying it to $y = s^3$ leaves a residual that is a multiple of $h^3$, with a factor depending on $\omega$ only (the cubic probe above). So, to leading order,

$$
\text{local error of BDF2} = C(\omega) h^3 y'''(t) + O(h^4)
$$

Local error $O(h^3)$ means global error $O(h^2)$ when the method is stable (section 7), so BDF2 has order 2.

**How the corpus measures it.** For a method of order $p$, error $\approx C h^p$, so halving $h$ divides the error by $2^p$. The corpus integrates the logistic equation on $[0, 5]$ with `Stepper.fixed (module Bdf2)` at `dt = 0.01` and `dt = 0.005` and compares the errors at $t = 5$. The expected output file records:

```
bdf2 logistic [0,5]: error dt=0.01 3.116e-06, dt=0.005 7.810e-07, ratio 3.99 in [3.5, 4.5]: true
```

A ratio of 3.99 means $p = \log_2 3.99 \approx 2.00$. The test accepts ratios from 3.5 to 4.5, orders from about 1.81 to 2.17, so first order (ratio 2) and third order (ratio 8) both fail it. Why the logistic equation? It is smooth and not stiff, so truncation error dominates; it has a closed-form solution to compare against; and its right-hand side is nonlinear, so Newton genuinely runs. The errors are thousands of times larger than Newton's tolerance of $10^{-10}$, so the ratio shows the truncation error of the method rather than the accuracy of the nonlinear solves.

You can reproduce the two numbers, and add a third step size, with:

```ocaml
open Vstiff
open Numerics

let () =
  let open Problems.Logistic in
  let error dt =
    match Stepper.fixed (module Bdf2) ~dt problem with
    | Ok y -> Vec.norm_inf (Vec.sub y (exact 5.))
    | Error e -> failwith (Fail.to_string e)
  in
  List.iter (fun dt -> Printf.printf "dt = %-7g error = %.3e\n" dt (error dt)) [ 0.01; 0.005; 0.0025 ]
```

`(module Bdf2)` hands the module to `Stepper.fixed` ([docs/ocaml.md](../ocaml.md) explains modular explicits); `Vec.norm_inf (Vec.sub y (exact 5.))` is the largest absolute entry of the difference between the computed and the exact state; `failwith` raises an exception, which is fine in a throwaway probe (library code returns `Error` instead). The first two lines match the corpus line above. Check that the third error is again about a quarter of the second.

## 5. The first step and the history

BDF2 needs $y_{n-1}$. At $n = 0$ we only have $y_0$. The method therefore carries a small memory from step to step, its **history** (the `history` type of the `Ode.Method` contract; [docs/architecture.md](../architecture.md) describes the contracts). Inside [`src/bdf2.ml`](../../src/bdf2.ml) it is a variant type with two cases, one carrying a record ([docs/ocaml.md](../ocaml.md) covers the syntax):

```text
type history = Start | After of { h_prev : float; y_prev : Vec.t }
```

`Start` means there is no previous point, so the step is a backward Euler step (`Bdf1.step`). `After { h_prev; y_prev }` holds the previous step and the state it started from, which is exactly what a BDF2 step needs. Every step returns the next history: `After` with `h_prev` equal to the step $h$ just taken and `y_prev` equal to the state $y_n$ it started from. The interface hides the constructors: callers get a history from `Bdf2.start` or from an earlier `step`, and pass it back. This probe shows the first step is exactly backward Euler and the second is not:

```ocaml
open Vstiff

let () =
  let open Problems.Logistic in
  let ok = function Ok v -> v | Error e -> failwith (Fail.to_string e) in
  let h = 0.1 in
  let at0 = { Ode.t = 0.; y = y0 } in
  let y1_be, _ = ok (Bdf1.step rhs h Bdf1.start at0) in
  let y1, history = ok (Bdf2.step rhs h Bdf2.start at0) in
  Printf.printf "first step identical to backward Euler: %b\n" (y1 = y1_be);
  let at1 = { Ode.t = h; y = y1 } in
  let y2_be, _ = ok (Bdf1.step rhs h Bdf1.start at1) in
  let y2, _ = ok (Bdf2.step rhs h history at1) in
  Printf.printf "second step differs from backward Euler: %b\n" (y2 <> y2_be)
```

Whoever drives the method chooses which history to pass. `Stepper.fixed` passes each step's history to the next. The adaptive driver passes the history of the last *accepted* step again after a rejection, so a retry uses the same previous point ([chapter 5](05-step-control.md)).

The first step is only first-order accurate: its local error is $O(h^2)$, one power of $h$ worse than BDF2's $O(h^3)$. That does not lower the order of the whole run. A zero-stable method (section 7) carries an early error forward without amplifying it by more than a bounded factor, so the single $O(h^2)$ error from the start adds an $O(h^2)$ term to the final error. The remaining $T / h$ steps each commit $O(h^3)$, which adds up to $O(h^2)$ as well. Everything has the same size, $h^2$, so the global error is still $O(h^2)$. The measured ratio of 3.99 in section 4 includes this startup step.

## 6. The stage equation: one shape for every BDF method

Move everything known to one side. Backward Euler and BDF2 both become

$$
x = \psi + \gamma\thinspace f(t_{n+1}, x), \qquad x = y_{n+1}
$$

| Method | $\psi$ | $\gamma$ |
|--------|-----|-------|
| BDF1 (backward Euler) | $y_n$ | $h$ |
| BDF2 | $a_1 y_n + a_0 y_{n-1}$ | $\beta h$ |

$\psi$ collects the history; $\gamma$ is the weight that multiplies the new slope. Newton's method solves the residual form

$$
G(x) = x - \psi - \gamma\thinspace f(t_{n+1}, x) = 0, \qquad \text{Jacobian of } G\text{:}\ I - \gamma J
$$

where $J$ is the Jacobian of $f$ (computed by finite differences, [chapter 3](03-jacobians-and-floating-point.md)). All of that lives in one function, `Stage.solve rhs { Stage.t; gamma; psi } guess` in [`src/stage.ml`](../../src/stage.ml), shared by all methods. Each method contributes only $\psi$, $\gamma$ and an initial `guess`; a higher BDF order would need new weights and a longer history, but the stage solve would not change. On the test equation $y' = -\lambda y$ the stage equation gives $x = \psi / (1 + \gamma\lambda)$: when $\gamma\lambda$ is huge the new value is close to zero whatever $\psi$ says, the same damping as backward Euler's $1 / (1 + h\lambda)$ in chapter 1.

**The initial guess.** Newton converges fast when it starts near the root. For backward Euler the guess is $y_n$ itself, which is $O(h)$ away from $y_{n+1}$. For BDF2 the code starts from the straight line through the last two points, extended to $t_{n+1}$:

$$
\mathrm{guess} = y_n + \omega (y_n - y_{n-1})
$$

(the line has slope $(y_n - y_{n-1}) / h_{\mathrm{prev}}$, and we move along it for time $h = \omega h_{\mathrm{prev}}$). For a smooth solution this guess is $O(h^2)$ away from $y_{n+1}$, so Newton starts closer to its root.

## 7. Stability

**A-stability.** Chapter 1 defined A-stability informally: the numerical solution of $y' = \lambda y$ decays for every positive $h$ whenever $\lambda$ has negative real part. Backward Euler is A-stable, and so is BDF2 with constant steps (a classical result). With variable steps the A-stability picture is more involved and this chapter does not claim it. What the code does rely on is zero-stability, below. The corpus pins the stiff behaviour with this line:

```
bdf2 canary t=1 dt=1e-3 (h lambda = 10): max error 1.53e-07 < 1e-06: true
```

The canary's fast component has $h\lambda = 10$, where explicit Euler would explode (chapter 1), and BDF2 ends with max error $1.53 \times 10^{-7}$.

Higher orders pay a price. BDF3 to BDF6 are not A-stable; they are only $\mathrm{A}(\alpha)$-stable, meaning stable for every positive $h$ when $\lambda$ lies within an angle $\alpha$ of the negative real axis, and $\alpha$ shrinks as the order grows. From BDF7 on the methods are not even zero-stable, so they do not converge. Dahlquist's second barrier says that no A-stable linear multistep method has order above 2. So order 2 is the highest order with full A-stability, and BDF2 is the natural stopping point for a small library ([chapter 5](05-step-control.md) mentions what production codes do instead).

**Zero-stability.** A multistep method has extra solutions that have nothing to do with the ODE, and they must die out. Take the simplest ODE, $y' = 0$, whose solutions are constants. Then $f = 0$ and a BDF2 step is $y_{n+1} = a_1 y_n + a_0 y_{n-1}$. This is a recurrence: a rule that computes each value from the two before it. With the coefficients held fixed, try $y_n = \zeta^n$. It satisfies the rule exactly when $\zeta^2 = a_1 \zeta + a_0$, and every solution is a combination of solutions of that form (when the two roots differ).

- Because $a_1 + a_0 = 1$, $\zeta = 1$ is always a root. It carries the true constant solution.
- The two roots multiply to $-a_0$, so the other root, the **parasitic root**, is $-a_0$.
- If its magnitude is below 1 the spurious solution fades. If it is above 1 it grows geometrically, however small $h$ is. A method with that property does not converge as $h \to 0$.

Requiring all roots to have magnitude at most 1 (and those of magnitude exactly 1 to be simple) is called **zero-stability**.

With variable steps the coefficients change at every step, so the picture with frozen coefficients is only the starting point. The classical result is that variable-step BDF2 is zero-stable as long as every step ratio stays below $1 + \sqrt{2}$. Above that, a constant ratio makes the parasitic root exceed 1, and the method is no longer guaranteed to be stable.

**Why the controller's $\omega \le 2$ is safe.** `Bdf2` does not check $\omega$; it computes $h / h_{\mathrm{prev}}$ and uses it. The guarantee comes from how the adaptive driver and the `Halving` controller work together ([chapter 5](05-step-control.md) has the rules in full). The history the driver passes to the method always describes the last *accepted* step, so $\omega$ is the new attempt over the last accepted step, and the controller can only change the step in these ways:

- After an accepted step it keeps its proposal or doubles it (at most once per three accepts). The step just taken was the previous proposal, so the next attempt is at most twice as long. (The one exception is the final step, which takes whatever is left to reach $t_{\mathrm{end}}$, and then the run ends.)
- After a rejection it halves the step that failed, and the history stays what it was. The failed step was at most twice the last accepted step, so the retry is at most as long as that: $\omega \le 1$.
- A fixed-step run, `Stepper.fixed (module Bdf2)`, has $\omega = 1$ on every BDF2 step up to rounding: each step is the difference of two grid times, which is $h$ to within an ulp of $t$ (its first step is backward Euler).

So every ratio of proposals is at most 2, safely below $1 + \sqrt{2} \approx 2.414$. (The driver also snaps each step to the floats. It rejects a non-final step below `Clock.resolution t`, $16\thinspace\varepsilon\thinspace\lvert t \rvert$, and above that the snap changes a step by at most $1/32$ of its length, so each step is within $1/32$ of its proposal and $\omega$ stays at most $2 \cdot \frac{33/32}{31/32} = 66/31$, about 2.13 ([chapter 5](05-step-control.md), section 9). Without the rejection a step of a few ulps could snap to a very different length and pass the limit.) Any change to step control that lets the ratio exceed 2 must re-check it against $1 + \sqrt{2}$ first, and the check belongs with the controller, because a different `Ode.Controller` can propose anything. Very small ratios are harmless: as $\omega \to 0$ BDF2 approaches backward Euler (exercise 3).

**Probe: watch the ratios.** This wrapper is an `Ode.Embedded` method made of `Bdf2` plus a note of the previous step in its history, and it records the largest $h / h_{\mathrm{prev}}$ it is ever asked to take. The driver drops the history of a rejected attempt, so the note always describes the last accepted step:

```ocaml
open Vstiff

let worst = ref 0.

module Spy : Ode.Embedded = struct
  type history = Bdf2.history * float

  let start = (Bdf2.start, 0.)
  let note h h_prev = if h_prev > 0. then worst := Float.max !worst (h /. h_prev)

  let step rhs h (history, h_prev) at =
    note h h_prev;
    Result.map (fun (y, next) -> (y, (next, h))) (Bdf2.step rhs h history at)

  let step_with_error rhs h (history, h_prev) at =
    note h h_prev;
    Result.map (fun (y, err, next) -> (y, err, (next, h))) (Bdf2.step_with_error rhs h history at)
end

let () =
  List.iter
    (fun (name, tol, problem) ->
      worst := 0.;
      let outcome = match Adaptive.integrate (module Spy) (module Halving) ~tol problem with Ok _ -> "Ok" | Error e -> Fail.to_string e in
      Printf.printf "%-12s %-4s largest ratio %g\n" name outcome !worst)
    [ ("logistic", 1e-6, Problems.Logistic.problem); ("van der Pol", 1e-4, Problems.VanDerPol.problem); ("Robertson", 1e-6, Problems.Robertson.problem) ]
```

OCaml notes: `module Spy : Ode.Embedded = struct ... end` is a module that satisfies the contract by delegating to `Bdf2`; `type history = Bdf2.history * float` pairs the real history with a number; `ref`, `:=` and `!` make, set and read a mutable cell, which is fine in a probe (the library keeps mutable state only in `Instrument`). No line should print a ratio above 2. If you write another controller, run it through the same wrapper before trusting it.

**Exercise 4.** Use your $a_0$ from exercise 1 to find the $\omega$ at which the parasitic root $-a_0$ reaches 1; you should get $1 + \sqrt{2}$. Then check it with the code: evaluate $-a_0$ from `Bdf2.coeffs` at $\omega = 2.4$ and $\omega = 2.5$, just below and above the crossing.

## In the code

| Idea | Where |
|------|-------|
| The weights $a_1$, $a_0$, $\beta$ for a given $\omega$ | `Bdf2.coeffs omega`, returning a record `{ a1; a0; beta }`, in [`src/bdf2.ml`](../../src/bdf2.ml). The only place where the numbers are defined |
| One variable-step BDF2 step | `Bdf2.step rhs h history at`, which returns the new state and the next history. Once a previous point exists the internal `bdf2` does the work: the ratio $\omega = h / h_{\mathrm{prev}}$, the weights from `coeffs`, $\psi$ and the initial guess of section 6 built with the helpers of [`src/numerics/vec.ml`](../../src/numerics/vec.ml) (each returns a fresh array; the library never writes into an array it was given), then `Stage.solve` at $t_n + h$ with $\gamma = \beta h$ |
| Backward Euler (BDF1), and BDF2's first step | `Bdf1.step` in [`src/bdf1.ml`](../../src/bdf1.ml) |
| The history | `Bdf2.start` is `Start`; every step returns `After { h_prev; y_prev }` (constructors internal to [`src/bdf2.ml`](../../src/bdf2.ml)) |
| The shared stage equation | `Stage.solve rhs { Stage.t; gamma; psi } guess` in [`src/stage.ml`](../../src/stage.ml) |
| Step ratios in the adaptive integrator | `Adaptive.integrate` ([`src/adaptive.ml`](../../src/adaptive.ml)) hands the next history back only after an accept; the halve and double rules are in `Halving` ([`src/halving.ml`](../../src/halving.ml)) |
| The order test | The `bdf2 logistic` line of [`test/corpus.ml`](../../test/corpus.ml); the logistic problem in [`test/problems/problems.ml`](../../test/problems/problems.ml); the expected line in [`test/corpus.expected`](../../test/corpus.expected) |
| The same test, rerun ten times | `LogisticOrder` in [`test/soak.ml`](../../test/soak.ml) |

## Check yourself

1. **Which polynomial does BDF2 differentiate, and which equation does it enforce?**
   The parabola through $(t_{n-1}, y_{n-1})$, $(t_n, y_n)$ and the unknown $(t_{n+1}, y_{n+1})$; it enforces $q'(t_{n+1}) = f(t_{n+1}, y_{n+1})$.

2. **How many conditions fix $a_1$, $a_0$ and $\beta$, and why does only $\omega$ appear?**
   Three: exactness on $1$, $s$ and $s^2$. After dividing by $h$ (or $h^2$) the time differences appear only through the ratio $r = 1 + 1/\omega$.

3. **On an ODE whose exact solution is a quadratic, what does the BDF2 formula give from exact data, and what do you expect on a cubic?**
   The exact new value up to rounding, for any $\omega$. On a cubic the residual is a multiple of $h^3$, so halving $h$ (with $\omega$ fixed) divides it by 8.

4. **Why can't the first step use BDF2? What does the code do instead, and why does the order survive?**
   There is no $y_{-1}$. The history is `Start`, so the step is backward Euler. Its $O(h^2)$ error is carried forward without growth and is the same size as the accumulated $O(h^2)$ of all later steps.

5. **State the BDF2 stage equation using $\psi$ and $\gamma$. What is the Jacobian Newton uses?**
   $x = \psi + \gamma\thinspace f(t_{n+1}, x)$ with $\psi = a_1 y_n + a_0 y_{n-1}$ and $\gamma = \beta h$. Newton solves $G(x) = x - \psi - \gamma\thinspace f(t_{n+1}, x) = 0$ with Jacobian $I - \gamma J$.

6. **What is Newton's initial guess for a BDF2 step, and how far is it from the answer?**
   The line through the last two points, $y_n + \omega (y_n - y_{n-1})$, which is $O(h^2)$ away for smooth solutions. Backward Euler starts from $y_n$, which is $O(h)$ away.

7. **What is the largest step ratio the adaptive integrator can produce, and why is that safe?**
   2 for the proposals, below the zero-stability limit $1 + \sqrt{2}$. `Halving` doubles at most once per three accepts and halves on every rejection, so the ratio of its proposals never exceeds 2. Snapping the steps to the floats changes each by at most $1/32$ of its length (the driver rejects non-final steps below $16\thinspace\varepsilon\thinspace\lvert t \rvert$, at least 16 ulps), which lifts the ratio of the steps themselves to at most $66/31$, about 2.13 (chapter 5, section 9). A controller with a larger growth factor would need the bound re-checked.

8. **The corpus reports a ratio of 3.99 for the logistic test. What would a first-order method report, and what does the test check?**
   About 2. The test passes only for ratios between 3.5 and 4.5, i.e. orders from about 1.81 to 2.17.

Next: [5. Error estimates and step-size control](05-step-control.md), how the integrator decides when to accept a step and what step to try next.
