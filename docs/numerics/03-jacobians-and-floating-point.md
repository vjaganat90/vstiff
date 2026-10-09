# 3. Jacobians and floating point

Chapter 3 of 6 in the numerical-methods track. Previous: [2. Newton's method](02-newton.md). Next: [4. BDF methods](04-bdf.md). Index and reading order: [docs/README.md](../README.md). Terms and symbols: [glossary](../glossary.md).

**Summary.** A floating-point number is scientific notation in base 2, so $0.1 + 0.2$ is not $0.3$ and subtracting nearly equal numbers loses digits (cancellation). That limits how well a derivative can be estimated from function values: a forward difference has a truncation error that shrinks with the step and a round-off error that grows as it shrinks, and the best step is near $\sqrt{\varepsilon}$. `Jac.forward` builds the Jacobian that Newton needs from such differences, with the step $10^{-8}(1 + \lvert y_j \rvert)$ and a division by the perturbation that was actually stored. A worked example shows a Jacobian entry of the corpus' canary off by $3 \times 10^{-5}$ although the function is linear, and why the corpus checks the Jacobian absolutely at the origin and relatively elsewhere. The chapter ends with NaN and infinity, and with the fused multiply-add, which can change the last digits of a result when arithmetic is rearranged.

The snippets run in the probe project of [chapter 1](01-odes-and-stiffness.md) ("Running the snippets"; [docs/ocaml.md](../ocaml.md) explains the OCaml). Notation: the Jacobian of the right-hand side $f$ is the matrix $J$ with entries $J_{ij} = \partial f_i / \partial y_j$ (`J.(i).(j)` in the code): row $i$ is the output component, column $j$ the input component. $\varepsilon = 2^{-52}$ is `Float.epsilon`, about $2.2 \times 10^{-16}$. $\lVert v \rVert_\infty$ is the largest absolute component of $v$. In this chapter $\delta$ is the size of a perturbation of an argument (the step of a finite difference); $h$ stays reserved for an integration step, which does not appear here.

## 1. Floating-point numbers

**What a float is.** OCaml's `float` is the 64-bit IEEE 754 number (binary64) that Python calls `float` and JavaScript `number`: a sign, an exponent $e$ and 52 fraction bits $f$, giving $\pm 1.f \times 2^e$ with $1.f$ read in binary. The 53 digits $1.f$ (the leading 1 is implicit) are the **significand**, about 16 decimal digits. Between two consecutive powers of two, $[2^e, 2^{e+1})$, the representable numbers are evenly spaced, $2^{e-52}$ apart: the **ulp** (unit in the last place) of those numbers. For $x$ it lies between $\varepsilon \lvert x \rvert / 2$ and $\varepsilon \lvert x \rvert$, so the *relative* precision is roughly constant while the *absolute* precision grows with the size of the number. $\varepsilon = 2^{-52}$ is the ulp of 1.

**Rounding.** Each $+$, $-$, $\times$, $/$ and $\sqrt{\cdot}$ returns the exact result rounded to the nearest float, so one operation has a relative error of at most the unit roundoff $u = \varepsilon / 2 \approx 1.1 \times 10^{-16}$ (unless the result exceeds the largest float, about $1.8 \times 10^{308}$, or is tinier than about $2.2 \times 10^{-308}$). Nothing is promised about a long computation, where errors can add up. Decimal fractions such as $0.1$ have no finite binary expansion, just as $1/3$ has no finite decimal one, so $0.1$, $0.2$ and $0.3$ are each replaced by the nearest float and the sum is rounded again:

```ocaml
let ulp x = Float.succ x -. x (* gap to the next float above x *)

let () =
  Printf.printf "0.1 + 0.2   = %.17g\n" (0.1 +. 0.2);
  Printf.printf "equals 0.3? %b\n" (0.1 +. 0.2 = 0.3);
  Printf.printf "epsilon     = %.17g\n" Float.epsilon;
  List.iter (fun x -> Printf.printf "ulp(%g) = %.17g\n" x (ulp x)) [ 1.; 1e4; 1e10 ];
  Printf.printf "1e10 + 1e-8 = 1e10? %b\n" (1e10 +. 1e-8 = 1e10);
  Printf.printf "(1 + 1e-15) - 1 = %.17g\n" ((1. +. 1e-15) -. 1.)
```

```text
0.1 + 0.2   = 0.30000000000000004
equals 0.3? false
epsilon     = 2.2204460492503131e-16
ulp(1) = 2.2204460492503131e-16
ulp(10000) = 1.8189894035458565e-12
ulp(1e+10) = 1.9073486328125e-06
1e10 + 1e-8 = 1e10? true
(1 + 1e-15) - 1 = 1.1102230246251565e-15
```

Two lessons. First, comparing floats for equality after arithmetic is rarely what you want, so the code compares with tolerances (Newton's stopping test, the corpus criteria); its few tests for exactly $0$, such as a zero pivot, are fragile too ([chapter 2](02-newton.md), section 5). Second, adding a small number to a large one can do nothing at all: the floats near $10^{10}$ are $1.9 \times 10^{-6}$ apart, so $10^{10} + 10^{-8}$ is $10^{10}$ in floating point. Section 4 uses this.

**Cancellation.** The last line is the third lesson. $(1 + 10^{-15}) - 1$ should be $10^{-15}$ and is $1.11 \times 10^{-15}$, eleven percent off. The sum $1 + 10^{-15}$ had to be rounded to a multiple of $2.2 \times 10^{-16}$ (five of them, $1.11 \times 10^{-15}$), and that rounding error is a tenth of the answer. The subtraction itself is exact, because two floats within a factor of two of each other have a representable difference; it merely exposes the error made when the sum was rounded. A difference of nearly equal numbers has a relative error of about

$$
\frac{\text{absolute error of the operands}}{\lvert \text{difference} \rvert}
$$

This is **cancellation**, and it is what finite differences fight against.

## 2. Derivatives from function values: forward differences

Newton needs $J$, and the library's right-hand sides are black boxes, so $J$ must come from values of $f$. By Taylor, $f(y + \delta) = f(y) + \delta f'(y) + \frac{\delta^2}{2} f''(\xi)$ for some $\xi$ between $y$ and $y + \delta$, so

$$
\frac{f(y + \delta) - f(y)}{\delta} = f'(y) + \frac{\delta}{2} f''(\xi)
$$

The quotient is the derivative plus an error of about $\frac{\delta}{2} \lvert f'' \rvert$, the **truncation error**. It shrinks with $\delta$, so take $\delta$ tiny. But $f(y + \delta)$ and $f(y)$ are nearly equal, so the subtraction cancels: each value carries a relative error up to $u$, the difference is off by up to about $\varepsilon \lvert f \rvert$, and dividing by $\delta$ gives the **round-off error** of about $\varepsilon \lvert f \rvert / \delta$. It grows as $\delta$ shrinks. Their sum, $\frac{\delta}{2} \lvert f'' \rvert + \varepsilon \lvert f \rvert / \delta$, is smallest where the two are about equal:

$$
\delta^* = \sqrt{\frac{2 \varepsilon \lvert f \rvert}{\lvert f'' \rvert}}, \qquad \text{about } \sqrt{\varepsilon} = 1.5 \times 10^{-8} \text{ when } \lvert f \rvert \text{ and } \lvert f'' \rvert \text{ are of order } 1
$$

and the best error is of order $\sqrt{\varepsilon}$ too: about half of the 16 digits survive. This program sweeps the step for $f(y) = \sqrt{y}$ at $y = 2$ ($f' = 0.5 / \sqrt{2}$) and prints the error next to the two estimates, $\frac{\delta}{2} \lvert f'' \rvert$ (truncation) and $\varepsilon \lvert f \rvert / \delta$ (round-off):

```ocaml
let () =
  let f = sqrt and y = 2. in
  let exact = 0.5 /. sqrt y in
  let f2 = 0.25 /. (y *. sqrt y) (* |f''(y)| *) in
  Printf.printf "%-7s %11s %11s %11s\n" "delta" "error" "truncation" "round-off";
  List.iter
    (fun delta ->
      let stored = y +. delta -. y in
      let q = (f (y +. delta) -. f y) /. stored in
      Printf.printf "%-7.0e %11.3e %11.3e %11.3e\n" delta (q -. exact) (0.5 *. delta *. f2)
        (Float.epsilon *. f y /. delta))
    [ 1e-2; 1e-4; 1e-6; 1e-8; 1e-10; 1e-12; 1e-14; 1e-16 ]
```

```text
delta         error  truncation   round-off
1e-02    -4.408e-04   4.419e-04   3.140e-14
1e-04    -4.419e-06   4.419e-06   3.140e-12
1e-06    -4.427e-08   4.419e-08   3.140e-10
1e-08    -1.877e-09   4.419e-10   3.140e-08
1e-10    -6.772e-07   4.419e-12   3.140e-06
1e-12    -8.980e-05   4.419e-14   3.140e-04
1e-14    -5.727e-03   4.419e-16   3.140e-02
1e-16           nan   4.419e-18   3.140e+00
```

For $\delta \ge 10^{-6}$ the error equals the truncation estimate and shrinks a hundredfold per row. Around $\delta = 10^{-8}$ it bottoms out at about $2 \times 10^{-9}$. Below that it grows like $1/\delta$, staying under the $\varepsilon \lvert f \rvert / \delta$ estimate: at $\delta = 10^{-14}$ the derivative is off by 1.6 percent. At $\delta = 10^{-16}$, $2 + 10^{-16}$ is $2$ (the floats near 2 are $4.4 \times 10^{-16}$ apart), the perturbation actually applied is zero, and the quotient is $0/0$, which is `nan`. Smaller is not better. (`stored` is `(y + delta) - y`, the perturbation that was really applied; section 4 explains why.)

## 3. From one derivative to a Jacobian

For $f$ with $n$ components, column $j$ of the Jacobian is approximated by perturbing only $y_j$:

$$
J_{ij} \approx \frac{f_i(y + \delta_j e_j) - f_i(y)}{\delta_j}, \qquad e_j = \text{the } j\text{-th unit vector}
$$

One evaluation $f(y)$ is shared by all columns and each column needs one more, so a Jacobian costs $n + 1$ evaluations of $f$. Newton asks for a new one at every iteration ([chapter 2](02-newton.md)), which is expensive; production codes reuse one across iterations. `Jac.forward f y` returns an array of rows: `Array.length (f y)` rows and `Array.length y` columns, so three outputs and two inputs give a $3 \times 2$ matrix. A check of the orientation, with $f_1 = y_1 + 2 y_2$ and $f_2 = 3 y_1 + 4 y_2$:

```ocaml
open Numerics

let () =
  let f y = [| y.(0) +. (2. *. y.(1)); (3. *. y.(0)) +. (4. *. y.(1)) |] in
  let j = Jac.forward f [| 1.; 1. |] in
  Array.iter
    (fun row ->
      Array.iter (Printf.printf "%8.4f") row;
      print_newline ())
    j
```

```text
  1.0000  2.0000
  3.0000  4.0000
```

(`Printf.printf "%8.4f"` with the number left out is a function that prints one float, and `Array.iter` applies it to every entry of a row.) Row $1$ is the first output, column $2$ the derivative with respect to the second input: $J_{12} = 2$ (`J.(0).(1)` in the code). A transposed matrix would print $3$ there. `Stage.solve` forms $I - \gamma J$ with the same indices, and `Linalg.solve` treats row $i$ as equation $i$, so orientation matters in all three places.

## 4. The step of `Jac.forward`

[`src/numerics/jac.ml`](../../src/numerics/jac.ml) has a private step function and the code that builds the perturbed point for column $j$:

```text
let step yj = 1e-8 *. (1. +. Float.abs yj)

let column j =
  let yp = Array.mapi (fun k yk -> if k = j then yk +. step yk else yk) y in
  (1. /. (yp.(j) -. y.(j)), f yp)
```

`yp` is a copy of `y` with component `j` replaced by $y_j + \delta_j$ (`step yj` is $\delta_j$; `Array.mapi` builds a new array; the library never writes into an existing one). The column is kept as the pair (reciprocal of the stored perturbation, `f yp`); entry `i` of the Jacobian is then `inv *. (fp.(i) -. fy.(i))`, built row by row.

**The size, $\delta_j = 10^{-8}(1 + \lvert y_j \rvert)$.** $10^{-8}$ is a round number near $\sqrt{\varepsilon} = 1.5 \times 10^{-8}$, the optimum of section 2 for $\lvert f \rvert$ and $\lvert f'' \rvert$ of order 1. The factor $1 + \lvert y_j \rvert$ makes the step relative to the size of $y_j$ when that is large and about $10^{-8}$ near zero. Relative, because floats have relative precision: at $y_j = 10^{10}$ a step of $10^{-8}$ is absorbed (section 1), the perturbation is zero, the divisor is $0$ and the column is NaN; with the factor the step is about $100$, far above the $1.9 \times 10^{-6}$ spacing. The $1$ in $1 + \lvert y_j \rvert$ keeps the step from collapsing to zero as $y_j$ approaches 0.

**Dividing by what was stored.** The sum $y_j + \delta_j$ is rounded to a float, so the perturbation actually applied, $(y_j + \delta_j) - y_j$, differs from $\delta_j$ by up to half an ulp of $y_j$: relative to the step, up to about $10^{-8}$ for large $y_j$, as big as the truncation error. The code therefore divides by `yp.(j) -. y.(j)`, which costs nothing. For $f(y) = y$, whose slope is exactly 1, the numerator is exactly the stored perturbation, so the stored divisor gives 1 up to the rounding of one reciprocal. This prints how far the nominal step is from the stored one:

```ocaml
let () =
  List.iter
    (fun y ->
      let nominal = 1e-8 *. (1. +. Float.abs y) in
      let stored = y +. nominal -. y in
      Printf.printf "y = %-8g relative difference between stored and nominal step: % .1e\n" y
        ((stored -. nominal) /. nominal))
    [ 0.; 1.; 123.456; 1e6; 1e10 ]
```

```text
y = 0        relative difference between stored and nominal step:  0.0e+00
y = 1        relative difference between stored and nominal step:  5.0e-09
y = 123.456  relative difference between stored and nominal step: -2.9e-09
y = 1e+06    relative difference between stored and nominal step:  2.1e-09
y = 1e+10    relative difference between stored and nominal step: -1.0e-10
```

Dividing by the nominal step would put those relative errors into the Jacobian. `Adaptive.integrate` applies the same idea to time: the step it takes is the $(t + \mathtt{dt}) - t$ that the clock really made ([chapter 5](05-step-control.md), section 9).

**Where the compromise is weak.** Near zero the step is absolute, so for a small variable it is a noticeable fraction of it, and the truncation error of a strongly curved function shows. Robertson's $y_2$ never exceeds about $4 \times 10^{-5}$. Its term $3 \times 10^7 y_2^2$ has derivative $6 \times 10^7 y_2 = 600$ at $y_2 = 10^{-5}$ and curvature $f'' = 6 \times 10^7$, so the truncation error $\frac{\delta}{2} f'' = 0.3$ already appears in the fourth digit: `Jac.forward (fun y -> [| 3e7 *. y.(0) *. y.(0) |]) [| 1e-5 |]` returns `600.3000` where the exact derivative is $600$. This costs a little speed, not the answer: a slightly wrong Jacobian does not change the equation Newton solves ([chapter 2](02-newton.md)).

## 5. Worked example: the canary and its rounding floor

The corpus' canary is $y' = -\Lambda y$ with $\Lambda = \mathrm{diag}(1, 100, 10^4)$ ([chapter 6](06-the-corpus.md)), so the true Jacobian is $\mathrm{diag}(-1, -100, -10^4)$. Take $y_0 = (1, 1, 1)$ and the entry `J.(2).(2)` in code, $J_{33}$ in the text (the text numbers components from 1, the code from 0):

1. The nominal step is $\delta_3 = 10^{-8}(1 + 1) = 2 \times 10^{-8}$. Floats just above 1 are $2^{-52}$ apart and $2 \times 10^{-8} / 2^{-52} = 90071992.5\ldots$, so the nearest float to $1 + 2 \times 10^{-8}$ is $1 + 90071993 \times 2^{-52}$ and the perturbation actually applied is $2.0000000100495186 \times 10^{-8}$, five parts in a billion above the nominal one.
2. $f_3 = -10^4 y_3$ is linear, so the quotient has no truncation error at all. The exact product is $-10^4 \times (1 + 2.0000000100495186 \times 10^{-8})$, about $-10000.000200000001$. Floats near $10^4$ are $1.8 \times 10^{-12}$ apart, so the product is rounded, with an error of $6.0 \times 10^{-13}$.
3. The numerator $f_3(y + \delta_3 e_3) - f_3(y)$ is $-0.00020000000040454324$ instead of the true $-0.00020000000100495186$: only $2 \times 10^{-4}$ out of values of size $10^4$, so the $6 \times 10^{-13}$ error is a visible fraction of it.
4. Dividing by the stored step gives $-9999.99997$ instead of $-10000$: an error of $6.0 \times 10^{-13} / (2 \times 10^{-8}) = 3.0 \times 10^{-5}$.

The error is pure round-off from one rounded product: the price of a forward difference with this step, not a bug. Its bound is half an ulp of $10^4$ over the step, $9.1 \times 10^{-13} / (2 \times 10^{-8}) = 4.5 \times 10^{-5}$; the rule of thumb $\varepsilon \lvert f \rvert / \delta$ gives $1.1 \times 10^{-4}$. Relative to the entry it is $3.0 \times 10^{-9}$. Check it (needs the `problems.ml` link):

```ocaml
open Numerics

let () =
  let open Problems.Canary in
  let j = Jac.forward (rhs 0.) y0 in
  Printf.printf "J_33 = %.10f\n" j.(2).(2);
  Array.iteri (fun i l -> Printf.printf "error of J_%d%d = % .2e\n" (i + 1) (i + 1) (j.(i).(i) -. -.l)) lambda
```

```text
J_33 = -9999.9999699796
error of J_11 =  1.11e-16
error of J_22 =  4.44e-08
error of J_33 =  3.00e-05
```

The smaller entries have smaller absolute errors, because the rounding error of a product is proportional to the product. $J_{11}$ is one ulp from $-1$: the code multiplies by the rounded reciprocal `1 / step` instead of dividing, which costs at most one more rounding. At the origin, in contrast, all three diagonal entries come out exact.

## 6. Why the corpus checks the origin absolutely and $y_0$ relatively

[`test/corpus.ml`](../../test/corpus.ml) compares `Jac.forward` with the analytic Jacobian at two points, with a different yardstick at each:

```
jac canary at origin, max entry error < 1e-6: true
jac canary at y0, max entry error relative to |J_ij| < 1e-6: true
```

- **At the origin** $f(0) = 0$ exactly and $f(0 + \delta_j e_j)$ is a single rounded product, so there is no cancellation: the subtraction $f(0 + \delta_j e_j) - 0$ loses nothing. The entries come out exact, and the absolute bound $10^{-6}$ is a strict test of layout, step and divisor: a relative mistake above $10^{-10}$ in the $10^4$ entry would break it.
- **At $y_0 = (1, 1, 1)$** $f(y_0) = (-1, -100, -10^4)$ is large and $f(y_0 + \delta_j e_j)$ is close to it, so cancellation applies and the floor of section 5 appears: an absolute error of $3.0 \times 10^{-5}$ on the $10^4$ entry. An absolute bound of $10^{-6}$ would fail a perfectly correct Jacobian. The test therefore divides each entry's error by $1 + \lvert J_{ij} \rvert$ (the label says `|J_ij|`; the $1$ also turns the zero off-diagonal entries into plain absolute errors), which gives $3.0 \times 10^{-9}$ for that entry, below $10^{-6}$ by a factor of about 300.

Two instruments, then: the first is a strict check of the structure, the second a check of realistic accuracy. When a bound seems too tight, compute where the floor of that case lies: it tells whether the code or the yardstick is wrong, and a bound is never loosened to make a line pass ([AGENTS.md](../../AGENTS.md), H2).

Both lines also pass with a step of $10^{-6}$, so a third line pins the step itself, line 41:

```
jac forward difference of y^3 at y = -2 (exact 12): 11.9999998224
```

A forward difference of a cubic is $12 - 6\delta + \delta^2$ plus round-off; with $\delta = 10^{-8}(1 + 2) = 3 \times 10^{-8}$ that is $11.99999982$, and the line prints `11.9999998224` with ten decimals. A step of $10^{-6}(1 + \lvert y \rvert)$ prints `11.9999820002`, a step without the `abs` (it is $-10^{-8}$ at this $y$) `12.0000000000`, and a backward difference `12.0000001776`.

## 7. NaN and infinity

**Where they come from.** Float arithmetic in OCaml never raises (integer division by zero does): `1. /. 0.` is `infinity`, `0. /. 0.` is `nan` ("not a number"), `exp 710.` overflows to `infinity`, `infinity -. infinity` and `0. *. infinity` are `nan`, `sqrt (-1.)` and `log (-1.)` are `nan`. In the integrator they come from a right-hand side that divides by a component that reaches zero, from a Newton trial point outside the function's domain (chapter 2's $\sqrt{x} - 1$ at $x = -3$), and from overflow when a trial state is huge (van der Pol squares $y_1$) or when a Newton step is added to an $x$ near the largest float ([chapter 2](02-newton.md), section 6). The last row of section 2's table is another source, a step too small for its point, which the relative step avoids. `Float.is_finite x` is false for `nan` and both infinities, and `Vec.finite v` applies it to every component.

**Comparisons with NaN are false.** `nan = nan`, `nan < 1.`, `nan <= 1.` and `nan >= 1.` are all `false`, and `<>` is simply the opposite of `=`: `nan <> x` is `true` for every `x`. So an acceptance test must be written as "accept if the comparison is true", so that NaN falls into the reject branch. The library does: `Halving.acceptable` keeps a step only if the weighted error is $\le \mathrm{tol}$, so a NaN estimate is rejected; Newton's damping accepts a trial point only if `norm_inf fx' <= ...`. The mirror-image test, "reject if $\mathrm{err}$ exceeds $\mathrm{tol}$, otherwise accept", accepts NaN. It is a classic bug.

**`max` and equality.** `Float.max` and `Float.min` return `nan` if either argument is `nan`, so one NaN component makes `Vec.norm_inf` (and the weighted error of `Halving.acceptable`) `nan`, which then fails the tests above. For the same reason the explicit `Vec.finite fx'` in Newton's line search is redundant: `norm_inf fx'` is already `nan` or `infinity` whenever an entry is not finite, and `nan <= bound` is false. The generic `max` (no `Float.` prefix) does not do this and is not even symmetric: `max nan 1.` is `1.` but `max 1. nan` is `nan`. A `norm_inf` built on it would lose a NaN that is followed by a larger entry, which is why the example below puts the NaN first. Equality: `=` says `nan = nan` is false, even for the same value, and `[| nan |] = [| nan |]` is false; `Float.equal nan nan` and `compare nan nan = 0` treat NaN as equal to itself; and `0. = -0.` is true. The soak test compares results with `=` ([docs/testing.md](../testing.md)), so a result containing NaN prints `identical: false`.

```ocaml
open Numerics

let () =
  let nan = Float.nan in
  Printf.printf "nan = nan: %b   nan <= 1: %b   nan <> nan: %b\n" (nan = nan) (nan <= 1.) (nan <> nan);
  Printf.printf "Float.max nan 1 = %f   Float.max 1 nan = %f\n" (Float.max nan 1.) (Float.max 1. nan);
  Printf.printf "max nan 1 = %f   max 1 nan = %f\n" (max nan 1.) (max 1. nan);
  Printf.printf "Vec.norm_inf [|nan; 1|] = %f\n" (Vec.norm_inf [| nan; 1. |]);
  Printf.printf "[|nan|] = [|nan|]: %b   Float.equal nan nan: %b\n" ([| nan |] = [| nan |]) (Float.equal nan nan);
  Printf.printf "1/0 = %f   0/0 = %f   inf - inf = %f   exp 710 = %f\n" (1. /. 0.) (0. /. 0.) (infinity -. infinity) (exp 710.)
```

```text
nan = nan: false   nan <= 1: false   nan <> nan: true
Float.max nan 1 = nan   Float.max 1 nan = nan
max nan 1 = 1.000000   max 1 nan = nan
Vec.norm_inf [|nan; 1|] = nan
[|nan|] = [|nan|]: false   Float.equal nan nan: true
1/0 = inf   0/0 = nan   inf - inf = nan   exp 710 = inf
```

**What the library does with them.** `Newton.solve` returns `Error Nan` when $x$ or $G(x)$ is not finite at the start of an iteration; after the first iteration that does not happen in practice, because damping refuses any trial point whose residual is not finite and shortens the step instead (if every shortened step is refused, the result is `Error Diverged`, and so it is for a non-finite Newton step; corpus lines 39 and 40 pin the names at the start and for a step that overflows, and line 47 pins the refusal: the full step of $\ln x$ from $x = 3$ lands where $\ln$ is `nan`). It returns `Error Nan` too when a converged step `x + dx` overflows to infinity, which the stopping test, being relative, cannot see. `Stepper.fixed` and `Adaptive.integrate` return `Error Nan` if $y_0$ or $f(t_0, y_0)$ is not finite (a fixed-step run with an empty span included; corpus lines 36 and 44); `Adaptive.integrate` treats a failed Newton solve like a bad error estimate: the step is rejected and halved ([chapter 5](05-step-control.md)). All numerical failures are `Fail.t` values, never exceptions.

## 8. Fused multiply-add and the last digits

Many CPUs can compute $a \times b + c$ as one operation with one rounding, the **fused multiply-add** (FMA). On arm64 (the 64-bit ARM processors of Apple-silicon Macs, phones and many servers), `ocamlopt` emits it when a product is a direct operand of an addition or subtraction (`a +. b *. c`, `a -. b *. c`): the product is not rounded before the sum. Binding the product with `let` first rounds twice and can give a different result. This has been checked on arm64 with OCaml 5.5.0 only; other targets may behave differently. On arm64 this program prints a different result for the same arithmetic:

```ocaml
let fused a b c = a -. (b *. c)

let unfused a b c =
  let p = b *. c in
  a -. p

let () =
  let b = 1. +. Float.epsilon and c = 1. -. Float.epsilon in
  Printf.printf "inline     %h\nlet-bound  %h\n" (fused 1. b c) (unfused 1. b c)
```

```text
inline     0x1p-104
let-bound  0x0p+0
```

The exact value of $1 - (1 + \varepsilon)(1 - \varepsilon)$ is $\varepsilon^2 = 2^{-104}$: the fused form keeps it, while in the other the product rounds to $1$ first (`%h` prints a float in hexadecimal). The IEEE rules fix the result of each single operation, but not whether the compiler fuses two of them, and library functions such as `exp` may round differently in different C libraries.

What this means for changing the code:

- Whether a product fuses depends on the shape of the expression, not on how it reads. In the current build on arm64 the products in `Vec.axpy` (`a *. x.(i) +. y.(i)`), in the stage Jacobian ($I - \gamma J$), in the perturbed point $y_j + \delta_j$ of `Jac.forward`, in the grid time $t_0 + k h$ of `Stepper.fixed` and twice in the van der Pol right-hand side of `test/problems.ml` fuse, as do $1 - c \lambda$ in `Newton` ($c$ is `armijo` in the code) and one product in `Bdf2.coeffs` (both are a float times a power of two, hence exact, so fusing them changes nothing); the elimination step of `Linalg` does not. To see for yourself, compile a standalone file with `ocamlopt -S -c file.ml` and look in `file.s` for the mnemonics `fmadd`, `fmsub`, `fnmadd` and `fnmsub`; for a library module, add `(ocamlopt_flags (:standard -S))` to the `library` stanza of `src/dune` (a solver module) or of `src/numerics/dune` (a kernel module: `Vec`, `Jac`, `Newton`, `Linalg`) in a scratch copy and look in `_build/default/src/.vstiff.objs/native/`, respectively `_build/default/src/numerics/.numerics.objs/native/`, after `dune build`.
- Moving a product into or out of a sum (binding it with `let`, factoring, reordering operands) can change the last bits of a result. That is allowed: checks are about correctness, not bits, and no check may depend on those bits ([AGENTS.md](../../AGENTS.md), H1). If an expected line moves after such an edit, find out whether behaviour changed or the line prints more than its check needs. Run `dune runtest` after any arithmetic edit.
- Today the recorded lines do not hinge on fusion: binding every fused product of both libraries and of `test/problems.ml` with `let` leaves the output of both test programs unchanged on arm64, and both expected files pass on Linux x86_64, where nothing fuses, and on Linux arm64 (OCaml 5.5.1, checked 2026-10-03). That is a fact about this code, not a promise about the next change, and no check may rely on it. A rounding difference is not a bug; a line that moves because of one prints more than its check needs (H1). What stays fixed is determinism within one build: the same binary, given the same inputs, returns the same outputs, and the soak test checks it by comparing ten runs with `=`.

## In the code

| Idea | Where |
|------|-------|
| The step $10^{-8}(1 + \lvert y_j \rvert)$ | the private function `step` in [`src/numerics/jac.ml`](../../src/numerics/jac.ml) |
| The forward-difference Jacobian: perturb one component, divide by the stored perturbation, build the rows | `Jac.forward`, local function `column` |
| The stage Jacobian $I - \gamma J$ | `Stage.solve` in [`src/stage.ml`](../../src/stage.ml): local `jacobian` calls `Jac.forward f x` and forms `(if i = j then 1. else 0.) -. gamma *. v` |
| Finite checks | `Vec.finite`, `Vec.norm_inf` in [`src/numerics/vec.ml`](../../src/numerics/vec.ml); `Newton.solve` ([`src/numerics/newton.ml`](../../src/numerics/newton.ml)) at the start of an iteration, on the step, on each trial residual and on a converged `x + dx` |
| The accept test that rejects NaN | `Halving.acceptable` in [`src/halving.ml`](../../src/halving.ml) |
| `Error Nan` at the start of a run | `Adaptive.integrate` and `Stepper.fixed`: `Vec.finite p.y0 && Vec.finite f0`, with `f0 = p.rhs p.t0 p.y0`, in [`src/adaptive.ml`](../../src/adaptive.ml) and [`src/stepper.ml`](../../src/stepper.ml); `Fail.Nan` in [`src/numerics/fail.ml`](../../src/numerics/fail.ml) |
| The Jacobian checks | the `jacobian`, `orientation` and `jacobian_step` cases of [`test/corpus.ml`](../../test/corpus.ml); the `jac` lines of [`test/corpus.expected`](../../test/corpus.expected) |

## Check yourself

1. **Why is $10^{10} + 10^{-8} = 10^{10}$ in floating point, and what would `Jac.forward` do at $y_j = 10^{10}$ if its step were the absolute $10^{-8}$?**
   The floats near $10^{10}$ are $1.9 \times 10^{-6}$ apart, so a change of $10^{-8}$ is rounded away. The perturbed point would equal `y`, the divisor `yp.(j) -. y.(j)` would be $0$, and the column would be infinity times $0$: NaN. The relative factor makes the step about $100$.

2. **In the sweep table, at which $\delta$ is the error smallest, and which error dominates at $\delta = 10^{-12}$?**
   Around $10^{-8}$, where it is about $2 \times 10^{-9}$. At $10^{-12}$ the round-off error $\varepsilon \lvert f \rvert / \delta = 3 \times 10^{-4}$ dominates (the measured error is $9 \times 10^{-5}$); the truncation error is $4 \times 10^{-14}$.

3. **The canary is linear, so why does `Jac.forward` not return its Jacobian exactly at $y_0 = (1, 1, 1)$?**
   There is no truncation error, but the product $-10^4 \times (1 + \delta)$ is rounded, with an error up to half an ulp of $10^4$ ($9.1 \times 10^{-13}$). Dividing by $\delta = 2 \times 10^{-8}$ turns that into up to $4.5 \times 10^{-5}$ (measured: $3.0 \times 10^{-5}$).

4. **Why does the corpus use the absolute bound $10^{-6}$ at the origin but the bound relative to $1 + \lvert J_{ij} \rvert$ at $y_0$?**
   At the origin there is no cancellation and the entries are exact, so a strict absolute bound checks the structure. At $y_0$ cancellation limits the $10^4$ entry to about $3 \times 10^{-5}$ absolute, so an absolute $10^{-6}$ would fail a correct Jacobian; relative to $1 + \lvert J_{ij} \rvert$ the error is $3 \times 10^{-9}$.

5. **What is wrong with `if not (est > tol) then accept else reject`, and what does the code write instead?**
   For a NaN `est`, `est > tol` is false, so `not (...)` is true and the step is accepted. The code writes `... <= tol`: a NaN estimate makes it false, so the step is rejected.

6. **You replace `a *. x.(i) +. y.(i)` by a version that binds the product first. Which test output could change, and why?**
   Only the last bits of results: on a machine that fuses, the original rounds once and the new form twice. No check may depend on those bits ([AGENTS.md](../../AGENTS.md), H1), and today no expected line moves (section 8); a line that did move would be printing more than its check needs.

Next: [4. BDF methods](04-bdf.md), where the stage equation $x = \psi + \gamma\thinspace f(t_{n+1}, x)$ that Newton solves comes from.
