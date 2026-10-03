# Exercises

Graded exercises, from running the tests to writing a second controller, then starter contributions. **A** warm-ups, **B** reading the code, **C** probes (small programs that measure), **D** breaking the code in a scratch copy to see which tests notice, **E** a capstone against the contracts, **F** starter contributions. [onboarding.md](onboarding.md) says which to do on which day. Each exercise has a goal, hints, where to look and a done-when line. None gives BDF coefficient values: they live in `Bdf2.coeffs`, which you can evaluate, or derive as [numerics/04-bdf.md](numerics/04-bdf.md) and day 6 of [onboarding.md](onboarding.md) ask. Exercises that add or change a corpus case (section F) are discussed first: new cases are added deliberately, one at a time, each with a stated purpose ([CONTRIBUTING.md](../CONTRIBUTING.md)).

## Setup

**The probe project.** A probe is a throwaway program that calls the library. It lives in a dune project outside the repository that links to `src/`, so nothing is written inside the repository. From the repository root, with your opam switch active:

```sh
REPO=$(pwd)
mkdir -p ../vstiff-scratch/p && cd ../vstiff-scratch
printf '(lang dune 3.0)\n' > dune-project
ln -sfn "$REPO/src" src
ln -sfn "$REPO/test/problems.ml" p/problems.ml
printf '(executable (name probe) (libraries vstiff))\n' > p/dune
```

The recipe leaves the shell in the probe project (`cd "$REPO"` returns to the repository root). Write a program in `p/probe.ml`, starting with `open Vstiff`, and run it with `dune build --root . ./p/probe.exe && ./_build/default/p/probe.exe` (`--root .` stops dune from picking an enclosing project as its root). The link to `problems.ml` provides `Problems.Canary`, `Problems.Logistic`, `Problems.VanDerPol` and `Problems.Robertson`, each with `rhs`, `y0` and a ready `problem` record; Canary and Logistic also have `exact`. Probes may use `ref` and print, which the library avoids. An unused top-level definition is an error in dune's dev profile (warning 32): use what you define. Run anything that might not finish under a time limit, applied to the built program and not to `dune runtest` (killing dune can leave the test programs running): `perl -e 'alarm 120; exec @ARGV' ./_build/default/p/probe.exe`. A killed program prints nothing it has not flushed: end a `Printf.printf` format with `%!` for the lines that must survive.

**The scratch copy.** A copy of the committed files that you can break and delete. From the repository root: `mkdir ../vstiff-mutation && git archive HEAD | tar -x -C ../vstiff-mutation`; then, inside the copy, `dune build --root . && dune runtest --root .` (silence means green). Never run `dune promote` in it to make a test pass. To run a probe against the modified library, point the probe project's `src` link at the copy (`ln -sfn ../vstiff-mutation/src src` in `../vstiff-scratch`) and restore it afterwards.

**Asking the compiler for a type.** Annotate a name with a type it cannot have, as in `let () = ignore (Vec.axpy : int)`. The error says `The value Vec.axpy has type float -> Vec.t -> Vec.t -> Vec.t but an expression was expected of type int`; other functions may print the library prefix, as in `Vstiff.Vec.t`. An editor with OCaml support shows the same on hover.

## A. Warm-ups

### A1. Run the tests
**Goal.** See a green run and know what its silence means. **Hints.** From the repository root run `dune build`, `dune runtest`, `dune exec ./test/corpus.exe` and `dune exec ./test/soak.exe`. `dune build` runs the two programs, when their output is out of date, to record it; `dune runtest` compares it with [corpus.expected](../test/corpus.expected) and [soak.expected](../test/soak.expected). `dune exec ./test/corpus.exe | diff test/corpus.expected -` compares by hand (`<` lines are expected, `>` lines were printed). **Look in.** [README.md](../README.md), [testing.md](testing.md). **Done when.** You can say what a silent `dune runtest` means, which program each expected file belongs to and which of the commands compares; what every line checks is the table of day 8 ([numerics/06-the-corpus.md](numerics/06-the-corpus.md)).

### A2. Break a test on purpose
**Goal.** Know what a failure looks like before you meet a real one. **Hints.** In a scratch copy change the last digit of `1.414213562373` in `test/corpus.expected` to 4, predict the diff, run `dune runtest --root .` and check `echo $?` at once. `-` lines come from the expected file, `+` lines from the program. Then restore the file (from the repository root, `git show HEAD:test/corpus.expected > ../vstiff-mutation/test/corpus.expected`). **Look in.** [testing.md](testing.md). **Done when.** You saw `-newton quadratic x0=1: Ok [1.414213562374]`, `+newton quadratic x0=1: Ok [1.414213562373]` and exit status 1, and can say why `dune promote` would be wrong here.

### A3. Newton by hand
**Goal.** Run Newton's method on paper until it feels mechanical. **Hints.** Solve `x^2 - 2 = 0` from `x = 1` with `x - G(x) / G'(x)`: three iterates. Then run this probe, which prints every point where the residual is evaluated:

```ocaml
open Vstiff

let () =
  let f x =
    Printf.printf "f(%.17g)\n" x.(0);
    [| (x.(0) *. x.(0)) -. 2. |]
  in
  let jac x = [| [| 2. *. x.(0) |] |] in
  match Newton.solve f jac [| 1. |] with
  | Ok x -> Printf.printf "Ok %.17g\n" x.(0)
  | Error e -> Printf.printf "Error %s\n" (Fail.to_string e)
```

The trial points are 1, 1.5, 1.4166666666666667, 1.4142156862745099, 1.4142135623746899, and the root is returned without another evaluation because the last step was already below the tolerance. Start from `0.` and the 1x1 Jacobian is singular. Then switch to `atan` (`f x = [| atan x.(0) |]`, `jac x = [| [| 1. /. (1. +. x.(0) *. x.(0)) |] |]`, start `3.`): the full step lands near -9.5, where the residual is larger, the half step (near -3.2) fails too, and the quarter step (near -0.12) is the first accepted. **Look in.** [src/newton.ml](../src/newton.ml), [numerics/02-newton.md](numerics/02-newton.md). **Done when.** Your iterates match, and you can name which exit of `Newton.solve` ends each of the three Newton corpus lines: the zero-residual test, the small-step test and the singular Jacobian.

### A4. Gaussian elimination by hand
**Goal.** See how `Linalg.solve` exchanges rows without writing to any array. **Hints.** Trace `Linalg.solve [| [| 1.; 2. |]; [| 4.; 5. |] |] [| 5.; 14. |]`: which row is the pivot, what one-row system is left, how the first unknown follows from the second. `swapped p i` names which original row sits at position `i` after rows 0 and `p` trade places. Check with a probe (`Some [|1.; 2.|]`), then try `[| [| 1e-20; 1. |]; [| 1.; 1. |] |]` with `[| 1.; 2. |]` and `[| [| 0.; 1. |]; [| 1.; 0. |] |]` with `[| 3.; 4. |]`: the answers are (1, 1) and (4, 3). **Look in.** [src/linalg.ml](../src/linalg.ml), [numerics/02-newton.md](numerics/02-newton.md). **Done when.** Your trace matches, and you can say why the row starting with `1e-20` must not be the pivot.

### A5. Floating-point surprises
**Goal.** See round-off yourself. **Hints.** Predict, then print in a probe: `0.1 +. 0.2` with `%.17g`; `0.1 +. 0.2 = 0.3`; `1. +. Float.epsilon /. 2. = 1.`; `Float.nan = Float.nan`; `(1e16 +. 1.) -. 1e16`. Neighbouring floats near 1 are `Float.epsilon` apart, near `1e16` they are 2 apart. **Look in.** [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md). **Done when.** You predicted `0.30000000000000004`, `false`, `true`, `false` and `0`, and can explain each.

## B. Reading the code

### B1. Trace one adaptive attempt
**Goal.** Know who calls whom and what is passed down. **Hints.** Start at `go` in `Adaptive.integrate` and follow one attempt, once on the first step (history `Start`) and once with history, down to `Linalg.solve`. Write each function with its file and, for every stage solve, `psi`, `gamma` and the starting guess. Then follow what the controller does with the outcome. **Look in.** [src/adaptive.ml](../src/adaptive.ml), [src/bdf2.ml](../src/bdf2.ml), [src/bdf1.ml](../src/bdf1.ml), [src/stage.ml](../src/stage.ml), [src/newton.ml](../src/newton.ml), [src/halving.ml](../src/halving.ml). **Done when.** Your chain shows one stage solve plus one extra `rhs` call on the first attempt (the explicit Euler point of the error estimate) and two stage solves afterwards, and each Newton iteration calling `Jac.forward` and then `Linalg.solve`.

### B2. Where does each failure come from?
**Goal.** Map every `Fail.t` constructor to the code that creates it. **Hints.** `grep -rnE 'Diverged|Nan|StepRejected' src`. Which can `Stepper.fixed` return, and which can `Adaptive.integrate` return with `Halving`? What happens to a failed Newton solve inside an attempt? **Look in.** [src/newton.ml](../src/newton.ml), [src/adaptive.ml](../src/adaptive.ml), [src/halving.ml](../src/halving.ml), [src/ode.mli](../src/ode.mli). **Done when.** You wrote the table before reading this answer and C4 reproduces each constructor. Answer: `Newton.solve` creates `Diverged` (four exits) and `Nan`, and `Stepper.fixed` passes on what its method returns; `Adaptive.integrate` creates `Nan` for a non-finite start; `Halving` creates `StepRejected`. A failed solve reaches the controller as `Solver e` and counts as a rejection, so with `Halving` the adaptive driver returns only `Nan` or `StepRejected`.

### B3. Count the steps of `Stepper.fixed`
**Goal.** Predict how many steps a fixed-step run takes. **Hints.** For `t0 = 0` and `t_end = 1`, how many steps, and how long, for `dt = 0.45`, `0.12` and `2`? For an empty span? For `dt = 0` (it raises, so the probe leaves it out; try it last)? This method does no numerics: it adds `h` to `y.(0)` and 1 to `y.(1)`, so the final state reports the total time and the count.

```ocaml
open Vstiff

module Ticker : Ode.Method = struct
  type history = unit

  let start = ()
  let step _rhs h () (at : Ode.point) = Ok ([| at.y.(0) +. h; at.y.(1) +. 1. |], ())
end

let run (t0, t_end, dt) =
  match Stepper.fixed (module Ticker) ~dt { Ode.rhs = (fun _ y -> y); t0; t_end; y0 = [| 0.; 0. |] } with
  | Ok y -> Printf.printf "[%g, %g] dt = %g: %g steps, total %g\n" t0 t_end dt y.(1) y.(0)
  | Error e -> print_endline (Fail.to_string e)

let () = List.iter run [ (0., 1., 0.45); (0., 1., 0.12); (0., 1., 2.); (1., 1., 0.1) ]
```

**Look in.** [src/stepper.ml](../src/stepper.ml), [src/check.ml](../src/check.ml). **Done when.** You predicted `n = round (span / dt)` (at least one step for a non-empty span) with `h = span / n`: two steps of 0.5, eight of 0.125, one of 1, none for the empty span, and `Invalid_argument` for `dt = 0`.

### B4. Run `Halving` by hand
**Goal.** Say what the controller proposes next after any run of accepts and rejects. **Hints.** From `dt0 = 0.08` with a large `dt_max`, predict the step tried and the next proposal for the outcomes A A A A R R A A A (A accepted, R rejected). The probe drives the controller directly through the `Ode.Controller` interface.

```ocaml
open Vstiff

let step c outcome =
  let h = Halving.proposal c in
  let c = if outcome = 'A' then Halving.accepted c else Result.get_ok (Halving.rejected c Ode.Too_large ~at:0. ~h) in
  Printf.printf "%c: tried %g, next proposal %g\n" outcome h (Halving.proposal c);
  c

let () =
  let c = Halving.init ~tol:1e-6 ~dt0:0.08 ~dt_max:1e6 ~max_rejects:50 in
  ignore (List.fold_left step c [ 'A'; 'A'; 'A'; 'A'; 'R'; 'R'; 'A'; 'A'; 'A' ])
```

**Look in.** [src/halving.ml](../src/halving.ml), [numerics/05-step-control.md](numerics/05-step-control.md). **Done when.** The steps tried are 0.08, 0.08, 0.08, 0.16, 0.16, 0.08, 0.04, 0.04, 0.04 and the proposals afterwards 0.08, 0.08, 0.16, 0.16, 0.08, 0.04, 0.04, 0.04, 0.08. Then list the ratios ω = h / h_prev BDF2 would see along the way (none exceeds 2) and say why a rejection halves the step that failed rather than the old proposal.

## C. Probes

### C1. Measure the order
**Goal.** Turn "halving `dt` divides the error by `2^p`" into a measurement of `p`. **Hints.** If `error ~ C dt^p` then `p ~ log2 (e(dt) / e(dt/2))`.

```ocaml
open Vstiff

let error (module M : Ode.Method) dt =
  let open Problems.Logistic in
  match Stepper.fixed (module M) ~dt problem with
  | Ok y -> Vec.norm_inf (Vec.sub y (exact 5.))
  | Error e -> failwith (Fail.to_string e)

let order name m dt =
  let e = error m dt and e2 = error m (dt /. 2.) in
  Printf.printf "%s dt = %-5g error %.3e, order %.2f\n" name dt e (Float.log (e /. e2) /. Float.log 2.)

let () =
  List.iter (order "bdf1" (module Bdf1)) [ 0.04; 0.02; 0.01 ];
  List.iter (order "bdf2" (module Bdf2)) [ 0.04; 0.02; 0.01 ]
```

**Look in.** The `bdf2 logistic` line of [corpus.expected](../test/corpus.expected), [numerics/04-bdf.md](numerics/04-bdf.md). **Done when.** The orders are close to 1 and 2, and your BDF2 error at `dt = 0.01` matches the corpus line (`3.116e-06`).

### C2. The Jacobian: rounding floor and orientation
**Goal.** Reproduce the rounding error of a forward difference and confirm which way `J` points. **Hints.** The canary is linear, so only round-off remains. At `y0 = (1, 1, 1)` the perturbation is `2e-8`; half an ulp of `1e4` divided by that bounds the error of the `-1e4` entry near 4.5e-5. The second function has off-diagonal entries, so a transpose shows.

```ocaml
open Vstiff

let () =
  let open Problems.Canary in
  let j = Jac.forward (rhs 0.) [| 1.; 1.; 1. |] and z = Jac.forward (rhs 0.) [| 0.; 0.; 0. |] in
  Printf.printf "entry (2,2) minus -1e4 at y0:   %.3e\n" (j.(2).(2) +. 1e4);
  Printf.printf "same entry at the origin:       %.3e\n" (z.(2).(2) +. 1e4);
  let f y = [| y.(0) *. y.(1); y.(0) +. (3. *. y.(1)) |] in
  Array.iter (fun row -> Array.iter (Printf.printf "%8.4f") row; print_newline ()) (Jac.forward f [| 2.; 5. |])
```

**Look in.** [src/jac.ml](../src/jac.ml), the `jac` cases in [test/corpus.ml](../test/corpus.ml), [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md). **Done when.** The first error is about 3.0e-5 (relative 3e-9), the origin entry is essentially exact, the matrix is `[[5, 2], [1, 3]]`, and you can say why the corpus bounds are absolute at the origin and relative at `y0`.

### C3. Count the right-hand-side calls
**Goal.** Measure what an implicit step costs. **Hints.** `Instrument.count rhs` returns a wrapped `rhs` and a function that reads its call count. Wrap `Problems.Canary.rhs`, take one step with `Bdf1.step counted 1e-3 Bdf1.start { Ode.t = 0.; y = Problems.Canary.y0 }` and print the count (`Bdf1.history` is abstract: match the history with `_`). Predict first: one residual evaluation to start, `n + 1` calls per Jacobian, one residual per accepted trial point. If every full step is accepted and the solve ends on the small-step test, `k` iterations cost `k (n + 2)` calls; invert that to find `k`. Then wrap `rhs` for a whole `Adaptive.integrate` run: `{ problem with Ode.rhs = counted }` is the problem record with the counted right-hand side. **Look in.** [src/instrument.ml](../src/instrument.ml), [src/newton.ml](../src/newton.ml), [src/stage.ml](../src/stage.ml). **Done when.** Your count is `k (n + 2)` for a whole `k`, and you can say what one adaptive attempt costs once there is history.

### C4. Provoke each failure
**Goal.** Produce `Nan`, `StepRejected` and `Diverged` on purpose. **Hints.** The lines should read `Error Nan`, `Error StepRejected 1` and `Error Diverged`. The first has a non-finite start. In the second, `max_rejects:0` makes the first rejection already too many, and the first attempt's estimate is far above `tol:1e-14`. The third is one backward Euler step of 2 on `y' = y^2` from 1, whose stage equation `x = 1 + 2 x^2` has no real root. Then raise `max_rejects` to 1000 on `y' = y^2` (`t_end = 2`, `y0 = 1`) and on a right-hand side that is NaN past `t = 0.5`: the count stops mattering, because the floor `16 eps |t|` ends the run.

```ocaml
open Vstiff

let show = function
  | Ok (s : _ Adaptive.solution) -> Printf.sprintf "Ok at t = %g" s.t
  | Error e -> "Error " ^ Fail.to_string e

let adaptive = Adaptive.integrate (module Bdf2) (module Halving)

let () =
  let open Problems.Logistic in
  print_endline (show (adaptive ~tol:1e-6 { problem with y0 = [| nan |] }));
  print_endline (show (adaptive ~max_rejects:0 ~tol:1e-14 problem));
  let square = { Ode.rhs = (fun _ y -> [| y.(0) *. y.(0) |]); t0 = 0.; t_end = 2.; y0 = [| 1. |] } in
  print_endline (match Stepper.fixed (module Bdf1) ~dt:2. square with Ok _ -> "Ok" | Error e -> "Error " ^ Fail.to_string e)
```

**Look in.** [src/fail.mli](../src/fail.mli), [src/halving.ml](../src/halving.ml), B2. **Done when.** You reproduced all three, and can explain why the `adaptive blow-up` and `adaptive rhs NaN` corpus lines say `StepRejected 1` and `StepRejected 46`, not 51.

### C5. Steps and rejections on van der Pol
**Goal.** Measure how the work grows as the tolerance tightens. **Hints.** Run `Adaptive.integrate (module Bdf2) (module Halving) ~tol Problems.VanDerPol.problem` for `tol` 1e-3, 1e-4 and 1e-5 and print `s.stats.accepted_steps` and `s.stats.rejected_steps`; they are not in the expected files, so they are yours to measure. Compute the share of rejected attempts, then add `~dt_max:1.` and see what changes. **Look in.** [numerics/05-step-control.md](numerics/05-step-control.md), the `vdp` line of [corpus.expected](../test/corpus.expected). **Done when.** You have a table of your own numbers and can say why the corpus line asserts only `rejected >= 1`.

## D. Breaking the code in a scratch copy

Make a scratch copy and check it is green. For each change predict which lines of the expected files differ before you run `dune runtest --root .`, then restore the file (from the repository root, `git show HEAD:src/jac.ml > ../vstiff-mutation/src/jac.ml`, say). Never promote. **Look in** [testing.md](testing.md) and [numerics/06-the-corpus.md](numerics/06-the-corpus.md) for what each line pins.

### D1. Swap two Jacobian rows
**Change.** In `src/jac.ml` replace the last expression of `forward`:

```diff
-  Array.mapi (fun i fyi -> Array.map (fun (inv, fp) -> inv *. (fp.(i) -. fyi)) cols) fy
+  let swap i = if Array.length fy = 3 && i <> 1 then 2 - i else i in
+  Array.mapi (fun i _ -> Array.map (fun (inv, fp) -> inv *. (fp.(swap i) -. fy.(swap i))) cols) fy
```

**Hints.** Five lines change: the two `jac canary` lines and the `jac non-symmetric` line turn `false`, `bdf2 canary t=1 dt=1e-3` becomes `Error Diverged`, and the `bdf1 canary` line stays `true` while its error moves from `3.68e-07` to `3.83e-07`: Newton's root does not depend on the Jacobian, only its speed does. **Done when.** You can say what kind of defect each of the five lines is sensitive to.

### D2. Break pivoting
**Change.** In `src/linalg.ml` replace `pivot` by `let pivot (_ : matrix) = 0`. **Hints.** Only the two `linalg` lines change, to `None` and to a wrong first component; no other line depends on the pivot choice. **Done when.** You can say why a pivot of `1e-20` loses the answer.

### D3. Change the doubling rule
**Change.** In `accepted` in `src/halving.ml` turn both `3`s (`streak = 3` and `streak mod 3`) into 2. **Hints.** Three lines change: `adaptive rhs NaN` and both `adaptive logistic` step-count lines; they are the pins on the step-control rules. **Done when.** You can say what the doubling rule guarantees about ω (never above 2) and why that matters ([numerics/04-bdf.md](numerics/04-bdf.md), zero-stability).

### D4. Keep backward Euler's answer
**Change.** In `Bdf2.step_with_error` return `be` instead of `y` as the first component of the BDF2 branch. **Hints.** The `robertson t=1e4 tol=1e-6 accuracy` line turns `false` with an error of `1.6e-04` against its bound of `1e-5`, and the logistic `dt0=0.5` line changes; the other Robertson line, with bound `1e-3`, still passes. **Done when.** You can say why the tighter line exists.

### D5. Transpose the Jacobian
**Change.** Make `forward` return columns instead of rows: `Array.map (fun (inv, fp) -> Array.mapi (fun i fyi -> inv *. (fp.(i) -. fyi)) fy) cols`. **Hints.** Build with `dune build --root . ./test/corpus.exe`, then limit the executable: `perl -e 'alarm 60; exec @ARGV' ./_build/default/test/corpus.exe`. Newton still reaches the same root but far more slowly, so the Robertson case stops finishing, a full corpus run stalls before the orientation line, and the limit kills it before it prints anything. Run the van der Pol and Robertson cases separately in a probe under the same limit, with C3's counter, to see which one suffers. Check orientation with C2's probe against the copy instead. **Done when.** You can say why a diagonal Jacobian cannot reveal a transpose, and why a stall is a poor way to notice it (F3).

## E. Capstone: write against the contracts

### E1. A second method: the trapezoidal rule
**Goal.** Add a method that `Stepper.fixed` can drive, without touching the library. **Hints.** `y_{n+1} = y_n + h/2 (f(t_n, y_n) + f(t_{n+1}, y_{n+1}))` fits the stage equation with `psi = y_n + h/2 f(t_n, y_n)` and `gamma = h/2`. Check second order like the `bdf2 logistic` line (error ratio near 4), then run it on the canary at `dt = 1e-3`, where `h λ = 10` on the fast component: its factor per step, `(1 - hλ/2) / (1 + hλ/2)`, is `-2/3`, so that component flips sign each step while it decays. Compare the error at `t = 1` with the `bdf2 canary` line. As a stretch make it an `Ode.Embedded` whose error estimate is its gap to `Bdf1.step` (`Bdf1.history` is abstract: ignore it with `_`) and pass it to `Adaptive.integrate` with `Halving`.

```ocaml
open Vstiff

module Trapezoid : Ode.Method = struct
  type history = unit

  let start = ()

  let step rhs h () (at : Ode.point) =
    let psi = Vec.axpy (h /. 2.) (rhs at.t at.y) at.y in
    Result.map (fun y -> (y, ())) (Stage.solve rhs { Stage.t = at.t +. h; gamma = h /. 2.; psi } at.y)
end

let () =
  let open Problems.Logistic in
  let error dt = Vec.norm_inf (Vec.sub (Result.get_ok (Stepper.fixed (module Trapezoid) ~dt problem)) (exact 5.)) in
  Printf.printf "error ratio for dt = 0.01 and 0.005: %.2f\n" (error 0.01 /. error 0.005)
```

**Look in.** [src/ode.mli](../src/ode.mli), [src/bdf1.ml](../src/bdf1.ml), [src/stage.ml](../src/stage.ml), [architecture.md](architecture.md). **Done when.** The ratio is near 4, the module has the signature `Ode.Method`, and you can say what `history` would hold for a multistep method. Making it a library module with an interface is a change to discuss first.

### E2. A second controller: separate failure counts
**Goal.** Write a controller against `Ode.Controller`. **Hints.** Keep the policy of `Halving` but report `stats = { accepted : int; too_large : int; solver_failures : int }`, telling the two `Ode.rejection` cases apart. Define the `stats` record inside your module, as `Halving` does, and give the module no signature: the compiler checks it against `Ode.Controller` where you pass it to `Adaptive.integrate`. Run it with `Adaptive.integrate (module Bdf2) (module Split)`; for a nonzero solver count try a large first step (`~dt0:1e4 ~dt_max:1e4`) on Robertson. A controller that scales the step by the error needs more than this contract offers (F5). **Look in.** [src/halving.ml](../src/halving.ml), [src/halving.mli](../src/halving.mli), [src/ode.mli](../src/ode.mli). **Done when.** The run reaches `t_end`, your counts add up to `Halving`'s on a run where no Newton solve fails, and you can say which contract items you did not need.

## F. Starter contributions

Most come from the known limitations in the [README](../README.md). Read [CONTRIBUTING.md](../CONTRIBUTING.md) first and run its checklist before handing anything over. *Discuss first* marks an item that adds or changes a corpus case: state its purpose in one sentence and agree it before writing code. *Design first* marks an item that changes signatures: write the design and have it reviewed before code.

- **F1. Verify a document.** Pick a page of `docs/`, run every snippet in the probe project, compare each number with the code, an expected file or a derivation, and fix or report what is wrong. *Look in* the page and [CONTRIBUTING.md](../CONTRIBUTING.md). *Done when* you have a list of claims checked, corrections made and claims you could not check.
- **F2. Pin van der Pol against a reference** (discuss first). Its accuracy is not checked. *Hints:* compute a reference outside vstiff with at least two independent solvers at tight tolerances, record how in [test/refs.ml](../test/refs.ml), measure the real error with a probe, and propose a case whose bound keeps a margin. *Look in* the `van_der_pol` group of [test/corpus.ml](../test/corpus.ml). *Done when* the line passes, fails when you break `Bdf2` in a scratch copy, and its purpose is one sentence.
- **F3. Turn a stall into a line** (discuss first). A transposed Jacobian stalls the run (D5). *Hints:* `Guard.budget` turns a runaway into a printed line; run the van der Pol and Robertson cases behind it (the `rhs` through `Guard.budget`, the case through `Guard.run`). Measure their real `rhs` calls with `Instrument.count` first and check the budget leaves a wide margin. *Look in* [test/guard.ml](../test/guard.ml), [test/corpus.ml](../test/corpus.ml). *Done when* both expected files are byte-identical and D5 prints a line instead of stalling.
- **F4. Reuse the Jacobian** (design first). `Newton.solve` calls `jac` at every iterate, `n + 1` calls of `rhs` each. *Hints:* build it once per solve and pass it down as a value (the library is pure). Convergence slows, so digits in corpus lines may move: find out why, never promote. *Look in* [src/newton.ml](../src/newton.ml), [src/stage.ml](../src/stage.ml). *Done when* C3 shows fewer `rhs` calls and every changed corpus line is explained.
- **F5. Scale the step by the error** (design first). Production controllers use a safety factor times `(tol / err)^(1/(p+1))`. *Hints:* `Controller.accepted` and `rejected` do not receive `err` today, so this changes the contract in [src/ode.mli](../src/ode.mli) and `Adaptive`; keep every ω at most 2 or redo the stability argument. *Look in* [numerics/04-bdf.md](numerics/04-bdf.md) (section 7), [numerics/05-step-control.md](numerics/05-step-control.md) (section 11), [src/adaptive.ml](../src/adaptive.ml). *Done when* a written design says which signatures change, why ω stays below `1 + √2` and which corpus lines may move.
- **F6. Separate `rtol` and `atol`** (design first). One weight `1 / (1 + |y_i|)` controls tiny components such as Robertson's `y2` only loosely. *Hints:* work out how a second tolerance reaches a controller through `Controller.init`. *Look in* [src/halving.ml](../src/halving.ml), [src/adaptive.mli](../src/adaptive.mli). *Done when* a written design shows the new signatures and what happens to the existing `~tol` calls.
- **F7. Dense output** (design first). The drivers return only the final state. *Hints:* design a way to observe the accepted points without breaking purity; the soak test compares results with `=`, so new fields must be deterministic and free of NaN. *Look in* [src/adaptive.ml](../src/adaptive.ml), [test/soak.ml](../test/soak.ml). *Done when* a written design says what the driver returns and how the soak test still compares results.
- **F8. Show Newton's iterations** (design first). `Instrument` counts `rhs` calls only. *Hints:* return the iteration count from `Newton.solve` and carry it through `Stage`; it changes several signatures, and a method has no place to report it today. *Look in* [src/newton.ml](../src/newton.ml), [src/stage.ml](../src/stage.ml), [src/ode.mli](../src/ode.mli). *Done when* a written design lists every signature that changes and where the count ends up.
- **F9. A cheaper canary pin** (discuss first). The `bdf1 canary` line takes 500,000 steps and dominates the test time. *Hints:* work out what that line catches, with the mutation table in [testing.md](testing.md), before proposing a faster case that catches the same defects; retiring the old line needs agreement. *Look in* the `backward_euler` group of [test/corpus.ml](../test/corpus.ml), [numerics/06-the-corpus.md](numerics/06-the-corpus.md) (section 4). *Done when* the new line fails under the mutations that fail the old one, and the discussion decided whether the old line stays.
