# Testing

vstiff is tested by two programs that print text and a file of expected text next to each: a test passes when the program prints exactly what the file says. This page explains the mechanics, the structure of the test code, the rules that keep the expected files meaningful, how to add a case, where reference values come from, how to run probes, and how to judge the strength of the tests by breaking the code on purpose. What each corpus case means mathematically is in [numerics/06-the-corpus.md](numerics/06-the-corpus.md). The OCaml needed to read the tests is in [ocaml.md](ocaml.md), the module structure in [architecture.md](architecture.md).

## What is in `test/`

| File | What it is |
|---|---|
| [test/dune](../test/dune) | declares the two test programs, `corpus` and `soak` |
| [test/corpus.ml](../test/corpus.ml), [test/corpus.expected](../test/corpus.expected) | the **corpus**: a table of cases, one printed line each, and the text it must print |
| [test/soak.ml](../test/soak.ml), [test/soak.expected](../test/soak.expected) | the **soak** test: four cases, ten rounds each, and its expected text |
| [test/problems.ml](../test/problems.ml) | the corpus problems: `Canary` (three independent decays at very different rates), `Logistic`, `VanDerPol`, `Robertson` |
| [test/refs.ml](../test/refs.ml) | reference values computed outside vstiff |
| [test/guard.ml](../test/guard.ml) | `Guard.run` and `Guard.budget`: the only test module that raises or catches |
| [test/report.ml](../test/report.ml) | `Report.lines`: the only module that prints |

`problems.ml`, `refs.ml`, `guard.ml` and `report.ml` are ordinary modules, not tests: dune compiles each once and links it into the programs that use it, so editing one can change an output. `corpus.ml` and `soak.ml` are main modules and cannot use each other (dune gives each an empty interface), which is why `soak.ml` carries its own copies of two small helpers.

## How an expect test works

`test/dune` declares `(tests (names corpus soak) (libraries vstiff))`. For each name dune builds `NAME.exe` from `NAME.ml`, runs it, captures standard output and compares it with `NAME.expected`. There are no assertions: a check fails by printing something other than the expected file says. The point of text: a failing check does not stop the run, so one run reports every line that changed, with the new value next to the old. The price is that the comparison is exact, which the rules below deal with.

- **Everything matches:** `dune runtest` is silent and exits with status 0. Dune remembers a pass, so a repeat with nothing changed does nothing; a failing comparison is reported again every time.
- **Output differs:** a diff and exit status 1 (next section).
- **A program crashes** with an uncaught exception: dune shows the exception and no diff.
- **A program never ends:** a case needs vastly more steps than it should. Interrupt it, and use probes to find which case.

| Goal | Command |
|---|---|
| run everything and compare | `dune runtest` |
| compare again although nothing changed | `dune runtest --force` (repeats the comparison, not the run) |
| run the programs again from scratch | `dune clean`, then `dune runtest` |
| see a program's raw output | `dune exec ./test/corpus.exe` (or `./test/soak.exe`) |
| compare by hand | `dune build`, then `./_build/default/test/corpus.exe \| diff test/corpus.expected -` (`<` lines are expected, `>` lines were printed) |
| time a program | `time ./_build/default/test/soak.exe` |

The recorded output of each program is `_build/default/test/corpus.exe.output` and `_build/default/test/soak.exe.output`.

### Reading a failure

Suppose `Jac.forward` builds its rows so that rows 0 and 2 trade places when the system has three components. `dune runtest` prints this (the `index` line with hashes, which follows the `diff` line, is left out):

```diff
File "test/corpus.expected", line 1, characters 0-0:
diff --git a/_build/default/test/corpus.expected b/_build/default/test/corpus.exe.output
--- a/_build/default/test/corpus.expected
+++ b/_build/default/test/corpus.exe.output
@@ -1,14 +1,14 @@
 newton linear 2d: Ok [2.000000000000; 3.000000000000]
 newton quadratic x0=1: Ok [1.414213562373]
 newton quadratic x0=0: Error Diverged
-jac canary at origin, max entry error < 1e-6: true
-jac canary at y0, max entry error relative to |J_ij| < 1e-6: true
-bdf1 canary t=1 dt=2e-6: max error 3.68e-07 < 1e-06: true
+jac canary at origin, max entry error < 1e-6: false
+jac canary at y0, max entry error relative to |J_ij| < 1e-6: false
+bdf1 canary t=1 dt=2e-6: max error 3.83e-07 < 1e-06: true
 bdf2 logistic [0,5]: error dt=0.01 3.116e-06, dt=0.005 7.810e-07, ratio 3.99 in [3.5, 4.5]: true
 vdp mu=1000 [0,2000] tol=1e-4: Ok at t=2000, rejected >= 1: true, y finite: true
 robertson t=1e4 tol=1e-6: |y1 - ref| < 1e-3: true, |y1 + y2 + y3 - 1| < 1e-8: true
-bdf2 canary t=1 dt=1e-3 (h lambda = 10): max error 1.53e-07 < 1e-06: true
-jac non-symmetric at (1, 2, 3), max entry error relative to |J_ij| < 1e-6: true
+bdf2 canary t=1 dt=1e-3 (h lambda = 10): Error Diverged
+jac non-symmetric at (1, 2, 3), max entry error relative to |J_ij| < 1e-6: false
 linalg zero leading pivot: [1.000000000000; 2.000000000000; 3.000000000000]
 linalg tiny leading pivot: [1.000000000000; 1.000000000000]
 adaptive blow-up y' = y^2 from y(0)=1 to t=2: Error StepRejected 1
```

The header names the two files compared: the expected file and what the program printed. Lines starting with `-` are in the expected file but were not printed, lines starting with `+` were printed but are not in the file, and unmarked lines are context. Ask first which lines changed and what they have in common: here the three lines about Jacobians and the stiff step that needs a good one. One printed number moved (3.68e-07 to 3.83e-07) although its criterion still holds: a root of the stage equation does not depend on the Jacobian, only the speed of Newton's iteration does, so the iteration stopped at a slightly different point. A number printed with three digits is a tight check, a boolean a loose one. The corpus runs from the lowest layer (Newton) to the highest (the adaptive integrator), so the first changed line usually points at the real problem.

A crash looks different: dune names the program in `test/dune`, where the caret marks `corpus` in the `names` field (the line number is wherever the field sits in `test/dune`), and prints the exception, here an array index out of bounds. A mistake deep in the library often crashes both programs, and dune then prints one such block for each:

```text
File "test/dune", line 5, characters 8-14:
5 |  (names corpus soak)
            ^^^^^^
Fatal error: exception Invalid_argument("index out of bounds")
```

## The case table, `Report` and `Guard`

[test/corpus.ml](../test/corpus.ml) is a table of cases, each a name and a function from `()` to the text printed after the name. The cases come in groups, one `let` per group, and the groups are joined in order:

```text
let newton = [ ("newton linear 2d", fun () -> show (Newton.solve ...)); ... ]
let corpus = List.concat [ newton; jacobian; backward_euler; ...; robertson_accuracy; newton_overflow; clock ]
let () = Report.lines corpus
```

`Report.lines` runs the cases in order and prints `name: text` for each. Computing the text is separate from printing it, so a case is a pure function and `Report` is the only module that prints. Reading the file from the top, the first six groups (Newton, Jacobian, backward Euler, BDF2 order, van der Pol, Robertson) run the library from the bottom layer to the top. The remaining groups are regression pins: each pins down specific mistakes, listed under "Judging test strength by mutation" below. The helpers at the top of the file (`pp_vec`, `failure`, `show`, `max_abs`, `max_error`, `error_below`, `bdf2_halving`) keep the cases short and the printed lines uniform.

`Guard.run ok compute` builds the function of a table entry so that every outcome a case can have on purpose becomes a line of text: `compute` returns a `result`, and `ok` renders an `Ok`. An `Error e` prints `Error` and the name of `e`; an `Invalid_argument m` prints `Invalid_argument m`; a right-hand side that exhausted its budget prints `no answer within 5e6 rhs calls`. `Guard.budget rhs` wraps a right-hand side (with `Instrument.count`) so that after 5,000,000 calls it raises instead of evaluating. The groups `give_up` and `clock` run their termination and argument-check cases through both, so a loop that never ends prints a line, and the diff shows it, instead of hanging the run.

## The rules

These keep the expected files worth trusting. Each has a reason.

1. **Never weaken or rewrite a passing expectation to get green.** The expected file records what correct output looks like. A run made green by editing the file, loosening a bound, deleting a case or skipping a test proves nothing. A changed line is a regression (fix the code) or a deliberate change of behaviour, which is explained in its commit message and reviewed as a diff of the `.expected` file.
2. **`dune promote` only records.** It is legitimate for the line of a newly added case, or for a deliberate, reviewed change of output format. It overwrites the `.expected` file with whatever the program printed, wrong answers included. Read the diff first, and read `git diff test/corpus.expected` afterwards.
3. **Reference values are never edited.** The values in [test/refs.ml](../test/refs.ml) are computed outside vstiff. A new reference is a new definition with its provenance; an existing one never changes to make a test pass.
4. **Every behaviour change comes with a test.** A new feature gets a case; a bug fix gets a case that failed before the fix.
5. **Cases are added deliberately, one at a time, each with a stated purpose** in a comment above it. A case costs run time and reading time.
6. **Corpus problems never supply an analytic Jacobian.** The solver must work for any black-box right-hand side and the corpus must keep exercising `Jac`. (The `newton` cases hand `Newton.solve` its own Jacobians, and the Jacobian cases compare `Jac.forward` with analytic ones; no integrator is ever given one.)
7. **Print stable lines.** No times, no counts that are not part of the contract, and no more digits than the check needs: a few (`%.3e`) or a boolean for a bound. The Newton lines print 12 decimals on purpose; that is still short of the 16 or 17 digits where results can differ between machines. Say in the line what is checked, as in `max error 3.68e-07 < 1e-06: true`.
8. **Effects are quarantined.** In the tests only `Guard` raises or catches and only `Report` prints; a case itself stays pure.
9. **A change that must not alter results leaves both expected files byte-identical.** Mind the floating-point trap in [architecture.md](architecture.md): moving a product into or out of a sum can change the last bits.

## Adding a case, step by step

The worked example is a case that does not exist yet: backward Euler is first order, so halving `dt` should roughly halve its error, and no line checks `Bdf1` at two step sizes.

**1. State the purpose in one sentence.** "Pin the order of backward Euler." Which mistake would this case catch that no existing line catches? If you cannot say, do not add the case. Raise it for discussion before you write it.

**2. Add the problem** to `test/problems.ml` unless one fits (here `Problems.Logistic` does): a module with `rhs`, `y0`, `problem` and, when a closed form exists, `exact`; no Jacobian. A reference value goes into `refs.ml` first (see "Reference values").

**3. Add a group to `test/corpus.ml`**, with a comment stating the purpose, above the `corpus` table, and append it to the table. The table's order is the order of the lines in `corpus.expected`, and a new line at the end leaves the line numbers that [numerics/06-the-corpus.md](numerics/06-the-corpus.md) cites for the existing ones unchanged:

```text
(* Backward Euler is first order: halving dt halves the error. *)
let bdf1_order =
  let open Problems.Logistic in
  let error dt = Result.map (max_error (exact 5.)) (Stepper.fixed (module Bdf1) ~dt problem) in
  [
    ( "bdf1 logistic [0,5]",
      fun () ->
        match (error 0.01, error 0.005) with
        | Ok coarse, Ok fine ->
            let ratio = coarse /. fine in
            Printf.sprintf "error dt=0.01 %.3e, dt=0.005 %.3e, ratio %.2f in [1.8, 2.2]: %b" coarse fine ratio
              (1.8 <= ratio && ratio <= 2.2)
        | Error e, _ | _, Error e -> failure e );
  ]
```

Use every helper you define: an unused top-level definition in the main module of a test program is a build error. Handle the `Error` case, so a failure prints a line instead of ending the run.

**4. Run `dune runtest`.** It fails, and the diff must show one `+` line and no `-` line:

```text
+bdf1 logistic [0,5]: error dt=0.01 <e1>, dt=0.005 <e2>, ratio <r> in [1.8, 2.2]: true
```

A `-` line means an existing line changed: stop, the new code affected something it should not have.

**5. Review the line as if someone else wrote it.** The ratio of a first-order method's errors at `dt` and `dt / 2` tends to 2 (a second-order method would give 4), so `<r>` should be close to 2. Is the bound loose enough to survive another machine and tight enough to fail when something is wrong? Say how you chose it.

**6. Record it.** Run `dune promote` from the repository root; it prints `Promoting _build/default/test/corpus.exe.output to test/corpus.expected.` Run `dune runtest` again: it must be silent, and `git diff test/corpus.expected` must show only your line.

**7. Check that the case can fail**, in a scratch copy (see "Judging test strength by mutation"). Break what it protects, for instance `gamma = h` to `gamma = h /. 2.` in `Bdf1.step`, and the new line must turn `false`. A case that cannot fail is decoration.

**8. Commit** the problem, the case and the one new expected line together as one logical change, such as `test(corpus): pin the first-order error of backward Euler`, and describe the case in [numerics/06-the-corpus.md](numerics/06-the-corpus.md). A case that also belongs in the soak test needs a `Case` module there and its line in `soak.expected`.

## Reference values

Where a problem has no closed form, the answer to compare with comes from outside vstiff. [test/refs.ml](../test/refs.ml) holds one: y1 of Robertson's problem at t = 1e4, `0.10730042854`, from SciPy's `solve_ivp` with Radau (rtol 1e-13), BDF and LSODA (rtol 1e-12), atol 1e-20 and finite-difference Jacobians, the three agreeing to 4e-12. Its doc comment is the model: an independent solver at tight tolerances, more than one method so that an error in one would show, and the provenance next to the value. With Python and SciPy installed, this script prints the value for each method:

```python
from scipy.integrate import solve_ivp

def robertson(t, y):
    a = 0.04 * y[0]
    b = 1e4 * y[1] * y[2]
    c = 3e7 * y[1] ** 2
    return [b - a, a - b - c, c]

for method, rtol in [("Radau", 1e-13), ("BDF", 1e-12), ("LSODA", 1e-12)]:
    sol = solve_ivp(robertson, (0.0, 1e4), [1.0, 0.0, 0.0], method=method, rtol=rtol, atol=1e-20)
    print(f"{method:6s} y1(1e4) = {sol.y[0, -1]:.11f}")
```

A new reference follows the same recipe and is added, never edited afterwards. Never compute a reference with vstiff: it could not catch that code's own mistakes. Choose the bound of a check from a measured error with a margin, and say how: a bound far looser than the observed error cannot catch a small regression. The corpus has both kinds of Robertson line, the loose `1e-3` and the tight `1e-5`; see [numerics/06-the-corpus.md](numerics/06-the-corpus.md).

## The soak test

[test/soak.ml](../test/soak.ml) runs four corpus cases ten times each in one process and prints one line per case, for instance `soak canary x10: passed 10/10, identical: true`. `passed` counts the rounds whose result satisfies the case's criterion (the same criterion as its corpus line) and `identical` is `true` when every round equals the first, compared with `=`. The code is a module type `Case` (a result type `t`, a `name`, `run : unit -> (t, Fail.t) result` and `pass : t -> bool`), a function `repeat (module C : Case)` that returns every round's result (a modular explicit, because its result type names `C.t`), and one module per case. To add a case, write a module with that signature and add `soak (module YourCase)` to the final list.

It is **a tripwire, not extra coverage.** The library has no hidden state, so `identical` can only turn `false` if someone adds hidden state, randomness or parallelism, or if a `nan` appears in a result (`nan` is never equal to itself, and `=` treats `0.` and `-0.` as equal, so "identical" is slightly weaker than bit for bit). The soak cases restate the corpus criteria in their own code, so changing a criterion means changing it in both places.

## Probes

A **probe** is a throwaway program that calls the library to answer a question: what does this return, how many steps does this case take, does that argument raise? Run one before you state a fact about behaviour in a comment or a document. Keep probes outside the repository: every `dune` file in the tree is part of the project, so a probe inside it is built for everyone. From the repository root:

```sh
REPO=$(pwd)
mkdir -p ../vstiff-scratch/p && cd ../vstiff-scratch
printf '(lang dune 3.0)\n' > dune-project
ln -sfn "$REPO/src" src
ln -sfn "$REPO/test/problems.ml" p/problems.ml
printf '(executable (name probe) (libraries vstiff))\n' > p/dune
```

(Add `ln -sfn "$REPO/test/refs.ml" p/refs.ml` for probes that use `Refs`.) Write `p/probe.ml`, starting with `open Vstiff`, then build and run it (`--root .` makes dune use the probe project even inside another dune project, and `dune exec --root . ./p/probe.exe` does both steps):

```sh
dune build --root . ./p/probe.exe && ./_build/default/p/probe.exe
```

A probe is the main module of an executable, so an unused top-level definition is a build error (warning 32): print what you define. A runaway probe belongs under a time limit, `perl -e 'alarm 120; exec @ARGV' ./_build/default/p/probe.exe`; a killed program prints only what it flushed, so end a `Printf.printf` format with `%!`. Two examples. Where does Newton evaluate the residual? Printing from inside `f` shows each point it tries:

```ocaml
open Vstiff

let () =
  let f x =
    Printf.printf "f called at x = %.17g\n" x.(0);
    [| (x.(0) *. x.(0)) -. 2. |]
  in
  let jac x = [| [| 2. *. x.(0) |] |] in
  match Newton.solve f jac [| 1. |] with
  | Ok x -> Printf.printf "result: %.17g\n" x.(0)
  | Error e -> print_endline (Fail.to_string e)
```

It prints the iterates of Newton's method for x² - 2 = 0 from x = 1, whose correct digits roughly double each time ([numerics/02-newton.md](numerics/02-newton.md)): `1`, `1.5`, `1.4166666666666667`, `1.4142156862745099`, `1.4142135623746899`, and then the result `1.4142135623730951`, the last iterate plus a final small step that is taken without evaluating `f` again. How many steps does a corpus case take? The counts are not in the expected files, so measure:

```ocaml
open Vstiff

let () =
  match Adaptive.integrate (module Bdf2) (module Halving) ~tol:1e-6 Problems.Robertson.problem with
  | Ok s -> Printf.printf "accepted %d, rejected %d\n" s.stats.accepted_steps s.stats.rejected_steps
  | Error e -> print_endline (Fail.to_string e)
```

## Judging test strength by mutation

A test suite is only as good as the mistakes it notices. You can measure that by hand: **break the code on purpose and see whether any line changes.** Silence after a plausible mistake means nothing checks that thing.

1. Work in a scratch copy of the committed files: from the repository root, `mkdir ../vstiff-mutation && git archive HEAD | tar -x -C ../vstiff-mutation`, then `cd ../vstiff-mutation` and run dune with `--root .` there.
2. Check that `dune runtest --root .` is silent.
3. Make one small, plausible mistake: a swapped argument, a changed constant, a removed check.
4. Run the tests. Classify: a **diff** (noticed; it shows which lines), a **crash**, a run that **never ends** (noticed clumsily), or **silence** (a gap).
5. Restore the file and try another.

Run what might not end under a time limit, building first and running the program directly so that nothing is left running: `dune build --root . ./test/corpus.exe`, then `perl -e 'alarm 120; exec @ARGV' ./_build/default/test/corpus.exe | diff test/corpus.expected -`. Output is printed when the program exits, so a killed run prints nothing.

The regression pins, the groups after the first six in `corpus.ml`, are worked examples: each row below can be repeated, and each pin fails when its mistake is made. "Silent" means no line of either expected file changed; [numerics/06-the-corpus.md](numerics/06-the-corpus.md) (sections 9 and 10) adds the line numbers and a few more changes.

| Mistake (one edit) | What notices it |
|---|---|
| `Jac.forward`: rows 0 and 2 swapped (three components) | the two canary Jacobian lines and the non-symmetric line (`false`), the stiff canary (`Error Diverged`); the backward Euler canary only changes digits |
| `Jac.forward`: the transpose | the non-symmetric line (`false`), but Robertson stops finishing, so a full run stalls before it reaches that line; van der Pol still finishes |
| `Jac.step`: `1e-14` for `1e-8` | the canary-at-y0 line and the non-symmetric line (`false`) |
| `Jac.step`: `1e-2` for `1e-8` | the non-symmetric line (`false`); the digits of the Robertson accuracy line move |
| `Linalg.pivot`: always row 0 | the two `linalg` lines |
| `Stage`: `I + γJ` for `I - γJ`; `Newton`: `max_iter = 1`; `Halving`: never double the step | the run does not finish (killed by the time limit) |
| `Bdf2.coeffs`: the sign of any one coefficient flipped | the BDF2 order line, the stiff canary, the Robertson lines and most adaptive lines |
| `Bdf2`: return backward Euler's result instead of BDF2's | the Robertson accuracy line (the loose `1e-3` line still passes), the `dt0 = 0.5` line and the blow-up count |
| `Bdf2`: start-up estimate factor `1e-3` for `0.5` | the `dt0 = 0.5` line |
| `Newton`: tolerance `1e-6` for `1e-10` | only `newton quadratic x0=1`, in its 12th decimal |
| `Newton`: the step test applied only after the line search accepts a step, not before it | the backward Euler canary, the BDF2 order line and the stiff canary print `Error Diverged`; van der Pol, the `dt0 = dt_max = 1e30` line and both logistic count lines print `Error StepRejected`, the blow-up and NaN-wall lines other counts; the overflow line prints `max_float` in full and the `from t=1e15` line `Error StepRejected 1`; `soak canary`, `soak logistic order` and `soak van der Pol` pass 0/10 |
| `Newton`: no finiteness check on the converged `x + dx` | the `newton step that converges into an overflow` line prints `Ok [inf]` |
| `Halving`: double after two accepts, not three | both logistic count lines and the counts on the blow-up line (`StepRejected 2`) and the NaN-wall line (the right-hand side that returns `nan` past `t = 0.5`) |
| `Halving`: no `dt_max` cap on doubling | the `dt_max = 1e-3` line |
| `Halving`: no `16 eps abs(t)` floor | the blow-up, NaN-wall and `from t=1e10` lines print `Error StepRejected 51` instead of 1, 46 and 1: only `max_rejects` is left to end the run |
| `Halving`: halve the proposal, not the step that failed | the `dt0 = dt_max = 1e30` line (`Error StepRejected 51`) |
| `Halving`: accept every step | the two Robertson lines, the `dt0 = 0.5` line, the blow-up and NaN-wall counts and `soak robertson` (it stops passing) |
| `Check`: remove the argument checks | the three `Invalid_argument` lines |
| `Adaptive`: do not shorten the last step (`h = dt` where it takes `t_end - t`) | the `dt0 = dt_max = 1e30` line, the Robertson accuracy line, and the `over one ulp` and `from t=1e15` lines, whose state advances by `dt` instead of the remainder |
| `Adaptive`: no snapping (`h = dt` for a step that is not the last) | the `from t=1e10` line prints the budget text and the `from t=1e15` line `y = 76.84` instead of `y = 100` |
| `Adaptive`: no remainder rule (`last = dt >= remaining`) | the `over one ulp` line prints `Error StepRejected 1` |
| `Adaptive`: neither snapping nor the remainder rule | the `over one ulp` and `from t=1e10` lines print the budget text and the `from t=1e15` line `y = 76` |
| `Bdf2`: start Newton from `y_n` instead of the extrapolation | silent: only the speed of Newton changes |
| `Newton`: Armijo constant `0.`; `min_damping = 1.`; no check that the step, the line-search residual or the iterate is finite | each silent, and so is dropping the line search altogether (always taking the full step) |
| `Adaptive`: default `dt_max` of `span` for `span / 10`; no finiteness check on `y0` and `rhs t0 y0` | each silent |
| `Stepper.fixed`: drop the `max 1` | silent |
| `Bdf1`, `Bdf2`: evaluate the stage at `at.t` instead of `at.t + h` | silent |

The corpus catches gross breakage (a wrong coefficient, a Jacobian that is plainly wrong, a step that is never allowed to grow back, which makes the run not finish), the pins catch particular mistakes, and it is nearly blind to anything that only changes how fast Newton converges.

### Known gaps

The silent rows are gaps, and starting points for contributions; [numerics/06-the-corpus.md](numerics/06-the-corpus.md) (section 10) suggests a case for most of them.

- **Newton's safeguards.** The line search (Armijo constant, minimum damping) and the non-finite checks at the start of an iteration, on the step and on the trial residual can be removed without changing any line (the check on a converged `x + dx` is pinned): `newton quadratic x0=1` converges with full steps, so no case needs damping. The order of the two tests is pinned, indirectly (table above); a direct pin would be `x² - 3` from 1, which converges today and gives `Error Diverged` when the line search comes first (check with a probe).
- **The time argument of `rhs`.** The four corpus problems are autonomous (their right-hand sides ignore `t`), and the NaN wall, the one right-hand side that reads `t`, does not notice a stage evaluated at the wrong time.
- **Defaults and guards of the drivers:** the default `dt_max`, the finiteness check at the start, `max 1` in `Stepper.fixed`.
- **The `Too_small` rejection.** Snapping alone gives the same lines: without the rejection the method is called with `h = 0`, the second such call fails on `ω = 0 / 0`, and the controller gives up with the same `StepRejected 1`. A case would wrap `Bdf2` and record the `h` of every call ([numerics/06-the-corpus.md](numerics/06-the-corpus.md), section 10).
- **van der Pol** checks only `Ok` at `t = 2000`, a rejection and a finite state, not the values against a reference.
- **A transposed Jacobian** is pinned by one line, but the mistake stalls the run before that line is reached.
- **The soak test** is a tripwire. **Run time:** the backward Euler canary takes (t_end - t0) / dt = 1 / 2e-6 = 500,000 fixed steps, once in the corpus and ten times in the soak test, and dominates the run; measure it with `time`.

Closing a gap is a self-contained contribution: write the case that would have noticed the mistake, add it by the steps above, and confirm with the mutation that it now fails. [exercises.md](exercises.md) turns several into starter contributions.
