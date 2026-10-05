# Testing

vstiff is tested by two programs that print text and a file of expected text next to each: a test passes when the program prints exactly what the file says. This page explains the mechanics, the structure of the test code, what the rules about the expected files protect, how to add a case, where reference values come from, how to run probes, and how to judge the strength of the tests by breaking the code on purpose. The rules themselves are in [AGENTS.md](../AGENTS.md): H1 to H4, H6 and H7, and its "Tests" section. What each corpus case means mathematically is in [numerics/06-the-corpus.md](numerics/06-the-corpus.md). [ocaml.md](ocaml.md) explains the OCaml the tests use and [architecture.md](architecture.md) the modules they call.

## What is in `test/`

| File | What it is |
|---|---|
| [test/dune](../test/dune) | declares the two test programs, `corpus` and `soak`, linked against both libraries |
| [test/corpus.ml](../test/corpus.ml), [test/corpus.expected](../test/corpus.expected) | the **corpus**: a table of cases, one printed line each, and the text it must print |
| [test/soak.ml](../test/soak.ml), [test/soak.expected](../test/soak.expected) | the **soak** test: four cases, ten rounds each, and its expected text |
| [test/problems.ml](../test/problems.ml) | the corpus problems: `Canary` (three independent decays at very different rates), `Logistic`, `VanDerPol`, `Robertson` |
| [test/refs.ml](../test/refs.ml) | a reference value computed outside vstiff; the file is never edited, and a new reference goes into a new module (H4) |
| [test/guard.ml](../test/guard.ml) | `Guard.run`, `Guard.budget` and `Guard.bounded`: the only test module that raises or catches |
| [test/report.ml](../test/report.ml) | `Report.lines`: the only module that prints |

`problems.ml`, `refs.ml`, `guard.ml` and `report.ml` are ordinary modules, not tests: dune compiles each once and links it into the programs that use it, so editing one can change an output. `corpus.ml` and `soak.ml` are main modules and cannot use each other (dune gives each an empty interface), which is why `soak.ml` carries its own copies of three small helpers (`max_error`, `bdf2_halving` and `on_budget`). Both start with `open Vstiff` and then `open Numerics`, because the cases call the solver and the kernel (`Newton`, `Jac`, `Linalg`, `Vec`) directly.

## How an expect test works

`test/dune` declares `(tests (names corpus soak) (libraries vstiff numerics))`. For each name dune builds `NAME.exe` from `NAME.ml`, runs it, captures standard output and compares it with `NAME.expected`. There are no assertions: a check fails by printing something other than the expected file says. The point of text: a failing check does not stop the run, so one run reports every line that changed, with the new value next to the old. The price is that the comparison is exact, so a line must print only what its check needs ("The rules behind the expected files", below).

- **Everything matches:** `dune runtest` is silent and exits with status 0. Dune remembers a pass, so a repeat with nothing changed does nothing; a failing comparison is reported again every time.
- **Output differs:** a diff and exit status 1 (next section).
- **A program crashes** with an uncaught exception: dune shows the exception and no diff.
- **A program never ends:** a case needs vastly more steps than it should. The long runs of the corpus are on a call budget and print a line instead (below); for the others (the two logistic count lines and the adaptive $y' = 2t$ line) interrupt the program and use probes to find which case. The soak test's adaptive cases are on the budget too, and its fixed-step cases always end.

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
@@ -40,7 +40,7 @@ newton residual nan at the start: Error Nan
 newton step that overflows in the linear solve: Error Diverged
 jac forward difference of y^3 at y = -2 (exact 12): 11.9999998224
 adaptive y' = 2t on [0, 1] tol=1e-6: max error 1.57e-12
-adaptive canary t=1 tol=1e-6: max error 4.66e-07 < 1e-6: true, rhs calls 106445
+adaptive canary t=1 tol=1e-6: max error 4.86e-08 < 1e-6: true, rhs calls 634048
 adaptive y0 = nan: Error Nan
 newton x^2 = 2 from x0=1e-3 (deep line search): Ok [1.414213562373]
 newton x^3 = 0 from x0=1 (needs more than 50 iterations): Error Diverged
```

The header names the two files compared: the expected file and what the program printed. Lines starting with `-` are in the expected file but were not printed, lines starting with `+` were printed but are not in the file, and unmarked lines are context. Ask first which lines changed and what they have in common: here the three lines about Jacobians, the stiff step that needs a good one and, far below, the adaptive canary. Two printed errors moved although their criteria still hold ($3.68 \times 10^{-7}$ to $3.83 \times 10^{-7}$, and $4.66 \times 10^{-7}$ to $4.86 \times 10^{-8}$ on the adaptive canary): a root of the stage equation does not depend on the Jacobian, only the speed of Newton's iteration does, so the iteration stopped at a slightly different point, and the adaptive canary shows the speed itself in its count of right-hand-side calls (106445 to 634048). A number printed with three digits is a tight check, a boolean a loose one. The corpus runs from the lowest layer (Newton) to the highest (the adaptive integrator), so the first changed line usually points at the real problem.

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
let corpus = List.concat [ newton; jacobian; backward_euler; ...; clock; too_small; fixed_clock; ...; adaptive_time; adaptive_canary; start_and_limits ]
let () = Report.lines corpus
```

`Report.lines` runs the cases in order and prints `name: text` for each. Computing the text is separate from printing it, so a case is a pure function and `Report` is the only module that prints. Reading the file from the top, the first six groups (Newton, Jacobian, backward Euler, BDF2 order, van der Pol, Robertson) run the library from the bottom layer to the top. The remaining groups are regression pins, 38 of the 47 lines: each pins down specific mistakes, listed under "Judging test strength by mutation" below. The helpers at the top of the file (`pp_vec`, `failure`, `show`, `max_abs`, `max_error`, `bdf2_halving`, `on_budget`, `error_below`) keep the cases short and the printed lines uniform.

`Guard.run ok compute` builds the function of a table entry so that every outcome a case can have on purpose becomes a line of text: `compute` returns a `result`, and `ok` renders an `Ok`. An `Error e` prints `Error` and the name of `e`; an `Invalid_argument m` prints `Invalid_argument m`; a right-hand side that exhausted its budget prints `no answer within 5e6 rhs calls`. `Guard.budget rhs` wraps a right-hand side (with `Instrument.count`) so that after 5,000,000 calls it raises instead of evaluating. The groups `give_up`, `clock`, `too_small`, `zero_floor`, `arguments` and `adaptive_canary` run their cases through both, and `on_budget` puts the van der Pol and Robertson lines on the budget too, so a loop that never ends, or a mistake that makes a long run crawl, prints a line, and the diff shows it, instead of stalling the run: a transposed Jacobian ends the corpus in about 2 s. `Guard.bounded run` is the same idea for the soak test: it returns `Some` of what `run ()` returns, or `None` when a budget inside it ran out, so a round can be counted as not passed without any module but `Guard` catching an exception.

## The rules behind the expected files

The rules for changing the tests are in [AGENTS.md](../AGENTS.md). Each paragraph below says what one of them protects in this suite, and what to do when the situation it covers comes up.

**The expected files record what correct output looks like (H2).** A run made green by editing a file, loosening a bound, deleting a case or skipping a test proves nothing. So a changed line is one of three things: a regression, which is fixed in the code; a deliberate change of algorithm, which re-pins the lines it moves (counts, digits) in the same commit, with the old and the new values in its message and the `.expected` diff reviewed; or a line that printed more than its check needs (H1), which is rewritten to print the check, in a commit of its own. The correctness bounds do not loosen in any case.

**`dune promote` only records.** It overwrites an `.expected` file with whatever the program printed, wrong answers included. It is for the line of a newly added case, for a deliberate re-pin and for a line rewritten to print its check, never for a failure. Read the diff before promoting, and `git diff test/corpus.expected` afterwards.

**A line prints the check, not the bits (H1).** A refactor may change the last bits of a result: fused multiply-add, another order of summation, another platform ([architecture.md](architecture.md), "Performance"). No check may depend on those bits. A line is therefore a verdict on a bound, a typed outcome, or a value printed to the digits its bound makes meaningful (a few, as in `%.3e`, or a boolean), and the same files pass on every platform. The Newton lines print 12 decimals on purpose; that is still short of the 16 or 17 digits where results can differ between machines. A line says what it checks, as in `max error 3.68e-07 < 1e-06: true`, and prints no times and no counts that are not part of the contract. A change that moves a line either changed behaviour, or exposed a line printed more precisely than its check needs: in the second case, make the line print the check, in a commit of its own whose message gives the old and the new line (H2). What stays exact is determinism within one build: the same binary returns the same output every time, which the soak test checks.

**The corpus exercises the solver the way users will (H7).** The solver must work for any black-box right-hand side, and the corpus must keep exercising `Jac`, so no integrator is ever given an analytic Jacobian. (The `newton` cases hand `Newton.solve` its own Jacobians, and the Jacobian cases compare `Jac.forward` with analytic ones; neither is an integrator.)

**Effects stay in `Guard` and `Report` (H6).** In the tests only `Guard` raises or catches and only `Report` prints; a case itself stays pure ("The case table, `Report` and `Guard`", above).

Where reference values come from (H4) is under "Reference values", and adding a case (H3 and "Tests" in AGENTS.md) is the next section.

## Adding a case, step by step

[AGENTS.md](../AGENTS.md) lists the steps under "Tests", "Adding a corpus case". This page walks through them with an example, and adds the check that the case can fail (H3) and the commit. The example is a case that does not exist yet: backward Euler is first order, so halving `dt` should roughly halve its error, and no line checks `Bdf1` at two step sizes.

**1. State the purpose in one sentence.** "Pin the order of backward Euler." Which mistake would this case catch that no existing line catches? If you cannot say, do not add the case: it would cost run time and reading time and protect nothing new.

**2. Add the problem** to `test/problems.ml` unless one fits (here `Problems.Logistic` does): a module with `rhs`, `y0`, `problem` and, when a closed form exists, `exact`; no Jacobian (H7). A reference value is computed outside vstiff and goes into a new module first (H4; see "Reference values").

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

Use every helper you define: an unused top-level definition in the main module of a test program is a build error. Handle the `Error` case, so a failure prints a line instead of ending the run. This case takes fixed steps and always ends; a run that might not finish, an adaptive one say, goes on the call budget: its right-hand side in `Guard.budget` and the case in `Guard.run` ("The case table, `Report` and `Guard`", above).

**4. Run `dune runtest`.** It fails, and the diff must show one `+` line and no `-` line:

```text
+bdf1 logistic [0,5]: error dt=0.01 <e1>, dt=0.005 <e2>, ratio <r> in [1.8, 2.2]: true
```

A `-` line means an existing line changed: stop, the new code affected something it should not have.

**5. Review the line as if someone else wrote it.** The ratio of a first-order method's errors at `dt` and `dt / 2` tends to 2 (a second-order method would give 4), so `<r>` should be close to 2. Is the bound loose enough to survive another machine and tight enough to fail when something is wrong? Say how you chose it.

**6. Record it.** Run `dune promote` from the repository root; it prints `Promoting _build/default/test/corpus.exe.output to test/corpus.expected.` Run `dune runtest` again: it must be silent, and `git diff test/corpus.expected` must show only your line.

**7. Check that the case can fail** (H3), in a scratch copy (see "Judging test strength by mutation"). Break what it protects, for instance $\gamma = h$ to $\gamma = h / 2$ in `Bdf1.step`, and the new line must turn `false`. A case that cannot fail is decoration.

**8. Commit** the problem, the case and the one new expected line together as one logical change ([AGENTS.md](../AGENTS.md), "Commits"), such as `test(corpus): pin the first-order error of backward Euler`, and describe the case in [numerics/06-the-corpus.md](numerics/06-the-corpus.md). A case that also belongs in the soak test needs a `Case` module there and its line in `soak.expected`.

## Reference values

Where a problem has no closed form, the answer to compare with comes from outside vstiff (H4). [test/refs.ml](../test/refs.ml) holds one: $y_1$ of Robertson's problem at $t = 10^4$, `0.10730042854`, from SciPy's `solve_ivp` with Radau ($\mathrm{rtol} = 10^{-13}$), BDF and LSODA ($\mathrm{rtol} = 10^{-12}$), $\mathrm{atol} = 10^{-20}$ and finite-difference Jacobians, the three agreeing to $4 \times 10^{-12}$. Its doc comment is the model for the provenance to keep next to a value: the solver and its version, the tolerances and how well the methods agree. With Python and SciPy installed, this script prints the value for each method:

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

`refs.ml` is never edited. A new reference is computed by two independent solvers that agree on every digit kept, and goes into a new module with the solver, version, tolerances and agreement recorded next to it (H4). Never compute a reference with vstiff: it could not catch that code's own mistakes. Choose the bound of a check from a measured error with a margin, and say how: a bound far looser than the observed error cannot catch a small regression. The corpus has both kinds of Robertson line, the loose $10^{-3}$ and the tight $10^{-5}$; see [numerics/06-the-corpus.md](numerics/06-the-corpus.md).

## The soak test

[test/soak.ml](../test/soak.ml) runs four corpus cases ten times each in one process and prints one line per case, for instance `soak canary x10: passed 10/10, identical: true`. `passed` counts the rounds whose result satisfies the case's criterion (the same criterion as its corpus line) and `identical` is `true` when every round equals the first, compared with `=`. The code is a module type `Case` (a result type `t`, a `name`, `run : unit -> (t, Fail.t) result` and `pass : t -> bool`), a function `repeat (module C : Case)` that returns every round's result, `None` for a round whose call budget ran out (a modular explicit, because its result type, `(C.t, Fail.t) result option list`, names `C.t`), and one module per case; `passed` counts the rounds that are `Some (Ok r)` with `C.pass r`. To add a case, write a module with that signature and add `soak (module YourCase)` to the final list.

It is **a tripwire, not extra coverage**, and what it watches is determinism within one build (H1). The library has no hidden state, so `identical` can only turn `false` if someone adds hidden state, randomness or parallelism, or if a `nan` appears in a result (`nan` is never equal to itself, and `=` treats `0.` and `-0.` as equal, so "identical" is slightly weaker than comparing every bit of the results). The soak cases restate the corpus criteria in their own code, so changing a criterion means changing it in both places. The adaptive cases, van der Pol and Robertson, run on the same call budget as the corpus lines (`on_budget` wraps their `rhs` with `Guard.budget`) and `repeat` runs every round through `Guard.bounded`, so a mistake that makes one of them crawl ends the round instead of stalling the soak test: with a transposed Jacobian the soak test ends in about 14 s and prints `soak robertson x10: passed 0/10, identical: true`, ten rounds that all ran out of budget (`None` equals `None`, which is why `identical` stays `true`). The fixed-step cases, the canary and the logistic order, have no budget and always end.

## Probes

A **probe** is a throwaway program that calls the library to answer a question: what does this return, how many steps does this case take, does that argument raise? Run one before you state a fact about behaviour in a comment or a document (H12). Keep probes outside the repository: every `dune` file in the tree is part of the project, so a probe inside it is built for everyone. From the repository root:

```sh
REPO=$(pwd)
mkdir -p ../vstiff-scratch/p && cd ../vstiff-scratch
printf '(lang dune 3.0)\n' > dune-project
ln -sfn "$REPO/src" src
ln -sfn "$REPO/test/problems.ml" p/problems.ml
printf '(executable (name probe) (libraries vstiff))\n' > p/dune
```

(The link to `src` brings both libraries, since `src/numerics` is inside it. Add `ln -sfn "$REPO/test/refs.ml" p/refs.ml` for probes that use `Refs`.) Write `p/probe.ml`, starting with the `open`s it needs, then build and run it (`--root .` makes dune use the probe project even inside another dune project, and `dune exec --root . ./p/probe.exe` does both steps):

```sh
dune build --root . ./p/probe.exe && ./_build/default/p/probe.exe
```

A probe that uses only the solver (`Adaptive`, `Stepper`, `Bdf1`, `Bdf2`, `Halving`, `Ode`, `Fail`, `Clock`, `Instrument`) needs just `(libraries vstiff)` and `open Vstiff`, as above. One that calls `Newton.solve`, `Jac.forward`, `Linalg.solve` or `Vec` needs `(libraries vstiff numerics)` in `p/dune` and `open Numerics` as well, as `test/corpus.ml` has it; an `open` that names nothing is itself a build error (warning 33). `Stage` and `Check` cannot be named from a probe at all (`Unbound module Vstiff.Stage`); to count the right-hand-side calls of one implicit step, call `Bdf1.step`. A probe is the main module of an executable, so an unused top-level definition is a build error (warning 32): print what you define. A runaway probe belongs under a time limit, `perl -e 'alarm 120; exec @ARGV' ./_build/default/p/probe.exe`; a killed program prints only what it flushed, so end a `Printf.printf` format with `%!`. Two examples. Where does Newton evaluate the residual? Printing from inside `f` shows each point it tries (a probe of the kernel: `numerics` goes into `p/dune`):

```ocaml
open Numerics

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

It prints the iterates of Newton's method for $x^2 - 2 = 0$ from $x = 1$, whose correct digits roughly double each time ([numerics/02-newton.md](numerics/02-newton.md)): `1`, `1.5`, `1.4166666666666667`, `1.4142156862745099`, `1.4142135623746899`, and then the result `1.4142135623730951`, the last iterate plus a final small step that is taken without evaluating `f` again. How many steps does a corpus case take? The counts are not in the expected files, so measure:

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

The regression pins, the groups after the first six in `corpus.ml`, are worked examples: each row below can be repeated, and each pin fails when its mistake is made. "Silent" means no line of either expected file changed; [numerics/06-the-corpus.md](numerics/06-the-corpus.md) (sections 9 and 10) adds the line numbers and a few more changes. Some lines are named by a word of their text, typeset where it is mathematics: $\sqrt{1 - t}$, `max_rejects = 5000` (the NaN wall at $t = 0$), `tol = 0`, `wrong length`, `span of 1e-320`, `atan`, `residual nan`, `overflows in the linear solve`, $y^3$, $y' = 2t$ (the fixed-step line, or the adaptive one where the row says so), `adaptive canary` (line 43) and, for lines 45 to 47, `deep line search`, `more than 50 iterations` and `leaves the domain`. `y0 = nan` is `bdf1 y0 = nan` (line 36) or `adaptive y0 = nan` (line 44).

| Mistake (one edit) | What notices it |
|---|---|
| `Jac.forward`: rows 0 and 2 swapped (three components) | the two canary Jacobian lines and the non-symmetric line (`false`), the stiff canary (`Error Diverged`) and the adaptive canary (`rhs calls 634048` instead of 106445); the backward Euler canary only changes digits |
| `Jac.forward`: rows 0 and 1, or rows 1 and 2, swapped | the same two canary lines and the non-symmetric line, the adaptive canary (`rhs calls 186690` or `551696`) and both Robertson lines, which print the budget text (`soak robertson` passes 0/10); the stiff canary only changes digits for rows 0 and 1 and prints `Error Diverged` for rows 1 and 2; the backward Euler canary does not change |
| `Jac.forward`: the transpose | the non-symmetric line (`false`) and both Robertson lines, which print the budget text: the corpus ends in about 2 s and van der Pol still passes; `soak robertson` passes 0/10, and the soak test ends in about 14 s |
| `Jac.step`: $10^{-14}$ for $10^{-8}$ | the canary-at-y0 line, the non-symmetric line (`false`), the $y^3$ line and the adaptive canary (`rhs calls 123846`) |
| `Jac.step`: $10^{-2}$ for $10^{-8}$ | the non-symmetric line (`false`), the $y^3$ line and both Robertson lines (the budget text); `soak robertson` passes 0/10 |
| `Jac.step`: no `abs`, $10^{-6}$ for $10^{-8}$, or a backward difference; `Jac.forward`: the nominal step as divisor instead of the stored one | only the $y^3$ line, in its digits |
| `Linalg.pivot`: always row 0 | the two `linalg` lines |
| `Stage`: $I + \gamma J$ for $I - \gamma J$ | the stiff canary (`Error Diverged`), van der Pol and both Robertson lines (the budget text), the adaptive canary (`rhs calls 872996`), a digit of the BDF2 order line; `soak van der Pol` and `soak robertson` pass 0/10 |
| `Newton`: `max_iter = 1` | van der Pol, both Robertson lines, the adaptive canary and the blow-up, NaN-wall and `dt0 = dt_max = 1e30` lines print the budget text, but both logistic count lines and the adaptive $y' = 2t$ line have no budget and run on: the corpus does not finish (killed by the time limit); the soak test ends, with all four lines at 0/10 |
| `Halving`: never double the step | van der Pol, both Robertson lines, the blow-up line and the adaptive canary print the budget text; both logistic count lines and the adaptive $y' = 2t$ line change; `soak van der Pol` and `soak robertson` pass 0/10 |
| `Bdf2.coeffs`: the sign of any one coefficient flipped | both BDF2 order lines (logistic and $y' = 2t$), the stiff canary, the Robertson lines, the adaptive canary and most other adaptive lines; `soak logistic order`, `soak van der Pol` and `soak robertson` pass 0/10 |
| `Bdf2`: return backward Euler's result instead of BDF2's | the Robertson accuracy line (the loose $10^{-3}$ line still passes), the `dt0 = 0.5` line, the blow-up count, the adaptive $y' = 2t$ line and the adaptive canary (`false`, error `3.45e-04`) |
| `Bdf2`: startup estimate factor $10^{-3}$ for $0.5$ | the `dt0 = 0.5` line and the adaptive canary's count (`rhs calls 106188`) |
| `Bdf2`, `Bdf1`: evaluate the stage at `at.t` instead of `at.t + h` | `Bdf2`: the $y' = 2t$ line (error ratio 1.75 instead of 4.00) and the adaptive one (`1.02e-03` instead of `1.57e-12`); `Bdf1`: the $\sqrt{1 - t}$, `max_rejects = 5000` and adaptive $y' = 2t$ lines |
| `Newton`: tolerance $10^{-6}$ for $10^{-10}$ | five lines: `newton quadratic x0=1` in its 12th decimal, the adaptive canary's count (`rhs calls 89820`), the 12th decimal of the $x^2 = 2$ and $\ln x = 0$ lines, and `Ok [0.000001545213]` instead of `Error Diverged` for $x^3 = 0$ |
| `Newton`: the step test applied only after the line search accepts a step, not before it | the backward Euler canary, the BDF2 order line, the stiff canary and the two fixed-step $y'$ lines print `Error Diverged`; van der Pol, the `dt0 = dt_max = 1e30` line, both logistic count lines and the adaptive $y' = 2t$ line print `Error StepRejected`, the blow-up and NaN-wall lines other counts, and the adaptive canary the budget text; the overflow line prints `max_float` in full and the `from t=1e15` line `Error StepRejected 1`; `soak canary`, `soak logistic order` and `soak van der Pol` pass 0/10 |
| `Newton`: no finiteness check on the converged `x + dx` | the `newton step that converges into an overflow` line prints `Ok [inf]` |
| `Newton`: no damping (the full step every time) | the `atan` line prints `Error Diverged` and the $\ln x = 0$ line `Error Nan` |
| `Newton`: a line search that accepts any finite trial point | the `atan` line prints `Error Diverged` |
| `Newton`: a `nan` residual at the start named `Diverged`, or a non-finite step named `Nan` | the `residual nan` line, respectively the `overflows in the linear solve` line |
| `Newton`: `min_damping = 1. /. 2.` | the `deep line search` line prints `Error Diverged` |
| `Newton`: an iteration limit of 100 for 50 | the `more than 50 iterations` line prints `Ok [0.000000000138]` |
| `Newton`: a `nan` trial residual accepted | the `leaves the domain` line prints `Error Nan` |
| `Halving`: double after two accepts, not three | both logistic count lines, the adaptive $y' = 2t$ line, the adaptive canary (`rhs calls 119661`) and the counts on the blow-up line (`StepRejected 2`) and the NaN-wall line (`StepRejected 47`) |
| `Halving`: no `dt_max` cap on doubling | the `dt_max = 1e-3` line |
| `Halving`: no $16\thinspace\varepsilon\thinspace\lvert t \rvert$ floor | the blow-up and NaN-wall lines print `Error StepRejected 7` and `Error StepRejected 51` instead of 1 and 46: the NaN wall is left to `max_rejects`, and the blow-up run ends when a halved step snaps to 0, after the driver's `Too_small` rejections of the steps below the resolution (the `from t=1e10` lines still end at once, because the halved step is 0) |
| `Halving`: no give-up on a halved step of 0 | the `max_rejects = 5000` line prints `Error StepRejected 5001` instead of 1055 |
| `Halving`: halve the proposal, not the step that failed | the `dt0 = dt_max = 1e30` line (`Error StepRejected 51`) |
| `Halving`: accept every step | the two Robertson lines, the `dt0 = 0.5` line, the blow-up and NaN-wall counts, the adaptive $y' = 2t$ line, the adaptive canary (`false`) and `soak robertson` (it stops passing) |
| `Check`: remove the argument checks | the three `Invalid_argument` lines of the termination group, and the `dt = 1e-7 at t=1e10`, `tol = 0` and two `wrong length` lines, which print `Ok`, a `StepRejected` or an index error instead |
| `Stepper.fixed`: a running sum for $t$ instead of the grid | the $\sqrt{1 - t}$ line prints `Error Nan` |
| `Stepper.fixed`: no finiteness check on `y0` and `rhs t0 y0` | the `bdf1 y0 = nan` line prints `Ok [nan]` |
| `Adaptive`: no fallback to the span when a default `dt0` or `dt_max` underflows to 0 | the `span of 1e-320` line prints `Invalid_argument Adaptive.integrate: dt0 and dt_max must be positive` |
| `Adaptive`: do not shorten the last step (`h = dt` where it takes `t_end - t`) | the `dt0 = dt_max = 1e30` line, the Robertson accuracy line, the adaptive $y' = 2t$ line, the adaptive canary (`false`), and the `over one ulp` and `from t=1e15` lines, whose state advances by `dt` instead of the remainder |
| `Adaptive`: no snapping (`h = dt` for a step that is not the last) | the `from t=1e15` line prints `y = 98.7` instead of `y = 100`, and the adaptive $y' = 2t$ line (`1.54e-12`) and the adaptive canary's count (`rhs calls 106441`) change; the `from t=1e10` and `never reaches the method` lines no longer notice, because the driver rejects those steps as below the resolution anyway |
| `Adaptive`: no remainder rule (`last = dt >= remaining`) | the `over one ulp` line prints `Error StepRejected 1` |
| `Adaptive`: neither snapping nor the remainder rule | the `over one ulp` line prints `Error StepRejected 1`, the `from t=1e15` line `y = 98.7`, and the adaptive $y' = 2t$ line and the adaptive canary's count change |
| `Adaptive`: no `Too_small` rejection (call the method with `h = 0`) | the `never reaches the method` line prints `rhs calls: 5` instead of `1` |
| `Adaptive`: no finiteness check on `y0` and `rhs t0 y0` | the `adaptive y0 = nan` line prints `Error StepRejected 51` instead of `Error Nan` |
| `Bdf2`: start Newton from $y_n$ instead of the extrapolation | silent: only the speed of Newton changes |
| `Newton`: Armijo constant $c = 0$; no check that $x$, the step or the line-search residual is finite | each silent |
| `Adaptive`: default `dt_max` of $t_{\mathrm{end}} - t_0$ for $(t_{\mathrm{end}} - t_0) / 10$; reject only a step with $h \le 0$, not one below the resolution of $t$; only one half of the start check (`y0` alone or `rhs t0 y0` alone) | each silent |
| `Halving`: no `max_rejects` give-up; a `nan` estimate accepted | each silent |
| `Stepper.fixed`: drop the `max 1`; steps of `dt` instead of $(t_{\mathrm{end}} - t_0) / n$; the last grid time $t_0 + n h$ instead of $t_{\mathrm{end}}$; $h$ instead of the difference of the end times; only one half of the start check | each silent |

The corpus catches gross breakage (a wrong coefficient, a Jacobian that is plainly wrong, a step that is never allowed to grow back, which exhausts the call budget of the long runs), the pins catch particular mistakes, and it is nearly blind to anything that only changes how fast Newton converges, except through the call count of the adaptive canary: that count moves with the Jacobian (under every swap of two rows, so each swap fails a line of its own) and with Newton's tolerance, but not with the starting guess or the Armijo constant.

### Known gaps

The silent rows are gaps, and starting points for contributions; [numerics/06-the-corpus.md](numerics/06-the-corpus.md) (section 10) suggests a case for most of them.

- **Newton's safeguards.** Damping is pinned (without it the `atan` line prints `Error Diverged`), and so are its depth (`min_damping` $= 1/2$ turns the `deep line search` line into `Error Diverged`), the iteration limit (a limit of 100 turns $x^3$ into `Ok [0.000000000138]`), the refusal of a `nan` trial point (`Error Nan` otherwise) and the failure names at the start of an iteration and for an infinite step. The Armijo constant is not pinned: no line tells $c = 10^{-4}$ from $c = 0$, because the root does not depend on it, so a case has to count evaluations. $x^2 - 5$ from 1 is one: the full step lands on 3, where $\lvert f \rvert$ is 4, the value it has at 1, a tie that only $c = 0$ accepts, and the solve takes 7 evaluations of $f$ with $c = 10^{-4}$ and 6 with $c = 0$. The checks on the step and on the trial residual are redundant in behaviour (the line search refuses the trial point of a non-finite step, and `nan <= bound` is false), and the check that $x$ itself is finite never fires alone in the corpus, so no line pins them. The order of the step test and the line search is pinned, indirectly (table above); a direct pin would be $x^2 - 3$ from 1, which converges today and gives `Error Diverged` when the line search comes first.
- **Defaults and guards of the drivers:** the default `dt_max`, the `max_rejects` give-up, a `nan` estimate accepted, the two halves of each driver's start check (in the two `y0 = nan` lines the right-hand side also returns `nan` for that `y0`, so dropping either half alone changes nothing), the rule that rejects a step below the resolution of $t$ when it would still move $t$ (lines 26 and 28 have steps that snap to 0, and line 27 steps above the resolution; the case would be `dt0 = dt_max = 0.19` from $t = 10^{15}$, which ends in `Error StepRejected 1` today and in `Ok` with `y = 100` without the rule), and three rules of `Stepper.fixed`: at least one step, steps of $(t_{\mathrm{end}} - t_0) / n$ rather than `dt`, and a last grid time that is $t_{\mathrm{end}}$ itself with each step the difference of its end times (in every corpus case `dt` divides the span and the grid reaches `t_end` by itself).
- **van der Pol** checks only `Ok` at $t = 2000$, a rejection and a finite state, not the values against a reference.
- **A transposed Jacobian** is pinned by one line, and the Robertson lines print the budget text within seconds; the soak test, whose adaptive cases are on the budget too, ends in about 14 s with `soak robertson x10: passed 0/10`.
- **The soak test** is a tripwire. **Run time:** the backward Euler canary takes $(t_{\mathrm{end}} - t_0) / \mathtt{dt} = 1 / (2 \times 10^{-6}) = 500{,}000$ fixed steps, once in the corpus and ten times in the soak test, and dominates the run; measure it with `time`.

Closing a gap is a self-contained contribution: write the case that would have noticed the mistake, add it by the steps above, and confirm with the mutation that it now fails. [exercises.md](exercises.md) turns several into starter contributions.
