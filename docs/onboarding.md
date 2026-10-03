# Onboarding: ten working days to a first contribution

This plan takes someone who has never seen this code, has never used OCaml and has not studied numerical analysis to a first small contribution in ten working days. It assumes you program in another language (Python or JavaScript, say), remember first-year calculus (derivatives, Taylor series) and a little linear algebra (matrices, solving `A x = b`). Everything else is taught by the documents it points to.

Days 1 to 8 are foundations: setup, then the code from the smallest module to the adaptive driver, each day with the mathematics it needs, then the tests. Day 9 is for changing the code safely in a copy and day 10 for preparing a real contribution. If you fall behind, shorten the reading before you skip the experiments; [README.md](README.md) has a one-hour path if you only need the idea.

## How the plan works

- Each day has a **goal**, things to **read**, things to **do** (exercises named like A1 or C3 live in [exercises.md](exercises.md)), a **self-check** with answers, and a **done when** line. Try the questions before you read the answers.
- Git and dune commands run from the repository root; probe commands run in the probe project (`../vstiff-scratch`) and scratch-copy commands in the scratch copy (`../vstiff-mutation`), as [exercises.md](exercises.md) shows. Never experiment inside the repository: use a **probe** (a throwaway program in the probe project) to measure things and a **scratch copy** to break things. The [glossary](glossary.md) defines both.
- Read a module's `.mli` first (it documents the API), then the `.ml` (implementation notes).
- When a probe does not build, read the first error message only: later ones often follow from it. "Asking the compiler for a type" in [exercises.md](exercises.md) is the quickest way to see what OCaml thinks something is.
- The numerics chapters are read in the order 2, 3, 1, 4, 5, 6 so that each day's code has its mathematics at hand: Newton's method and the Jacobian belong to the lowest layers of the library, and chapter 1 then ties them to ODEs and stiffness. Each chapter's introduction says what it assumes; if the README left you unsure what an ODE is, read the first pages of chapter 1 on day 1. The "Next" links at the foot of the chapters follow their numbering, not this order.
- Never run `dune promote` to make a failing test pass. Keep a notebook of what surprised you and what you predicted wrongly; days 9 and 10 use it. Look up unfamiliar words and symbols in the [glossary](glossary.md).

## The plan at a glance

| Day | Topic | Exercises ([exercises.md](exercises.md)) |
|---|---|---|
| 1 | Setup, a tour, the first run | A1, A2, the README program |
| 2 | OCaml through `fail` and `vec` | asking the compiler for a type |
| 3 | Linear solves and Newton's method | A3, A4 |
| 4 | Floating point and Jacobians | A5, C2 |
| 5 | ODEs, stiffness, backward Euler, the contracts | B3, C3 |
| 6 | BDF2, history and order | C1, the probes of chapter 4 |
| 7 | Step control and the adaptive loop | B1, B2, B4, C4, C5 |
| 8 | The corpus, the effect quarantine, the test policy | D1 |
| 9 | Breaking the code, and a first method | D2 to D5, E1 |
| 10 | A starter contribution | one of F, or E2 |

## Day 1. Setup and a guided tour

- **Goal.** A working toolchain, a green test run, a map of the repository and a program of your own that calls the library.
- **Read.** The [README](../README.md) at the repository root; the index in [README.md](README.md) of this folder; the tooling part of [ocaml.md](ocaml.md); the overview of [architecture.md](architecture.md) (skim the rest).
- **Do.** Install OCaml 5.5.0 and dune as [ocaml.md](ocaml.md) describes. Clone, `git switch night/wip`, run `dune build` and `dune runtest`. Do A1 and A2, set up the probe project and run the README program. Skim every file of `src/` for a minute, one sentence each, then draw the call chain: `Adaptive.integrate` calls `Bdf2.step_with_error`, which calls `Bdf1.step` and (once there is history) a BDF2 step; both call `Stage.solve`, which hands `Newton.solve` a residual and a Jacobian built on `Jac.forward`; every Newton iteration calls that Jacobian and `Linalg.solve`; `Halving` judges the outcome, and `Clock.resolution` tells `Halving`, `Adaptive` and `Check` how short a step `t` can still resolve.
- **Self-check.**
  1. What does a silent `dune runtest` mean? *Every program's output matched its expected file; a failure prints a diff and exits with status 1.*
  2. What compares the test output with the expected files, `dune build` or `dune runtest`? *Only `dune runtest`. `dune build` runs the two programs, when their output is out of date, to record it and compares nothing.*
  3. Which function integrates with automatic step sizes, and which with equal steps? *`Adaptive.integrate` and `Stepper.fixed`.*
  4. What does `open Vstiff` do? *Dune wraps the library, so its modules are `Vstiff.Vec`, `Vstiff.Newton` and so on; `open` allows the short names.*
  5. Why do experiments happen in a probe or a scratch copy? *Nothing is written in the repository, so there is nothing to undo, and a stray file cannot break the build for others.*
- **Done when.** The tests are green, the README program printed `t = 10000`, and you can draw the call chain from memory.

## Day 2. OCaml basics through `fail` and `vec`

- **Goal.** Read the two smallest modules fluently, ask the compiler for a type, and know what an `.mli` is for.
- **Read.** The language section of [ocaml.md](ocaml.md) from "How the code reads" through "Modules, `.mli` files and abstraction", and its "Warnings are errors". [src/fail.mli](../src/fail.mli), [src/fail.ml](../src/fail.ml), [src/vec.mli](../src/vec.mli), [src/vec.ml](../src/vec.ml). Glossary: variant, result, option, record, labelled argument, pure.
- **Do.** From the repository root start `ocaml` (end each phrase with `;;`, leave with `#quit;;`), type `#use "src/vec.ml";;`, then `axpy 2. [|1.; 2.|] [|10.; 20.|];;` and `add [|1.; 2.|] [|10.|];;`. Ask the compiler for the types of `Vec.axpy`, `Newton.solve` and `Jac.forward` (see "Setup" in [exercises.md](exercises.md); it may spell out the library prefix, `Vstiff.Vec.t`: question 4 of day 1). In a probe write `norm_1` (sum of absolute values, with `Array.fold_left`) and `normalized v = Vec.scale (1. /. Vec.norm_inf v) v`, and check that `v` is unchanged. In a scratch copy add a constructor at the end of `Fail.t` in both `fail.ml` and `fail.mli` and run `dune build @check`.
- **Self-check.**
  1. Why does OCaml have `+.` as well as `+`? *No operator overloading and no silent int-to-float conversion: `+` is for integers, `+.` for floats.*
  2. What does `Vec.add [|1.; 2.|] [|10.|]` do, and with the arguments swapped? *The first raises `Invalid_argument`: the operations index by the first argument's length. The second returns `[|11.|]` and ignores the extra entry.*
  3. What does "never mutated" mean concretely? *No `a.(i) <- x`: every operation allocates its result.*
  4. Which constructor of `Fail.t` carries data, and what happens when you add a constructor? *`StepRejected of int`. `Fail.to_string` stops being exhaustive, warning 8, an error under dune's dev profile: the compiler lists every match that needs a case.*
  5. What is an `.mli` for? *It is the module's API and its documentation: only what it lists can be used from outside.*
  6. What does `let* x = r in e` mean (see `Fail.Syntax`)? *`Result.bind r (fun x -> e)`: continue with `x` when `r` is `Ok x`, otherwise return the `Error` unchanged.*
- **Done when.** You can read `fail` and `vec` without help, your helpers print `norm_1 = 7, norm_inf = 4` and `[0.75; -1]` for `[|3.; -4.|]`, and you have seen the compiler print a type.

## Day 3. Linear solves and Newton's method

- **Goal.** Understand `Linalg.solve` and `Newton.solve` and run both by hand.
- **Read.** [numerics/02-newton.md](numerics/02-newton.md), the whole chapter (its last section, on Newton inside the integrator, previews day 5: skim it now and come back to it then); [src/linalg.mli](../src/linalg.mli), [src/linalg.ml](../src/linalg.ml), [src/newton.mli](../src/newton.mli), [src/newton.ml](../src/newton.ml); in [ocaml.md](ocaml.md) the sections on exceptions and mutable cells and on tail recursion. Glossary: partial pivoting, damping, Armijo condition, residual, inf-norm.
- **Do.** A3 and A4. Read `Newton.solve` case by case (every `if` branch and `match` case) and write down which result each produces and why.
- **Self-check.**
  1. Newton on `x^2 - 2` from `x = 1`: the first two iterates? *1.5 and 1.4166666666666667.*
  2. Why does it return `Error Diverged` from `x = 0`? *The derivative is 0, so the 1x1 Jacobian is singular, `Linalg.solve` returns `None`, and `Newton.solve` turns that into `Diverged`.*
  3. Which damping factors does `Newton.solve` try, and what if none passes? *1, 1/2, 1/4, down to 1/1024, accepting the first that makes the residual fall enough (the Armijo test); if none does, `Error Diverged`.*
  4. When does `Linalg.solve` return `None`, and does `Some` guarantee a good answer? *Only when the best pivot is exactly zero. No: a singular matrix can survive rounding, which is why Newton checks that its step is finite and judges it with the line search.*
  5. How are rows exchanged without mutation? *`pivot` picks the row, `swapped` says which original row sits where after the exchange, and the recursion works on a new, smaller array of rows.*
- **Done when.** Your hand iterates agree with the probe and you can explain every case of `Newton.solve`.

## Day 4. Floating point and Jacobians

- **Goal.** Know why a computer cannot differentiate exactly, how `Jac.forward` approximates the Jacobian and where its error comes from.
- **Read.** [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md); [src/jac.mli](../src/jac.mli), [src/jac.ml](../src/jac.ml); the `jacobian` cases in [test/corpus.ml](../test/corpus.ml) and the comment above them. Glossary: eps, ulp, round-off, cancellation, forward difference, Jacobian, NaN.
- **Do.** A5 and C2.
- **Self-check.**
  1. What are the two errors of a forward difference? *Truncation error, which shrinks with the perturbation, and round-off, which grows like `eps |f| / perturbation`.*
  2. What does `J.(i).(j)` mean, and how many `rhs` calls does `Jac.forward` make? *The derivative of output `i` with respect to input `j`; `n + 1`.*
  3. At `y0 = (1, 1, 1)` how far is the canary's `1e4` entry from exact, and why? *About 3.0e-5: round-off on values near `1e4`; half an ulp of `1e4` divided by the perturbation `2e-8` bounds it near 4.5e-5.*
  4. Why is the corpus bound absolute at the origin and relative at `y0`? *At the origin nothing cancels, so entries are essentially exact; at `y0` the rounding error scales with the entry.*
  5. Why does `Jac.forward` divide by `yp.(j) -. y.(j)` and not by the nominal perturbation? *The sum is rounded; the quotient must match the points where `f` was actually evaluated.*
  6. What does `Vec.norm_inf` return when an entry is NaN, and what follows? *NaN. `Newton.solve` checks `Vec.finite` at the start of every iteration and on a converged `x + dx`, and returns `Error Nan`; `Halving.acceptable` rejects a NaN estimate because `nan <= tol` is false.*
- **Done when.** Your C2 probe prints an error near 3.0e-5 and the matrix `[[5, 2], [1, 3]]`, and you can explain the first with half an ulp.

## Day 5. ODEs, stiffness, backward Euler and the contracts

- **Goal.** Know what an ODE is, why stiffness defeats explicit methods, how backward Euler escapes, how one step is built, and what a contract is.
- **Read.** [numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md), the whole chapter; [src/ode.mli](../src/ode.mli), [src/stage.mli](../src/stage.mli), [src/stage.ml](../src/stage.ml), [src/bdf1.mli](../src/bdf1.mli), [src/bdf1.ml](../src/bdf1.ml), [src/stepper.mli](../src/stepper.mli), [src/stepper.ml](../src/stepper.ml), [src/check.ml](../src/check.ml); in [ocaml.md](ocaml.md) the sections on module types and on modular explicits. Glossary: stiffness, stage equation, A-stability, contract, module type, modular explicit, history.
- **Do.** Run the explicit Euler snippet of chapter 1. B3 and C3. Derive the Jacobian of `G(x) = x - psi - gamma f(t, x)` on paper and find the lines of `Stage.solve` that build it.
- **Self-check.**
  1. What is the stage equation of a backward Euler step? *`x = psi + gamma f(t_{n+1}, x)` with `psi = y_n`, `gamma = h`; its residual has Jacobian `I - gamma J`.*
  2. What is the largest stable step of explicit Euler on `y' = -1e4 y`? *`2e-4`: the factor `1 - 1e4 h` must stay at most 1 in magnitude.*
  3. Why does the corpus take 500,000 backward Euler steps although the method is stable for any step? *Accuracy: the error at `t = 1` is about `0.184 h`, and the line demands less than 1e-6.*
  4. What must a module provide to be an `Ode.Method`? *`type history`, `start` and `step`. `Bdf1.history` is abstract outside `Bdf1`, which is how callers are kept from depending on it.*
  5. Why does `Stepper.fixed` take a module? *So one loop drives any method; the module is passed as a modular explicit.*
- **Done when.** You can derive the explicit and backward Euler factors on `y' = -lambda y` and explain `Stage.solve` and `Bdf1.step` line by line.

## Day 6. BDF2, history and order

- **Goal.** Understand how BDF2 is derived, what order two means, and how the code handles the first step and variable steps.
- **Read.** [numerics/04-bdf.md](numerics/04-bdf.md); [src/bdf2.mli](../src/bdf2.mli), [src/bdf2.ml](../src/bdf2.ml). Glossary: BDF, multistep method, order, variable step, zero-stability, embedded method, history, ω.
- **Do.** The probes of chapter 4 (they evaluate `Bdf2.coeffs`; never copy its numbers into a comment or a document). C1. Trace the first two calls of `Bdf2.step` by hand: which stage solve each makes and what the history holds afterwards. Derive `a1`, `a0` and `beta` on paper from the conditions in chapter 4 and compare with `Vstiff.Bdf2.coeffs`.
- **Self-check.**
  1. Which polynomial does BDF2 differentiate, and what does it enforce? *The quadratic through `(t_{n-1}, y_{n-1})`, `(t_n, y_n)` and the unknown `(t_{n+1}, y_{n+1})`; its slope at `t_{n+1}` must equal `f(t_{n+1}, y_{n+1})`.*
  2. Why is ω = h / h_prev the only input of `Bdf2.coeffs`? *The exactness conditions involve the time differences only as ratios.*
  3. What does the history hold, and why is the first step backward Euler? *`Start` or `After { h_prev; y_prev }`; with `Start` there is no previous point. Chapter 4 explains why order two survives.*
  4. What does the corpus ratio 3.99 show? *Halving `dt` divided the error by about 4: second order.*
  5. Where are the values of `a1`, `a0` and `beta`? *In `Bdf2.coeffs`, the single source of truth: evaluate it or derive them.*
- **Done when.** Your C1 orders are near 1 and 2 and your derivation agrees with `Bdf2.coeffs`.

## Day 7. Step control and the adaptive loop

- **Goal.** Narrate `Adaptive.integrate` from start to finish: choose and snap the step, error estimate, accept, reject, grow, give up.
- **Read.** [numerics/05-step-control.md](numerics/05-step-control.md); [src/halving.mli](../src/halving.mli), [src/halving.ml](../src/halving.ml), [src/clock.mli](../src/clock.mli), [src/clock.ml](../src/clock.ml), [src/adaptive.mli](../src/adaptive.mli), [src/adaptive.ml](../src/adaptive.ml); the data-flow and error-flow parts of [architecture.md](architecture.md). Glossary: local truncation error, scaled error, rejection, step floor, snapping, `tol`, controller.
- **Do.** B4, B1, B2, C4 and C5. Close `adaptive.ml` and write `go` from memory in plain words, then fix your version against the file.
- **Self-check.**
  1. What does `Halving` compare with `tol`? *`max_i |err_i| / (1 + |y_i|)`.*
  2. What does the error estimate measure on the first step, and later? *First step: half the gap between backward Euler and explicit Euler. Later: the gap between BDF2 and backward Euler.*
  3. When does the step double, and what does a rejection do? *After three accepts in a row, up to `dt_max`. A rejection, for any reason, halves the step that failed and restarts the count.*
  4. Which failures can `Adaptive.integrate` return with `Halving`? *`Nan` before the first step, from `Adaptive`; `StepRejected n` from `Halving`. A failed Newton solve is just a rejection.*
  5. Why does the NaN-wall corpus line say `StepRejected 46` when `max_rejects` is 50? *The floor `16 eps |t|` ends the run before the count does (C4).*
  6. What are the three reasons for a rejection? *`Too_large`, the estimate exceeded `tol`; `Solver e`, the method could not take the step (Newton failed); and `Too_small`, the driver found the step too short to move `t` and did not call the method. The controller sees which, and `Halving` treats them alike.*
  7. What does the driver do with a step shorter than half an ulp of `t`? *It snaps the step to `h = (t + dt) - t = 0` and rejects it as `Too_small`; `Halving` halves 0, which is below the floor, so the run ends with `StepRejected 1`.*
- **Done when.** Your B4 table matches, you have measured counts for three tolerances, and you can retell `go` without looking.

## Day 8. The corpus, the effect quarantine and the test policy

- **Goal.** Know what every test line proves and cannot catch, which modules may have effects, and the rules for changing expectations.
- **Read.** [numerics/06-the-corpus.md](numerics/06-the-corpus.md); [testing.md](testing.md); the contracts and effects parts of [architecture.md](architecture.md); the conventions in [CONTRIBUTING.md](../CONTRIBUTING.md); [test/problems.ml](../test/problems.ml), [test/refs.ml](../test/refs.ml), [test/guard.ml](../test/guard.ml), [test/report.ml](../test/report.ml), [test/corpus.ml](../test/corpus.ml), [test/soak.ml](../test/soak.ml) and both expected files. Glossary: corpus, canary, expect test, soak test, effect quarantine.
- **Do.** Make a table with a row per corpus and soak line: problem, library code exercised, one defect it would catch and one it would not. Do D1 in a scratch copy. Read two commits that added corpus cases: `git log --oneline -- test/corpus.expected` lists them and `git show <hash> -- test/corpus.expected` shows the line one added. The recent ones add a fix together with its pins and use the current API; the older ones use an older version of the API and of `corpus.ml`, so read those for the purpose of the case and its expected line, not for the code.
- **Self-check.**
  1. Why are reference values computed outside vstiff and never edited? *A reference computed by vstiff agrees with vstiff whatever it does.*
  2. What is the soak test for? *A tripwire: the library has no hidden state, so `identical` turns false only if hidden state, randomness, parallelism or a NaN appears.*
  3. When is `dune promote` legitimate? *To record the line of a newly added case, or a deliberate, reviewed format change; never to turn red green.*
  4. Which modules may raise, mutate or print? *Library: `Check` raises, `Instrument` has mutable state. Tests: `Guard` raises and catches, `Report` prints. The rest is pure.*
  5. Which lines notice swapped Jacobian rows (D1)? *Five change. The `bdf1 canary` line stays `true` and only its digits move, because Newton's root does not depend on the Jacobian, only its speed does.*
- **Done when.** Your table is complete, you saw D1's five changed lines, and you can state the policy for expected files.

## Day 9. Breaking the code in a scratch copy, and a first method

- **Goal.** Change the code safely, predict what the tests say, and check.
- **Read.** The mutation part of [testing.md](testing.md); sections D and E of [exercises.md](exercises.md); the code and comment style in [CONTRIBUTING.md](../CONTRIBUTING.md).
- **Do.** D2, D3, D4 and D5 (under a time limit), then E1 up to its order check. For each experiment write down what you changed, what you predicted, what the tests said and which line would have caught it.
- **Self-check.**
  1. After breaking pivoting, which lines changed? *Two `linalg` lines; nothing else depends on the pivot choice.*
  2. How do you check that a new test can fail? *Break the thing it protects in a scratch copy and see the line change.*
  3. Why a scratch copy rather than `git stash`? *There is nothing to undo, your tree stays untouched, and builds cannot interfere.*
  4. After D4, which Robertson line failed and which did not? *The line with bound `1e-5` failed; the other one, with bound `1e-3`, passed. That is why the tighter line exists.*
- **Done when.** You have written up four experiments, and your `Trapezoid` method builds and passes its order check in a probe.

## Day 10. Picking and preparing a starter contribution

- **Goal.** Choose a contribution, plan it and prepare it so that a reviewer can accept it.
- **Read.** All of [CONTRIBUTING.md](../CONTRIBUTING.md), checklist included; section F of [exercises.md](exercises.md); your notebook.
- **Do.** Pick one item you understand: F1 needs no discussion; F2 adds a corpus case and F9 changes the corpus, so they are discussed first; F3 to F8 are design-first. Write a plan of about five sentences: goal, files, the check that shows it works, documents to update, what you will not do. Branch from `night/wip` (`git switch -c docs/short-slug night/wip`), commit in small Conventional Commits steps and run the checklist item by item. Hand the branch over by pushing it (`git push -u origin <branch>`, to your fork if you cannot push to the project) and asking your reviewer to take it, as [CONTRIBUTING.md](../CONTRIBUTING.md) says; do not merge it.
- **Self-check** (the answers are the items of the checklist). Does the change alter behaviour, and which test covers it or which gap remains? Did an expected file change, and only by the lines of a case discussed first? Are effects still quarantined? Do the documents you touched still tell the truth, with every snippet compiled in a probe? Is every commit one logical change?
- **Done when.** A branch exists whose diff you would defend line by line and every checklist item holds.

## After ten days

You can build the project and run its tests, read a test diff and explain why `dune promote` does not fix a failing test; read any module, say what it does and find the chapter behind it; trace an adaptive step from `Adaptive.integrate` to `Linalg.solve`; check a claim with a probe; judge a test by breaking the code; and make a small change that follows the conventions. Changes to stability-sensitive logic (the step-ratio bound, Newton's safeguards, the error estimate), and new methods or controllers, start as a design discussion with a reviewer: the capstones E1 and E2 and the design-first items of section F are where to practise them.
