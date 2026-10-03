# Contributing to vstiff

Work happens on the `night/wip` branch and nothing is merged for now (see [README.md](README.md)). If OCaml or numerical methods are new to you, start with [docs/onboarding.md](docs/onboarding.md); the [glossary](docs/glossary.md) explains the terms below and [docs/README.md](docs/README.md) indexes all documentation.

## Setting up

You need OCaml 5.5.0 and dune 3.x on macOS, Linux or WSL (Windows Subsystem for Linux); [docs/ocaml.md](docs/ocaml.md) shows how to install them and how to get type-on-hover in your editor. The project depends on nothing else. `<repository-url>` is the address of the project's repository, which you were given with your access.

```sh
git clone <repository-url> vstiff && cd vstiff
git switch night/wip
dune build
dune runtest
```

`dune build` and `dune runtest` print nothing and exit with status 0 when everything is in order.

## The edit, build, test loop

| Command | What it does |
|---|---|
| `dune build @check` | Type-checks everything without linking or running anything: the fastest way to see compile errors. |
| `dune build` | Compiles both libraries and the tests and runs the two test programs, when their output is out of date, to record it. It compares nothing. |
| `dune runtest` | Compares the recorded output of `corpus` and `soak` with [test/corpus.expected](test/corpus.expected) and [test/soak.expected](test/soak.expected). Silent on success; on failure it prints a diff (`-` is the expected line, `+` what was printed) and exits with status 1. |
| `dune exec ./test/corpus.exe` | Runs the corpus and prints its output unfiltered; `./test/soak.exe` does the same for the soak test. |
| `dune promote` | Overwrites expected files with the last output. Legitimate only as described under "Expect files are the contract" below. |

A round: edit; run `dune build @check` until it compiles (dune's dev profile turns many warnings into errors: an unused value, a `match` that misses a case, a documentation comment in the wrong place; fix them, never loosen the flags); run `dune runtest` and read the diff if there is one; commit when it is green. The backward Euler canary takes 500,000 steps and dominates the test time, so while you iterate on something that cannot change numerical behaviour, such as a comment, `dune build @check` is the cheap check. Run the full `dune runtest` before you commit. If a change might make a run never finish (an edit to `Newton`, `Halving` or `Adaptive`, say), run `dune build ./test/corpus.exe ./test/soak.exe` first (it compiles without running anything), then put the time limit on the executable, `perl -e 'alarm 120; exec @ARGV' ./_build/default/test/corpus.exe`, not on `dune build` or `dune runtest`: killing dune can leave the test programs running. The corpus runs van der Pol, Robertson, the adaptive canary and its termination cases on a call budget, so a run that crawls there prints `no answer within 5e6 rhs calls` on its line within seconds; the soak test's van der Pol and Robertson cases are on the same budget, and a round that exhausts it counts as not passed.

To check a claim about the library, or to try a change you are not ready to commit, leave the repository alone: [docs/exercises.md](docs/exercises.md) describes the probe project (a throwaway program that links to `src/`) and the scratch copy (a copy of the repository you can break).

## Project conventions

Each rule comes with its reason, because a rule you understand is one you can apply to a case it does not mention.

**Effects are quarantined.** In the library only `Check` raises and only `Instrument` has mutable state; in the tests only `Guard` raises or catches and only `Report` prints. Methods, controllers, drivers and the whole kernel (`Newton`, `Jac`, `Linalg`, `Vec`, `Fail`) are pure. A new effect goes into one of the four effect modules (the effects map below) or into a new dedicated module, never into a pure one. *Why:* a pure function's result depends only on its arguments, so the adaptive driver can retry a rejected step from unchanged data, the soak test can demand identical reruns, and effect types stay easy to annotate. (Standard-library functions can still raise on misuse, for example an out-of-bounds index when two vectors differ in length; that is a bug, not behaviour to rely on.)

**Contracts are module types in `Ode`.** Methods and controllers are modules implementing `Ode.Method`, `Ode.Embedded` and `Ode.Controller`, passed to the drivers as modular explicits: `Stepper.fixed (module Bdf1) ~dt problem`. Prefer modular explicits to functors and to first-class modules packed into values ([docs/ocaml.md](docs/ocaml.md) explains all three), and share types between signatures (`with type`) only when necessary (`halving.mli` shares `stats` so callers can read its fields). *Why:* a driver depends on what a method provides, not on which method it is, and the result type of `Adaptive.integrate` can depend on the controller passed in.

**Numerical failures are values; invalid arguments raise.** A solve or a step that cannot succeed returns a `Fail.t` in a `result`; arguments that make no sense (`dt <= 0`, `t_end < t0`, `tol <= 0`, an `rhs t0 y0` of the wrong length) raise `Invalid_argument` from `Check`. *Why:* a failed Newton solve or a rejected step is an expected event that the controller handles by halving the step, and a `result` forces every caller to decide what to do; an invalid argument is a bug in the caller and should stop it loudly.

**Two libraries, one public API.** The kernel in `src/numerics/` (`Fail`, `Vec`, `Linalg`, `Newton`, `Jac`) knows nothing about ODEs: nothing in it names `Ode`, and dune would reject a dependency of `numerics` on `vstiff`. The solver in `src/` offers only what [src/vstiff.ml](src/vstiff.ml) and [src/vstiff.mli](src/vstiff.mli) list: a module is public if and only if it is listed there, and `Stage` and `Check` are not (they are also `private_modules`). Public signatures say `float array`, never `Vec.t`, and a solver module names what it uses from the kernel at its top (`module Vec = Numerics.Vec`). *Why:* the public surface stays small enough to document and keep stable, a user who writes a method or a controller needs the solver library alone, and the kernel can be read and tested with no ODE in sight.

**Naming.** Module types are CamelCase, values snake_case. `dt` is a requested step size (`~dt`, `dt0`, `dt_max`, `Controller.proposal`); `h` is a step taken or attempted (`h_prev`, the `h` of `Method.step`). Labels appear only where two arguments of one type could be swapped (`~at ~h`), or to name a bare literal at a call site (`~tol:1e-6`). *Why:* a reader tells a request from an outcome at a glance, and labels stay informative instead of decorating every call.

**Jacobians are always forward differences inside the solver.** `Jac.forward` approximates the Jacobian from function values, and a corpus problem never supplies an analytic one (`corpus.ml` writes analytic Jacobians only as the reference for `Jac.forward` and to drive `Newton` alone). *Why:* the solver has to work for any black-box right-hand side, and the corpus has to keep exercising `Jac`.

**Expect files are the contract.** Never weaken an expectation or edit one to get green. A deliberate change of algorithm (a new controller, estimate or method) re-pins the lines it moves, such as step counts and digits, in the same commit, and the commit message gives the old and the new values; the bounds that define correctness (`< 1e-6`, `< 1e-3`, `< 1e-8`, the order window `[3.5, 4.5]`) are never loosened. `dune promote` only records a newly added case, such a deliberate re-pin, or a reviewed change of format. Checks are about correctness, not bits. A refactor may change the last bits of a result (fused multiply-add, another order of summation, another platform), and the expected files print only the digits their checks need, so they do not move. A change that moves a line either changed behaviour, or exposed a line printed more precisely than its check needs; the commit message says which. *Why:* the expected files are the project's memory of what correct output looks like; a wrong answer that was promoted stays silent forever.

**References come from outside.** The values in [test/refs.ml](test/refs.ml) are computed outside vstiff and never edited; record how a new one was computed. *Why:* a reference computed by vstiff agrees with vstiff whatever it does.

**No coefficient tables.** No comment or document contains BDF coefficient values, tables or the formulas for them. `Bdf2.coeffs` is the single source of truth: explain how the coefficients are derived and point readers to `Vstiff.Bdf2.coeffs`. *Why:* a copy can be mistyped or go stale, while the code can always be asked.

**Standard library only.** `src/numerics/dune` lists no libraries, `src/dune` lists only `numerics` and `test/dune` lists `vstiff` and `numerics`. If you think a dependency is justified, state the reason in the change description. *Why:* the project builds with a bare OCaml and dune installation, and each dependency is one more thing that can break.

**Small commits, work on `night/wip`, nothing merged.** One logical change per commit, in Conventional Commits style (below). *Why:* a reviewer checks one idea at a time, and a commit that turns out wrong can be reverted without losing the rest.

## The effects map

| Module | Effect |
|---|---|
| `Check` (library, internal) | raises `Invalid_argument` for bad driver arguments |
| `Instrument` (library) | a mutable call counter wrapped around an `rhs` |
| `Guard` (tests) | raises `Exhausted` past its call budget; catches `Exhausted` and `Invalid_argument` and turns them into text (`run`), or `Exhausted` into `None` (`bounded`) |
| `Report` (tests) | prints one line per case |
| everything else | pure |

## Code and comment style

- **Match the surrounding code.** Read the neighbouring files first. There is no formatter configuration, so keep indentation, naming and line width (the code stays near 120 columns) consistent by hand.
- **Comments are concise and say why.** `.mli` files are the API documentation: a module header of 1 to 4 lines, then 1 to 3 lines for each type, constructor, field and `val` (one whose name and type say it all, like `Vec.add`, may leave it to the module header): what it means, its inputs and outputs, how it fails. `.ml` files carry implementation notes only: why this approach, a numerical reason, an invariant, or a non-obvious OCaml construct at its first appearance (one line, pointing to [docs/ocaml.md](docs/ocaml.md)). Each comment is 1 to 3 lines and a module header at most 5, naming the idea and the chapter that explains it. Do not repeat the `.mli`, narrate the code or write tutorials: depth belongs in `docs/`, and a module's comment lines (`.ml` and `.mli` together) should not outnumber its code lines. In tests, a short header per file and one short comment per group of cases saying what it proves and why its threshold.
- **Math in code comments is ASCII-friendly:** `y_{n+1}`, `|x|_inf`, `omega`, `d f_i / d y_j`. In Markdown use plain text or Unicode math in prose or fenced blocks, never LaTeX.
- **Documentation comments** `(** ... *)` go immediately before the item they document, or right after a record field or a constructor. Inside a function body use `(* ... *)`: a documentation comment there is warning 50, an error here.
- **Comment pitfalls** that break the build: comments nest, so an opening parenthesis directly followed by a star inside comment text opens a nested comment (write `( *. )` with spaces); string literals are lexed inside comments, so avoid double quotes, and `{|` likewise; a star directly followed by a closing parenthesis ends the comment.
- **Documents:** every sentence true; a number comes from an expected file, `refs.ml`, a derivation shown in the text, or a probe the reader can run, never from memory. Use the symbols of the [glossary](docs/glossary.md). Link to files with relative paths, never to `#anchors`. Every OCaml snippet meant to compile must compile in a probe (one that names `Vec`, `Newton`, `Jac` or `Linalg` needs `numerics` in its `libraries` and `open Numerics`). No coefficient values.

## Tests

Behaviour changes come with a test; if a gap remains, say so in the change description. [docs/testing.md](docs/testing.md) explains how the tests are built and run, and [docs/numerics/06-the-corpus.md](docs/numerics/06-the-corpus.md) what each corpus problem proves.

**Expect files.** A test program prints and dune compares the output with the `.expected` file next to it. If an existing line changes when you did not mean to change behaviour, you found a regression: fix the code, not the expectation. The rules for `dune promote` are under "Project conventions".

**Adding a corpus case.** New cases are added deliberately, one at a time.

1. State the purpose in one sentence: which defect would this case catch that no existing line catches? If you cannot, do not add it. Raise the case for discussion before you write it.
2. If it needs a new problem, add a module to [test/problems.ml](test/problems.ml) with `rhs`, `y0`, a `problem` record and, when the solution is known, `exact`. No analytic Jacobian.
3. If it needs a reference value, compute it outside vstiff and add it to [test/refs.ml](test/refs.ml) with how it was computed.
4. Add `(name, fun () -> text)` in a new group at the end of the `corpus` list in [test/corpus.ml](test/corpus.ml) ([docs/testing.md](docs/testing.md) says why at the end: a case placed elsewhere shifts the line numbers that [docs/numerics/06-the-corpus.md](docs/numerics/06-the-corpus.md) cites, which would then need renumbering). Print a stable line: a few digits or a boolean for a bound. A case that might not terminate runs its `rhs` behind `Guard.budget` and the case behind `Guard.run`.
5. Run `dune runtest`: the diff must contain your line and nothing else. If another line changed, fix that first.
6. Check the new line by reasoning, not because the program printed it. Then run `dune promote` and commit the expected file with the case.
7. Check that the case can fail: in a scratch copy, break what it protects and see the line change.

**Adding a soak case.** Write a module with the `Case` signature in [test/soak.ml](test/soak.ml) and add `soak (module YourCase)` to the list at the end of the file; its line goes into [test/soak.expected](test/soak.expected) by the same steps. A case that runs the adaptive driver wraps its problem with `on_budget`, as `VanDerPol` and `Robertson` do: a bug that makes it crawl then ends the round (`Guard.bounded` turns the exhausted budget into `None`, a round that did not pass) instead of stalling the soak.

**Adding a method.** Write `src/your_method.ml` with `type history`, `start` and `step` (and `step_with_error` for an `Ode.Embedded`), and an `.mli` that says `include Ode.Method` (or `Ode.Embedded`), as `bdf1.mli` and `bdf2.mli` do: the history stays abstract. Dune picks up new modules in `src/` by itself. If the method fits the stage equation `x = psi + gamma f(t_{n+1}, x)`, call `Stage.solve`. Keep it pure. To make it public, so that the tests and other programs can use it, add `module Your_method = Your_method` to [src/vstiff.ml](src/vstiff.ml) and the same line with a one-line doc comment to [src/vstiff.mli](src/vstiff.mli); until then it is internal, like `Stage`. Then drive it with `Stepper.fixed (module Your_method) ~dt problem` and add corpus cases as above (the `bdf2 logistic` line, an order check, is the model). A method written outside the library, in a probe or another project, implements `Ode.Method` in the same way but cannot see `Stage`: it builds its equation with `Newton` and `Jac` from the kernel ([docs/exercises.md](docs/exercises.md), E1). See [docs/architecture.md](docs/architecture.md).

**Adding a controller.** Write `src/your_controller.ml` with `t`, `stats`, `init`, `proposal`, `acceptable`, `accepted`, `rejected` and `stats`, and an `.mli` like [src/halving.mli](src/halving.mli): declare `stats` concretely, then `include Ode.Controller with type stats := stats`. List it in `src/vstiff.ml` and `src/vstiff.mli` to make it public, and run it with `Adaptive.integrate (module Bdf2) (module Your_controller)`; a controller outside the library needs the solver library alone. The contract passes the error estimate only to `acceptable`; a policy that needs it elsewhere needs a contract change, which is discussed first. `rejected` must end every run of rejections with an `Error`, as `Halving` does with `max_rejects`, the step floor and a halved step of 0: the termination of `Adaptive.integrate` rests on it ([docs/numerics/05-step-control.md](docs/numerics/05-step-control.md), section 9). See [docs/architecture.md](docs/architecture.md) and [docs/testing.md](docs/testing.md).

## Commits and branches

Commits follow Conventional Commits: `type(scope): imperative summary`, where the type is `feat`, `fix`, `perf`, `refactor`, `test`, `docs` or `chore` and the scope names the module or area, such as `adaptive`, `newton`, `corpus` or `docs`. No trailing period; use the body for why. For example `fix(halving): stop when a halved step falls below the floor` or `docs: explain the stage equation in the glossary`. One logical change per commit: do not mix a refactor with a behaviour change, and keep a fix and the test that pins it together.

Base your work on the latest `night/wip` and give each branch (one per change) a name made of the commit type and a short subject, for example `docs/glossary-typos`. Hand the branch over by pushing it to the remote you cloned from (`git push -u origin <branch>`), or to your fork if you cannot push there, and asking your reviewer to take it; they will say whether they want a pull request against `night/wip`. Do not merge it, and do not merge anything into `main`.

## Review checklist

Run through this before you ask for review.

- [ ] `dune build @check`, `dune build` and `dune runtest` pass: no warnings, no diff.
- [ ] No expected file changed except by the lines of a new case that was discussed first, and no existing expectation was weakened: no looser bound, no rewritten line, no removed case.
- [ ] A refactor leaves both expected files unchanged, or its commit message says why a line moved.
- [ ] Effects are still quarantined; numerical failures are `Fail.t` values and only `Check` raises.
- [ ] No corpus problem supplies an analytic Jacobian, and no new dependency was added.
- [ ] The public API is deliberate: a module users should see is listed in `src/vstiff.ml` and `src/vstiff.mli`, an internal one is not, public signatures say `float array`, and nothing in `src/numerics/` names `Ode`.
- [ ] A new reference value was computed outside vstiff, with how written next to it; no existing one was edited.
- [ ] Every behaviour change has a test, or the description names the gap that remains.
- [ ] `.mli` comments say what, inputs, outputs and failure in 1 to 3 lines; `.ml` comments say why and do not narrate; no coefficient values, tables or formulas.
- [ ] Documents affected by the change are updated, relative links work and none uses an anchor, every number comes from an expected file, `refs.ml`, a derivation or a probe, and every OCaml snippet was compiled in a probe.
- [ ] One logical change per commit in Conventional Commits style, no unrelated formatting churn, the branch is based on the current `night/wip`, and nothing was merged.
