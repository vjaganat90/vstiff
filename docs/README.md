# vstiff documentation

These pages explain vstiff to someone who has never seen the code, has never used OCaml and has not studied numerical analysis.

Throughout, $t$ is time, $y$ the state vector, $f$ the right-hand side of $y' = f(t, y)$ (called `rhs` in code), $h$ a step taken or attempted and `dt` a step size someone asked for. The [glossary](glossary.md) lists every term and symbol.

## Before you start

What you need before reading these pages or the code is stated here, once, for every document.

**What to know.** You program in some language (Python or JavaScript, say), remember first-year calculus (derivatives and Taylor series) and a little linear algebra (matrices, and solving $A x = b$). The pages teach everything else from the start: ODEs and stiffness, floating point, Newton's method, finite differences, BDF methods, step-size control, OCaml and dune, and how the tests work. Eigenvalues and Jacobians are explained where they first appear.

**What to install.** Only the snippets, the exercises and the code need a toolchain; the numerics chapters can be read without one. The toolchain is OCaml 5.5 and dune 3, usually installed through opam; the libraries and the tests use only the OCaml standard library, so no other OCaml package is needed. [ocaml.md](ocaml.md) shows how to install opam, OCaml and dune. The shell recipes in these pages also use `git`, `ln`, `tar` and `perl`, and run on macOS, Linux or WSL (Windows Subsystem for Linux). The scipy scripts of the bench, in [bench/compare/](../bench/compare/), are outside the build and need Python 3.13 with the packages that [bench/compare/requirements.txt](../bench/compare/requirements.txt) pins; only a reader who runs the bench's Python side (its references, scipy rows or transcription check) installs them.

**How to get the code.**

```sh
git clone https://github.com/vjaganat90/vstiff
cd vstiff
dune build
dune runtest
```

`dune build` compiles the libraries, the tests and the bench and also runs the test programs to record their output, so it takes a little while. `dune runtest` compares that output with [test/corpus.expected](../test/corpus.expected), [test/soak.expected](../test/soak.expected) and [bench/test/check.expected](../bench/test/check.expected): silence means it matches, otherwise dune prints a diff ("Reading a failure" in [testing.md](testing.md) explains it).

**Which document needs which.**

- **The numerics chapters** build on each other, so read them in order, 1 to 6. Chapter 3 needs chapter 2, chapter 4 needs chapter 1 and a first look at chapter 2, chapter 5 needs chapters 1 and 4, and chapter 6 needs chapters 1 to 5 and [testing.md](testing.md), which says how the tests run. Chapters 1 and 2 do not need each other, except that the last section of chapter 2, on Newton inside the integrator, reads best after chapter 1. The programs in all six chapters are optional; the chapters can be read without running OCaml.
- **The snippets and the exercises** run in one probe project outside the repository, which you set up once, as "Running the snippets" in [chapter 1](numerics/01-odes-and-stiffness.md) or "Setup" in [exercises.md](exercises.md) shows. The two recipes differ only in `p/dune`: chapter 1 lists the `numerics` library too, which a probe that calls the kernel needs.
- **The code** is written in OCaml, so [ocaml.md](ocaml.md) comes before it: its first part gets the toolchain running, and its tour of the language is read alongside the first modules.
- **A first change**, to the code or to a document, comes after [AGENTS.md](../AGENTS.md), which holds the rules for every change.

## Index

| Document | What it is |
|---|---|
| [README.md](../README.md) (repository root) | The front page: what vstiff is, a first program, quick start, status and known limits, layout. |
| [AGENTS.md](../AGENTS.md) (repository root) | The rules for every change: the hard rules H1 to H17, the design defaults, the build, test and bench commands, how to add a corpus case, a soak case, a method or a controller, the effects map, commits, pull requests and the checklist for a finished change. [CONTRIBUTING.md](../CONTRIBUTING.md) points to it. |
| [docs/README.md](README.md) | This page: what to know before you start, the index, reading orders and a one-hour path. |
| [onboarding.md](onboarding.md) | A ten-working-day plan from zero to a first contribution, with readings, exercises, self-checks and a done-when line for each day. |
| [glossary.md](glossary.md) | Terms and symbols, each with a short definition and where it is explained. |
| [exercises.md](exercises.md) | Graded exercises (warm-ups, reading, probes, breaking the code in a scratch copy), a capstone against the contracts, and starter contributions. |
| [ocaml.md](ocaml.md) | OCaml and tooling primer: opam, dune, the commands you will type, and the language features the code uses, libraries and main modules included. |
| [architecture.md](architecture.md) | The two libraries and the public API, the contracts, the modules, how data flows between them, and where effects live. |
| [testing.md](testing.md) | How the tests are built and run, how to read a failure, the policy for changing expectations, adding a case, probes. |
| [bench/README.md](../bench/README.md) | The work-precision bench: how to run it, its golden table and gate, the references and scipy rows, the transcription check. |
| [numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md) | ODEs, explicit and backward Euler, order, stability and stiffness. |
| [numerics/02-newton.md](numerics/02-newton.md) | Newton's method, damping and the linear solves inside it. |
| [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md) | Floating point, finite-difference Jacobians, round-off and NaN. |
| [numerics/04-bdf.md](numerics/04-bdf.md) | BDF1, BDF2, variable steps, the stage equation and stability. |
| [numerics/05-step-control.md](numerics/05-step-control.md) | Error estimates, step-size control, the clock that steps are snapped to (and the grid of the fixed-step driver), and why every run ends. |
| [numerics/06-the-corpus.md](numerics/06-the-corpus.md) | The test problems and what each one proves. |

## Reading orders

**Following the onboarding plan.** [onboarding.md](onboarding.md) interleaves everything below over ten working days. It is the recommended route if you will contribute.

**The numerics track.** Read the six chapters in `numerics/` in order, 1 to 6; "Before you start" says which needs which. Each ends with questions to check yourself. The OCaml snippets are optional but worth running. ([onboarding.md](onboarding.md) visits the chapters in the order 2, 3, 1, 4, 5, 6, which those dependencies allow, to match the code it asks you to read each day.)

**The code track.** Learn OCaml from [ocaml.md](ocaml.md) while reading the code from the smallest pieces towards the drivers: first the kernel, the library `numerics`, then the solver, the library `vstiff`. For each module read the `.mli` first, which documents the API, then the `.ml`, which has implementation notes:

1. [src/numerics/fail.mli](../src/numerics/fail.mli): the named failures and the `let*` operators.
2. [src/numerics/vec.mli](../src/numerics/vec.mli): vectors that are never mutated.
3. [src/numerics/linalg.mli](../src/numerics/linalg.mli): the linear solve.
4. [src/numerics/newton.mli](../src/numerics/newton.mli): damped Newton's method.
5. [src/numerics/jac.mli](../src/numerics/jac.mli): forward-difference Jacobians.
6. [src/vstiff.mli](../src/vstiff.mli) and [src/fail.mli](../src/fail.mli): the public API of the solver, and how it re-exports the kernel's failures.
7. [src/ode.mli](../src/ode.mli): the contracts every method and controller implements.
8. [src/stage.mli](../src/stage.mli): the implicit equation every BDF step solves.
9. [src/check.mli](../src/check.mli) and [src/instrument.mli](../src/instrument.mli): the only effects in the library.
10. [src/bdf1.mli](../src/bdf1.mli), [src/bdf2.mli](../src/bdf2.mli): the methods.
11. [src/clock.mli](../src/clock.mli) and [src/halving.mli](../src/halving.mli): the resolution of time, and the step-size controller.
12. [src/stepper.mli](../src/stepper.mli), [src/adaptive.mli](../src/adaptive.mli): the two drivers.

Then the tests: [test/problems/problems.ml](../test/problems/problems.ml), [test/refs.ml](../test/refs.ml), [test/guard.ml](../test/guard.ml), [test/report.ml](../test/report.ml), [test/corpus.ml](../test/corpus.ml) with [test/corpus.expected](../test/corpus.expected), and [test/soak.ml](../test/soak.ml) with [test/soak.expected](../test/soak.expected).

**Where each module is explained.**

| Source | Read about it in |
|---|---|
| `vstiff`, `numerics` (the libraries) | [architecture.md](architecture.md) (the two libraries, the public API), [ocaml.md](ocaml.md) (libraries, main modules, aliases, re-exports) |
| `fail`, `vec` | [ocaml.md](ocaml.md) (variants, arrays, `let*`), [architecture.md](architecture.md) |
| `linalg`, `newton` | [numerics/02-newton.md](numerics/02-newton.md) |
| `jac` | [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md) |
| `ode`, `stepper` | [architecture.md](architecture.md) (contracts), [ocaml.md](ocaml.md) (module types, modular explicits), [numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md) |
| `stage`, `bdf1` | [numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md), [numerics/04-bdf.md](numerics/04-bdf.md) |
| `bdf2` | [numerics/04-bdf.md](numerics/04-bdf.md) |
| `clock`, `halving`, `adaptive` | [numerics/05-step-control.md](numerics/05-step-control.md) |
| `check`, `instrument` | [architecture.md](architecture.md) (effects) |
| `test/` | [testing.md](testing.md), [numerics/06-the-corpus.md](numerics/06-the-corpus.md) |

## If you only have an hour

If the toolchain is installed ("Before you start" above says what you need), start `dune build` in a terminal first: it also runs the test programs, so it takes a little while.

1. About 5 minutes: the [README](../README.md), for what the library does and how the repository is laid out.
2. About 15 minutes: [architecture.md](architecture.md), for the contracts, the modules and how data flows between them. Skim the rest.
3. About 15 minutes: the parts of [numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md) on stability, stiffness and backward Euler. They explain why the library exists.
4. About 15 minutes: [src/ode.mli](../src/ode.mli), then [src/stage.ml](../src/stage.ml), [src/bdf1.ml](../src/bdf1.ml) and [src/adaptive.ml](../src/adaptive.ml). Between them they show a contract, an implicit step and the loop that drives it.
5. About 5 minutes: [test/corpus.expected](../test/corpus.expected), then `dune runtest`, which prints nothing when everything passes.

Keep the [glossary](glossary.md) open for the terms you meet. To carry on, start [onboarding.md](onboarding.md) at day 1.

## Finding what you need

| You want to | Go to |
|---|---|
| Know what to learn and install first | "Before you start", above |
| Know what a word or symbol means | [glossary.md](glossary.md) |
| Know what is public API and what is internal | [architecture.md](architecture.md), [src/vstiff.mli](../src/vstiff.mli) |
| Run a piece of the library and see what it does | the probe project in [exercises.md](exercises.md) |
| Understand a failing test | [testing.md](testing.md) |
| Know the rules for a change | [AGENTS.md](../AGENTS.md) |
| Add a test | [testing.md](testing.md) and the Tests section of [AGENTS.md](../AGENTS.md) |
| Add a method or a controller | the Tests section of [AGENTS.md](../AGENTS.md) (Adding a method or a controller) and [architecture.md](architecture.md) |
| Find something small to work on | the starter contributions in [exercises.md](exercises.md) |
| Fix or add documentation | the writing rules, H11 to H13, and the Style defaults, in [AGENTS.md](../AGENTS.md) |
| See what works today, what is missing and what is planned | "Status and known limits" in the [README](../README.md) |
