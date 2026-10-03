# vstiff documentation

These pages explain vstiff to someone who has never seen the code, has never used OCaml and has not studied numerical analysis. They assume you program in another language (Python or JavaScript, say), remember first-year calculus (derivatives, Taylor series) and a little linear algebra (matrices, solving `A x = b`). Everything else is explained from the start: ODEs and stiffness, floating point, Newton's method, finite differences, BDF methods, step-size control, OCaml and dune, and how the tests work.

Throughout, `t` is time, `y` the state vector, `f` the right-hand side (called `rhs` in code), `h` a step taken or attempted and `dt` a step size someone asked for. The [glossary](glossary.md) lists every term and symbol.

## Index

| Document | What it is |
|---|---|
| [README.md](../README.md) (repository root) | The front page: what vstiff is, a first program, quick start, layout, known limitations. |
| [CONTRIBUTING.md](../CONTRIBUTING.md) (repository root) | How to contribute: setup, the build and test loop, the project's rules with their reasons, code and comment style, adding a test, a method or a controller, commits, a review checklist. |
| [docs/README.md](README.md) | This page: the index, reading orders and a one-hour path. |
| [onboarding.md](onboarding.md) | A ten-working-day plan from zero to a first contribution, with readings, exercises, self-checks and a done-when line for each day. |
| [glossary.md](glossary.md) | Terms and symbols, each with a short definition and where it is explained. |
| [exercises.md](exercises.md) | Graded exercises (warm-ups, reading, probes, breaking the code in a scratch copy), a capstone against the contracts, and starter contributions. |
| [ocaml.md](ocaml.md) | OCaml and tooling primer: opam, dune, the commands you will type, and the language features the code uses, libraries and main modules included. |
| [architecture.md](architecture.md) | The two libraries and the public API, the contracts, the modules, how data flows between them, and where effects live. |
| [testing.md](testing.md) | How the tests are built and run, how to read a failure, the policy for changing expectations, adding a case, probes. |
| [numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md) | ODEs, explicit and backward Euler, order, stability and stiffness. |
| [numerics/02-newton.md](numerics/02-newton.md) | Newton's method, damping and the linear solves inside it. |
| [numerics/03-jacobians-and-floating-point.md](numerics/03-jacobians-and-floating-point.md) | Floating point, finite-difference Jacobians, round-off and NaN. |
| [numerics/04-bdf.md](numerics/04-bdf.md) | BDF1, BDF2, variable steps, the stage equation and stability. |
| [numerics/05-step-control.md](numerics/05-step-control.md) | Error estimates, step-size control, the clock that steps are snapped to (and the grid of the fixed-step driver), and why every run ends. |
| [numerics/06-the-corpus.md](numerics/06-the-corpus.md) | The test problems and what each one proves. |

## Reading orders

**Following the plan.** [onboarding.md](onboarding.md) interleaves everything below over ten working days. It is the recommended route if you will contribute.

**The numerics track.** Read the six chapters in `numerics/` in order, 1 to 6. Each begins by saying what it assumes from the earlier ones and ends with questions to check yourself. You need calculus, plus the little linear algebra that is explained where it appears. The OCaml snippets are optional but worth running. ([onboarding.md](onboarding.md) visits the chapters in the order 2, 3, 1, 4, 5, 6, which their prerequisites allow, to match the code it asks you to read each day.)

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

Then the tests: [test/problems.ml](../test/problems.ml), [test/refs.ml](../test/refs.ml), [test/guard.ml](../test/guard.ml), [test/report.ml](../test/report.ml), [test/corpus.ml](../test/corpus.ml) with [test/corpus.expected](../test/corpus.expected), and [test/soak.ml](../test/soak.ml) with [test/soak.expected](../test/soak.expected).

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

If the toolchain is installed (the [README](../README.md) quick start lists what you need), start `dune build` in a terminal first: it also runs the test programs, so it takes a little while.

1. About 5 minutes: the [README](../README.md), for what the library does and how the repository is laid out.
2. About 15 minutes: [architecture.md](architecture.md), for the contracts, the modules and how data flows between them. Skim the rest.
3. About 15 minutes: the parts of [numerics/01-odes-and-stiffness.md](numerics/01-odes-and-stiffness.md) on stability, stiffness and backward Euler. They explain why the library exists.
4. About 15 minutes: [src/ode.mli](../src/ode.mli), then [src/stage.ml](../src/stage.ml), [src/bdf1.ml](../src/bdf1.ml) and [src/adaptive.ml](../src/adaptive.ml). Between them they show a contract, an implicit step and the loop that drives it.
5. About 5 minutes: [test/corpus.expected](../test/corpus.expected), then `dune runtest`, which prints nothing when everything passes.

Keep the [glossary](glossary.md) open for the terms you meet. To carry on, start [onboarding.md](onboarding.md) at day 1.

## Finding what you need

| You want to | Go to |
|---|---|
| Know what a word or symbol means | [glossary.md](glossary.md) |
| Know what is public API and what is internal | [architecture.md](architecture.md), [src/vstiff.mli](../src/vstiff.mli) |
| Run a piece of the library and see what it does | the probe project in [exercises.md](exercises.md) |
| Understand a failing test | [testing.md](testing.md) |
| Add a test | [testing.md](testing.md) and [CONTRIBUTING.md](../CONTRIBUTING.md) |
| Add a method or a controller | [CONTRIBUTING.md](../CONTRIBUTING.md) and [architecture.md](architecture.md) |
| Find something small to work on | the starter contributions in [exercises.md](exercises.md) |
| Fix or add documentation | the documentation rules in [CONTRIBUTING.md](../CONTRIBUTING.md) |
| See where the project is meant to go, or how parts of it could be proved correct | [plans/plan.md](plans/plan.md), [plans/roadmap.md](plans/roadmap.md), [plans/formal-verification.md](plans/formal-verification.md) |

## Plans

Three documents look ahead instead of describing the code as it is:

| Document | What it is |
|---|---|
| [plans/plan.md](plans/plan.md) | The top-level plan: what done means, and the order of the solver, property-testing, proof and benchmark work. |
| [plans/roadmap.md](plans/roadmap.md) | The roadmap from today's small BDF2 integrator to a general stiff solver, with milestones and decisions. |
| [plans/formal-verification.md](plans/formal-verification.md) | How parts of vstiff could be verified in Rocq with MathComp, with a ranked list of theorems and a pilot plan. |
