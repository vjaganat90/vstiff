# vstiff

GPL-3.0-only. Work happens on night/wip. Do not merge.

vstiff is a small OCaml library (OCaml 5.5, standard library only, built with dune) that integrates stiff ordinary differential equations with a variable-step BDF2 method.

An ordinary differential equation (ODE) `y' = f(t, y)` says how fast a state `y` (a vector of numbers, for example three concentrations) changes with time `t`; integrating it means computing the state at a later time from the state at the start. It is *stiff* when some parts of the state change on very short time scales while the part you care about changes slowly. Explicit methods then need tiny steps to stay stable, even when accuracy alone would allow long ones. vstiff uses implicit methods, backward Euler and BDF2 (backward differentiation formulas): every step solves an equation for the new state with Newton's method, and stays stable at long steps. The ideas are explained from scratch in [docs/numerics/01-odes-and-stiffness.md](docs/numerics/01-odes-and-stiffness.md) and the five chapters after it; [docs/glossary.md](docs/glossary.md) defines the terms.

## What it provides

- **Contracts** (OCaml module types) in [src/ode.mli](src/ode.mli): `Ode.Method` (one step), `Ode.Embedded` (a step plus an error estimate) and `Ode.Controller` (a step-size policy), with the types `Ode.problem`, `Ode.point` and `Ode.rhs`. Methods and controllers are modules that implement them.
- **Methods:** `Bdf1` (backward Euler, an `Ode.Method`) and `Bdf2` (variable-step BDF2, an `Ode.Embedded` whose first step is backward Euler).
- **Controller:** `Halving` rejects a step whose scaled error estimate exceeds `tol` (or whose solve fails, or that is too short to move `t`), halves the step that failed, and doubles the step after three accepts in a row, up to `dt_max`.
- **Drivers:** `Stepper.fixed` takes equal steps with any `Ode.Method`; `Adaptive.integrate` takes adaptive steps with any `Ode.Embedded` and `Ode.Controller`, lands exactly on `t_end` and snaps every other step to the floats, `h = (t + dt) - t`, so that the state advances by exactly what the clock does.
- **Named failures:** numerical trouble comes back as `Error` of `Fail.t` (`Diverged`, `StepRejected n` or `Nan`), never as an exception. Only `Check` raises on purpose: `Invalid_argument`, for arguments that make no sense.
- **Building blocks:** `Newton`, `Jac` (forward-difference Jacobians), `Linalg`, `Stage`, `Vec`, `Clock` (how short a step `t` can still resolve) and `Instrument` to count right-hand-side calls.

## A first program

Robertson's chemical kinetics is a classic stiff problem: three species with rate constants from 0.04 to 3e7. This program integrates it from `t = 0` to `t = 1e4`:

```ocaml
open Vstiff

let rhs _t y =
  let a = 0.04 *. y.(0) and b = 1e4 *. y.(1) *. y.(2) and c = 3e7 *. y.(1) *. y.(1) in
  [| b -. a; a -. b -. c; c |]

let problem = { Ode.rhs; t0 = 0.; t_end = 1e4; y0 = [| 1.; 0.; 0. |] }

let () =
  match Adaptive.integrate (module Bdf2) (module Halving) ~tol:1e-6 problem with
  | Ok s ->
      Printf.printf "t = %g\n" s.t;
      Array.iteri (fun i v -> Printf.printf "y%d = %.5f\n" (i + 1) v) s.y;
      Printf.printf "accepted %d, rejected %d\n" s.stats.accepted_steps s.stats.rejected_steps
  | Error e -> Printf.printf "failed: %s\n" (Fail.to_string e)
```

`Adaptive.integrate` takes a method and a controller as modules, then `~tol` and the problem, and returns `Ok solution` or `Error failure`; the solution holds the final `t`, `y` and the controller's `stats`. Optional `?dt0`, `?dt_max` and `?max_rejects` have defaults ([src/adaptive.mli](src/adaptive.mli)). Only the final state comes back, not the path. If the syntax is new, [docs/ocaml.md](docs/ocaml.md) explains each construct.

To run it, create the probe project: a throwaway dune project outside the repository that links the library's `src/` and the corpus problems, with your program saved as `p/probe.ml`. From the repository root, with your opam switch active:

```sh
REPO=$(pwd)
mkdir -p ../vstiff-scratch/p && cd ../vstiff-scratch
printf '(lang dune 3.0)\n' > dune-project
ln -sfn "$REPO/src" src
ln -sfn "$REPO/test/problems.ml" p/problems.ml   # the corpus problems, which the exercises use
printf '(executable (name probe) (libraries vstiff))\n' > p/dune
# save the program above as p/probe.ml, then build and run it
dune build --root . ./p/probe.exe && ./_build/default/p/probe.exe
```

(`--root .` keeps dune from treating an enclosing directory as the project.) The first line is `t = 10000`, and `y1` prints as `0.10730`, the reference value in [test/refs.ml](test/refs.ml) to five decimals (the corpus pins the error below `1e-5`). The step counts on the last line are yours to measure: they are not in the expected files and depend on the tolerance. [docs/exercises.md](docs/exercises.md) builds on this setup.

## Quick start

You need OCaml 5.5.0 and dune 3.x (`ocaml -version` and `dune --version` show what you have; [docs/ocaml.md](docs/ocaml.md) explains how to install them) on macOS, Linux or WSL (Windows Subsystem for Linux): the shell recipes in these documents use `ln`, `tar`, `git` and `perl`. `Stepper.fixed` and `Adaptive.integrate` take modules as arguments (modular explicits, explained in [docs/ocaml.md](docs/ocaml.md)), which older OCaml releases reject. Nothing else is needed. `<repository-url>` is the address of the project's repository, which you were given with your access.

```sh
git clone <repository-url> vstiff && cd vstiff
git switch night/wip
dune build
dune runtest
```

`dune build` compiles the library and the tests and also runs the two test programs to record their output, so it is not instant. `dune runtest` compares that output with [test/corpus.expected](test/corpus.expected) and [test/soak.expected](test/soak.expected): silence and exit status 0 mean it matches, otherwise dune prints a diff (`-` is the expected line, `+` what was printed) and exits with status 1. A failing test is information: never run `dune promote` to silence it ([CONTRIBUTING.md](CONTRIBUTING.md) says when promoting is legitimate).

## Layout

```text
src/        the library; every module has an .mli, ode.mli is interface only
  fail  vec  linalg  clock   failures, vectors, the linear solve, the resolution of time
  newton  jac  stage         the equation inside every implicit step
  ode.mli                    the contracts
  bdf1  bdf2  halving        methods and the step-size controller
  stepper  adaptive          the two drivers
  check  instrument          the only effects: raising, counting
test/       corpus.ml (the cases) with corpus.expected, soak.ml with soak.expected,
            problems.ml, refs.ml (reference values), guard.ml, report.ml
docs/       documentation; docs/README.md is the index
```

## Where to start reading

- **New to OCaml or numerical analysis?** Follow [docs/onboarding.md](docs/onboarding.md), ten working days from setup to a first contribution.
- **Looking for something?** [docs/README.md](docs/README.md) indexes every document, with reading orders and a one-hour path.
- **Want the code first?** Read the `.mli` files of `src/` in the order of [docs/README.md](docs/README.md), then `test/corpus.ml` and `test/corpus.expected`.
- **Ready to change something?** Read [CONTRIBUTING.md](CONTRIBUTING.md), then pick an item from [docs/exercises.md](docs/exercises.md).

## Status and limitations

vstiff is a work in progress. The known gaps double as starter contributions ([docs/exercises.md](docs/exercises.md)):

- Newton rebuilds the finite-difference Jacobian at every iteration (`n + 1` right-hand-side calls each time); production codes reuse it across iterations and steps.
- `Halving` is halve and double. Production controllers scale the step by a safety factor times `(tol / err)^(1/(p+1))`, where `p` is the order of the method whose error is estimated; that needs the error estimate when a step is accepted, which `Ode.Controller.accepted` does not receive, so it needs a contract change. The gap between backward Euler and BDF2 is a conservative estimate for BDF2.
- One mixed absolute and relative weight `1 / (1 + |y_i|)`, no separate `rtol` and `atol`: tiny components such as Robertson's `y2` are controlled only loosely.
- `Halving` counts every rejection together, whatever its reason; `Instrument` counts right-hand-side calls only, so Newton iterations per step are invisible from outside.
- No dense output: the drivers return only the final state. `Linalg.solve` is dense Gaussian elimination, meant for small systems.
- Tests: the backward Euler canary takes 500,000 steps and dominates the test time; van der Pol accuracy is not checked against a reference; a transposed Jacobian stalls the corpus run instead of failing one line; the soak test is a tripwire, not extra coverage.

## License

GNU General Public License, version 3 only (SPDX license identifier `GPL-3.0-only`). The full text is in [LICENSE](LICENSE).
