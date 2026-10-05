# vstiff

A stiff ODE solver in OCaml that puts correctness first, using only the standard library.

An ordinary differential equation (ODE) $y' = f(t, y)$ says how fast a state $y$ (a vector of
numbers, such as three concentrations) changes with time $t$; solving it means computing the
state at a later time from the state at the start. It is *stiff* when part of the state changes on
a very short time scale while the part you care about changes slowly, so an explicit method, which
computes each new state directly from values it already has, needs tiny steps just to stay stable.
vstiff uses implicit methods, which solve an equation for each new state and stay stable at long
steps, and it estimates the derivative matrices that equation needs (Jacobians) from $f$ alone,
so $f$ is the only code you write.

vstiff is early work: the API will change, it returns only the final state, and it is slow.
Newton's method rebuilds the Jacobian at every iteration, at a cost of $n + 1$ calls of $f$ for
a state of $n$ components, and the controller only halves and doubles the step instead of sizing
it from the error estimate. It suits small systems and learning how a stiff solver works; for speed
or large systems use a mature solver such as SciPy (`solve_ivp`) or SUNDIALS. The status table
below says what is missing and what is planned.

## What it does today

- **Method and drivers.** Variable-step BDF2 (a backward differentiation formula of order 2) with a
  backward Euler first step; each step solves its equation by damped Newton's method with a
  forward-difference Jacobian. `Adaptive.integrate` chooses step sizes with a controller
  (`Halving`); `Stepper.fixed` takes equal steps.
- **Long steps on stiff problems.** The test corpus has a stiff canary: three independent decays
  with rates $1$, $100$ and $10^4$, so explicit Euler is stable only for steps up to
  $2 / 10^4 = 2 \times 10^{-4}$. BDF2 with a fixed step of $10^{-3}$ ends within $10^{-6}$ of
  the exact solution at $t = 1$. On Robertson's chemical kinetics (the first program below) at
  `~tol:1e-6`, $y_1$ at $t = 10^4$ is within $10^{-5}$ of a reference computed with three
  SciPy solvers that agree to $4 \times 10^{-12}$. Van der Pol with $\mu = 1000$ runs to
  $t = 2000$ at `~tol:1e-4` and ends in a finite state; no reference checks its accuracy yet.
- **Failures are values.** A solve that cannot finish returns `Error` of `Fail.t` (`Diverged`,
  `StepRejected n` or `Nan`); arguments that make no sense raise `Invalid_argument`.

## Quick start

Install the toolchain first ([docs/README.md](docs/README.md), "Before you start"), then:

```sh
git clone https://github.com/vjaganat90/vstiff
cd vstiff
dune build     # compiles both libraries and the tests, and records what the test programs print
dune runtest   # compares that output with the .expected files; silent when every check passes
```

`dune runtest` prints nothing when every check passes. If it prints a diff, a check failed:
[docs/testing.md](docs/testing.md) ("Reading a failure") shows how to read it. Do not run
`dune promote` to silence it; [AGENTS.md](AGENTS.md) says when promoting is right.

## A first program

Robertson's chemical kinetics is a classic stiff problem: three species with rate constants from
$0.04$ to $3 \times 10^7$. This program follows it from $t = 0$ to $t = 10^4$:

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

`rhs t y` must return a new array and leave `y` alone ([src/ode.mli](src/ode.mli) says why).
`Adaptive.integrate` takes a method and a controller as modules (modular explicits, see
[docs/ocaml.md](docs/ocaml.md)), then `~tol` and the problem; a method or controller of your own
plugs in the same way, by implementing a contract from `src/ode.mli`. It returns `Ok` with the final
`t`, the state `y` and the step counts, or `Error`; `~tol` limits each step's error estimate, not
the final error.

To run it, make a throwaway dune project next to your clone that links the clone's `src/`, so that
nothing is written inside the repository. From the repository root, with your opam switch active:

```sh
mkdir -p ../vstiff-scratch/p && ln -sfn "$PWD/src" ../vstiff-scratch/src && cd ../vstiff-scratch
printf '(lang dune 3.0)\n' > dune-project
printf '(executable (name probe) (libraries vstiff))\n' > p/dune
```

Save the program above as `p/probe.ml`, then build and run it (`--root .` makes dune ignore any
enclosing dune project):

```sh
dune build --root . ./p/probe.exe && ./_build/default/p/probe.exe
```

It prints `t = 10000`, `y1 = 0.10730`, `y2 = 0.00000` and `y3 = 0.89270`, then the numbers of
accepted and rejected steps, which no test pins: change `~tol` and watch them move. $y_1$ agrees
to five decimals with the reference in [test/refs.ml](test/refs.ml).

## How correctness is checked

No solver is right on every input, so vstiff aims to state what it promises and to check each
promise. The checks are *expect tests*: three programs, the corpus, the soak test and the properties, print one line
per case, and `dune runtest` compares what they print with the `.expected` file beside each. The
corpus ([test/corpus.expected](test/corpus.expected)) has 47 lines, one per case (Newton's method,
Jacobians, the order of BDF2, the stiff problems above, step control, the named failures); a line
reads like `max error 3.68e-07 < 1e-06: true`. The soak test ([test/soak.ml](test/soak.ml)) runs
four of the cases ten times; every round must pass and equal the first. The properties
([test/props.ml](test/props.ml)) run claims such as "the backward error of `Linalg.solve` stays
below $`n\varepsilon`$" on thousands of generated cases, and shrink a failing case. A check is a
bound, a typed outcome, a value printed to the digits its bound needs, or a pin (a count or a digit
string that fixes one rule of the algorithm), so that a compiler or a platform changing the last
bits of a result fails no test (rule H1 in [AGENTS.md](AGENTS.md)).

## Status and known limits

The Planned column says what is intended, not what exists.

| Area | Today | Planned |
|---|---|---|
| Method, step size | variable-step BDF2; halve a rejected step, double after three accepts | variable-order BDF (orders 1 to 5); Radau IIA (an implicit Runge-Kutta method); steps sized from the error estimate |
| Tolerance | one `~tol` on each step's error estimate (absolute for small components, relative for large), not the final error | separate relative and absolute tolerances, $\mathrm{rtol}$ and $\mathrm{atol}$; an optional global error estimate |
| Jacobian, linear algebra | forward differences rebuilt at every Newton iteration; dense Gaussian elimination | reuse; a Jacobian you supply; banded and sparse |
| Problems, result | $y' = f(t, y)$ forward in time; the final state only | index-1 differential-algebraic equations (DAEs); backward in time; dense output; events |
| Assurance | the expect tests above | property tests, mutation testing, proofs in Rocq (a theorem prover), comparisons with SciPy and SUNDIALS |

Two gaps in the tests are known: van der Pol's accuracy is not checked against a reference, and no
line pins the constant of Newton's line search. There is no opam package yet, so use vstiff from a
clone.

## Layout

```text
src/numerics/   kernel library numerics: Fail, Vec, Linalg, Newton, Jac; knows nothing about ODEs
src/            solver library vstiff: ode.mli holds the contracts, vstiff.mli the public API
test/           corpus.ml, soak.ml and props.ml with their .expected files, and their helper modules
docs/           README.md is the index; numerics/ explains the mathematics
```

## Where to go next

- **Learning?** [docs/README.md](docs/README.md) indexes the documents, which explain ODEs,
  numerical methods and OCaml from the start, and says what to know before you begin. Start with
  [ODEs and stiffness](docs/numerics/01-odes-and-stiffness.md) and the [glossary](docs/glossary.md),
  or follow the ten-working-day plan in [docs/onboarding.md](docs/onboarding.md).
- **Reading or changing the code?** Every library module has an `.mli` that documents it.
  [docs/architecture.md](docs/architecture.md) maps the modules, and
  [docs/testing.md](docs/testing.md) explains the tests. Contributions are welcome:
  [AGENTS.md](AGENTS.md) has the rules every change follows, and starter contributions end
  [docs/exercises.md](docs/exercises.md). To propose a change, fork the repository, branch from
  `main` and open a pull request against it; questions and bug reports go in the
  [issue tracker](https://github.com/vjaganat90/vstiff/issues).

## License

GNU General Public License, version 3 only (`GPL-3.0-only`); the text is in [LICENSE](LICENSE).
