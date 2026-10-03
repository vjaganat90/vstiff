# bench

The work-precision bench of the [roadmap](../docs/plans/roadmap.md) (4.3 and milestone M1): vstiff's adaptive driver on Robertson, HIRES, van der Pol (mu = 1000, t = 2000) and the Brusselator on 40 cells (n = 80), at the tolerances of its Appendix B. Around it sit the references computed with scipy, the scipy rows of the same tables, and a check that the OCaml and the Python problems are the same. Everything runs from the repository root. `dune build @runtest` does not run the bench; it runs `bench/test`, which checks the pure modules in a fraction of a second.

## Run vstiff

```sh
dune build ./bench/bench.exe
./_build/default/bench/bench.exe
```

Each row is one problem and tolerance (`rtol = atol = tol`), run with `Adaptive.integrate (module Bdf2) (module Halving)` and the other arguments at their defaults. It prints the accepted and rejected steps, the right-hand-side calls (every call, the finite-difference Jacobians included, counted through an `Instrument`-wrapped rhs), the error `max_i |y_i - ref_i| / (1 + |ref_i|)` of Appendix B, and `scd = -log10 max_i |y_i - ref_i| / |ref_i|` against the stored reference, then a verdict against the golden table. The whole run takes about 4 s of CPU on an idle Apple M5; more on a slower or busier machine.

| Flag | Effect |
|---|---|
| `--cpu` | adds the CPU time of each row; it is not reproducible, so nothing compares it |
| `--only NAME` | runs one problem: `robertson`, `hires`, `van_der_pol` or `brusselator_80` |
| `--csv` | prints the rows as CSV, with the columns of `bench/results/scipy.csv` (`nfev` and `njev` empty), and no verdicts |
| `--pin` | prints the source of `bench/golden.ml` for this run |
| `--transcription` | prints the right-hand sides at seeded points, for the transcription check below |

## The golden table

`bench/golden.ml` pins the figures of one run, and every later run is judged against it as the roadmap's gate has it: the scd must not be lower, and the cost, the right-hand-side calls, not more than 10 % higher. Both are bands, so they hold from one platform to the next: the scd may lose 0.05 digits (the pin is rounded), and a count may move with the platform's rounding (arm64 fuses a multiplication into the addition after it, x86_64 does not: with the fusion compiled away in a scratch copy, two of the ten rows moved, by 3 and 10 calls). The exit status is 1 unless every verdict is `ok`; a row the table does not have is `UNPINNED`, and a run that fails is a regression. Steps, rejections and the error are printed for the reader and not compared.

A deliberate change of behaviour re-pins the table in the same commit, with the old and the new figures in its message (CONTRIBUTING, "Expect files are the contract"):

```sh
dune build ./bench/bench.exe
./_build/default/bench/bench.exe --pin > bench/golden.ml
```

Rebuild afterwards: until `dune build ./bench/bench.exe` runs again the executable still carries the old table.

## References

`bench/reference.ml` holds the state of each problem at `t_end`, computed outside vstiff with scipy's Radau (rtol 1e-13), BDF and LSODA (rtol 1e-12), atol 1e-20 and finite-difference Jacobians, keeping only the digits on which all the solvers that succeeded agree (roadmap 4.2). The module says, next to every value, the solvers, versions and tolerances, how many digits agree and the largest difference between the solvers. scipy's BDF fails on van der Pol at that tolerance, so Radau and LSODA cross-check each other there. `bench/compare/references.json` keeps every solver's raw values. Both files are generated; never edit them by hand:

```sh
python3 bench/compare/references.py            # recompute, rewrite both files
python3 bench/compare/references.py --check    # recompute, verify the stored files, write nothing
```

## The scipy rows

```sh
python3 bench/compare/scipy_rows.py > bench/results/scipy.csv
```

scipy's BDF, Radau and LSODA on the same problems and tolerances, with `rtol = atol` and finite-difference Jacobians. `rhs` counts every call of the right-hand side through a wrapper, the finite-difference Jacobian included, as vstiff's counter does; `nfev` is what scipy reports, which for BDF and Radau leaves the Jacobian out (roadmap 2.0). Columns: `problem, solver, rtol, steps, rejected, rhs, nfev, njev, error, scd`; `rejected` stays empty, scipy does not report it. The counts are those of macOS on arm64, where the roadmap's Appendix B was measured: scipy's BDF and LSODA pick their steps from the last bits of the right-hand side, so another libm can move a row by a step or more (`bench/compare/problems.py` says where).

## The transcription check

Every problem exists in OCaml (`bench/hires.ml`, `bench/brusselator.ml`, `test/problems/problems.ml`) and in Python (`bench/compare/problems.py`), written apart. Before any comparison the two right-hand sides must agree at 100 seeded points per problem, or a mistranscribed problem would look like a solver failure:

```sh
dune build ./bench/bench.exe
./_build/default/bench/bench.exe --transcription | python3 bench/compare/transcription.py
```

The script requires 1e-15, relative to the magnitude of the terms of each component; its docstring says why that and not the result. A coefficient that is off by 1e-10 of itself, or more, is caught.

## Python

Python 3.13 with the two packages of `bench/compare/requirements.txt`, which pins exact versions. Dune ignores `bench/compare/`.

```sh
python3 -m venv bench/compare/.venv
bench/compare/.venv/bin/pip install -r bench/compare/requirements.txt
bench/compare/.venv/bin/python bench/compare/scipy_rows.py
```
