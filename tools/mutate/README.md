# mutate

Mutation testing for vstiff ([plan](../../docs/plans/plan.md), section 3.4). `mutate` parses a module with the compiler's own parser and lists its mutants: copies of the module with one small edit, such as `+.` turned into `-.`, a comparison flipped or a constant moved a decade. `mutate run` builds each mutant in a scratch copy of the repository and runs `dune build @runtest`. A mutant that a test notices is killed; one that no test notices survives, and then either a test is missing or the mutant behaves like the original (it is equivalent). [baseline.md](baseline.md) is a run on the kernel and the solver, with a verdict on every survivor.

## Run it

Build the tool once, from the repository root, with dune on your `PATH`:

```sh
dune build ./tools/mutate/mutate.exe
```

On one module (the report goes to standard output, progress to standard error):

```sh
./_build/default/tools/mutate/mutate.exe run src/numerics/newton.ml
```

On all of them:

```sh
./_build/default/tools/mutate/mutate.exe run --results results.tsv src/numerics/*.ml src/*.ml > report.md
```

`--results` writes each outcome to the file as soon as it is known and skips the mutants already in it, so repeating the command resumes a run that was interrupted (Ctrl-C stops it cleanly); delete some lines of the file to run those mutants again. The scratch copies are made in the system's temporary directory (`TMPDIR`), never in the repository, and removed at the end. [baseline.md](baseline.md) records how long the run over twelve modules took.

| Option | Meaning |
|---|---|
| `-j N` | N scratch copies, and as many mutants at a time (default 4) |
| `--root DIR` | the repository (default `.`); the files are relative to it |
| `--results FILE` | record outcomes in FILE, and skip those already there |
| `--equivalent FILE` | mutants known to be equivalent: an id, then the reason, per line; `#` starts a comment. They are reported, not run |
| `--only TEXT` | only the mutants whose id contains TEXT |
| `--timeout-factor X` | a mutant may take X times the baseline run (default 3) |
| `--keep` | leave the scratch copies in place |

## What a run does

1. It copies the repository once per worker, leaving out what dune ignores (`_build`, `.git`).
2. It runs `dune build @runtest` on each unmodified copy, which must pass; the slowest run is the baseline. The time limit of a mutant is three times the baseline, and at least 5 s. It then runs the tests on the modules printed unchanged, to check that the printer alone changes nothing.
3. It applies one mutant at a time to a copy. `dune build @check` fails: stillborn. `dune build @runtest` fails, or runs past the limit: killed. It passes: survived.

Mutants are built with `--profile release`, so a warning (an unused variable that the edit leaves behind) does not make one stillborn, and without dune's cache.

## Look at a mutant

Every mutant has an id, `file:line:col:operator`, which stays valid while the file is unchanged. `list` prints the mutants of a module, with the expression before and after the edit, and `apply` prints the module with one mutant, to build in a scratch copy (the compiler prints it, so its layout is not the source's):

```sh
./_build/default/tools/mutate/mutate.exe list src/halving.ml
./_build/default/tools/mutate/mutate.exe apply src/halving.ml:36:10:if_swap
```

`mutant.mli` names the operators. Swapping two arguments of one type, which plan 3.4 also lists, is not among them: it needs types, and `mutate` reads only the syntax.

## Resolving a survivor

Add a test line that kills it, or list it as equivalent with the reason. The ids hold a line and a column, so an entry goes stale when its module changes; run again to see which.

## The tool's own test

`dune build @runtest` runs `selftest`, which prints the mutants of [fixture.ml](fixture.ml) and the other pure parts of the tool, and dune compares that with [selftest.expected](selftest.expected). `mutate run` is not part of `dune build @runtest`. The tool's effects are in `Files` (the file system), `Child` (processes and the clock) and the two main modules; the rest is pure.
