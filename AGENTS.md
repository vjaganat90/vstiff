# AGENTS.md: how to change vstiff

These rules are for anyone who changes vstiff, coding agents and people alike. Read this file
before your first change, and check every change against it.

There are two kinds of rule:

- **Hard rules** (H1 to H17) are absolute. A change that breaks one is not merged. The sections
  after the design defaults say how to meet them, and they bind the same way.
- **Design defaults** describe how vstiff is built. Depart from one when the result is clearly
  cleaner, and say why in the pull request.

When a hard rule stands in the way of the work, stop and ask a maintainer to change the rule.
Never work around it.

## The project

vstiff is a stiff ODE solver in OCaml that puts correctness first.

- [README.md](README.md) says what vstiff does today and what is planned, and its "Layout" section
  maps the repository.
- [docs/README.md](docs/README.md) lists, once, what to know and to install before reading the
  code or the documents.

## Hard rules

### Correctness

- **H1. Correct, not bug-compatible.** Every check is one of four things:
  - a correctness bound;
  - a typed outcome;
  - a value printed to the digits its bound makes meaningful: a boolean, a few significant
    digits, or as many decimals as a tolerance fixes, always short of the last digits where
    machines differ;
  - a pin: a count or a digit string that fixes one rule of the algorithm, such as step and
    rejection counts or right-hand-side calls.

  A refactor, a compiler or a platform may change the last bits of a result, through fused
  multiply-add, the order of a sum or the C library. No check may depend on those bits.
  Determinism within one build is required: the same binary, given the same inputs, returns the
  same outputs.
- **H2. Never weaken a check.**
  - Never loosen a bound. Never delete a passing expect line. Never `dune promote` a failure.
  - A failing test is information: find out why it fails before you touch it.
  - A deliberate change of algorithm, such as a new method, estimate or controller, re-pins the
    lines it moves, in the same commit. The commit message gives the old and the new values, and
    the correctness bounds stay as they are.
  - A line may also move only because the last bits of a result changed. Then it printed more
    than its check needs (H1). Rewrite it to print the check, in a commit of its own whose message
    gives the old and the new line. The bound stays.
- **H3. Test every behaviour change.**
  - Every behaviour change ships with a test, and every fixed bug with the line that would have
    caught it.
  - Watch a new test fail before you rely on it. Break what it protects in a scratch copy of the
    working tree, outside the repository ([docs/exercises.md](docs/exercises.md), "Setup"), and
    see the test go red.
  - A fix's commit message gives the output its new line shows without the fix.
- **H4. References come from outside.**
  - [test/refs.ml](test/refs.ml) is never edited.
  - A new reference value is computed outside vstiff, never with it. At least two independent
    solvers (different algorithms or implementations) must agree on every digit kept.
  - The value goes into a new module next to the code that uses it (for the corpus,
    `test/refs_<problem>.ml`). Record the solvers, their versions, the tolerances and the
    agreement next to the value.
- **H5. Numerical failures are values; invalid arguments raise; every run ends.**
  - A solve or a step that cannot succeed returns `Error` of `Fail.t`. Nothing returns `Ok` with a
    non-finite value.
  - The drivers check their arguments with `Check`, which raises `Invalid_argument`. The kernel
    does not check: its `.mli` files state each precondition.
  - Every run ends. A controller ends every run of rejections with an `Error`, and the termination
    of `Adaptive.integrate` rests on it ([docs/numerics/05-step-control.md](docs/numerics/05-step-control.md)).
- **H6. Effects live in named modules; mutation stays inside the call.**
  - In the library, only `Check` raises on purpose and only `Instrument` keeps state from one call
    to the next. In the tests, only `Guard` raises or catches and only `Report` prints. The effects
    map below has the details.
  - A new effect goes into one of these modules, or into a new module that exists for it.
  - Mutation inside a call is allowed: `ref`s, loops, writes to arrays the call allocated, work on
    a private copy. It must not be visible outside: no function writes an argument, an array it
    returned earlier, or anything its caller can see, and a result never changes once returned.
  - Mutate only when it pays: a pure version is replaced by a mutating one only if the mutating one
    is at least 5% faster, averaged over 3 runs of each. Anything under 5% is run noise. The pull
    request gives both averages.
- **H7. The corpus exercises the solver the way users will.**
  - No integrator is given an analytic Jacobian: the problems in `test/problems.ml` supply none,
    the solver works from function values alone, and the corpus keeps exercising `Jac`.
  - The `newton` cases hand `Newton.solve` their own Jacobians, and the `jacobian` cases compare
    `Jac.forward` with analytic ones. Neither runs an integrator.
- **H8. A method's coefficients have one source.**
  - No comment or document gives the value of a coefficient, a table of them, or the solved
    closed form of a weight.
  - The source is the code that computes them, today `Bdf2.coeffs`. A document states the
    equations the coefficients satisfy, and points to that code.

### Structure

- **H9. Standard library only.**
  - The libraries and the tests use the OCaml standard library and nothing else.
  - A tool, meaning an executable that reads the sources (a mutation tester, a tripwire), may also
    use `compiler-libs`, which ships with the compiler.
  - A program outside the dune build may run another solver as an oracle, to compute a reference
    or a comparison. What it produces enters the repository as data, under H4.
  - Any other dependency needs a reason stated in its pull request, and a maintainer's agreement.
- **H10. Two libraries, one public API.**
  - The kernel in `src/numerics/` depends on the standard library alone, and calls no solver
    module.
  - The solver's public API is exactly what [src/vstiff.ml](src/vstiff.ml) and
    [src/vstiff.mli](src/vstiff.mli) list.
  - Public signatures write vectors as `float array`, never `Vec.t`. The one kernel type they
    name is `Fail.t`, which the solver re-exports as `Vstiff.Fail`.
- **H17. Code is written for vstiff.**
  - No source is copied from another solver or library.
  - An algorithm taken from a publication is implemented afresh, and cited where it is used: in an
    `.mli` header or in `docs/`.
  - Adapting someone else's code needs a maintainer's agreement first, and keeps its licence
    notice. vstiff is GPL-3.0-only.

### Writing

- **H11. Write about the project, for its readers.**
  - Documents, comments, commit messages and pull requests describe the code and the reasons for
    it.
  - They never describe how a change was requested or produced: no conversations, no "as
    discussed", no session or tool notes, no credit to an agent, an assistant or a tool, and no
    instructions to whoever writes the text.
  - No commit hash or other raw commit reference appears in a document, a comment or code. Name
    the change by what it did, or point to the test line that pins it.
- **H12. Every statement is checked, and every number has a source.**
  - A statement about the code was checked against the code. A statement about behaviour was
    seen in a run.
  - A number comes from one of five sources:
    - an expected file;
    - a reference module (H4);
    - a derivation shown in the text;
    - a probe the reader can run;
    - a publication the text cites.

    A figure with none of these is marked as unverified.
  - Links between project files are relative paths to files, never `#anchors`. An external source
    is linked by its full URL.
  - Every OCaml snippet in an `ocaml` block compiles. A sketch that does not compile goes in a
    `text` block.
- **H13. Prerequisites are stated once.**
  - What a reader needs to know and to install lives in [docs/README.md](docs/README.md),
    "Before you start".
  - A document links there. It neither opens with its own list of prerequisites nor repeats one
    further down.

### Git

- **H14. Protect `main`.**
  - Never push to `main`, and never rewrite it.
  - Cut a branch from the latest `main`; in a stack, from the branch below.
  - Only a maintainer merges, after review. A coding agent never merges a pull request.
- **H15. Use Conventional Commits** (below), one logical change per commit.
  - Every commit passes the checks of "Done means".
  - A fix and the test that pins it go in the same commit.
- **H16. Title pull requests with an area** (below), never with a commit type.

## Design defaults

The goal is clean, powerful modularity, and every tool OCaml offers is available for it. These
are the defaults. Break one when the alternative is clearly cleaner, and say why in the pull
request.

- **Arrays and speed.** Vectors are `float array`, which OCaml stores as one flat block of
  unboxed doubles. Where it is faster, use the mutation H6 allows: fill a freshly allocated result
  in a loop rather than through `Array.init` or `Array.map` with a closure, and factor a matrix in
  place on a private copy, when the gain clears H6's 5% bar. Keep every output within its checks
  (H1).
- **Contracts and implementations.**
  - Contracts are module types; implementations are modules.
  - Methods and controllers implement `Ode.Method`, `Ode.Embedded` and `Ode.Controller`. They are
    passed to the drivers as modular explicits: `Adaptive.integrate (module Bdf2) (module Halving)
    ~tol problem`. A modular explicit lets a result type depend on the module passed in.
  - Use a functor, a first-class module in a data structure, a GADT or an effect handler when it
    gives the clearer design.
  - A change to a contract or to the public API shows the old and the new signatures in its pull
    request.
- **Sharing types between signatures.** Use `with type` only where callers need the type.
- **Names.**
  - Module types are CamelCase; values are snake_case.
  - `dt` is a step size someone asked for. `h` is a step taken or attempted.
- **Labels.** Use one where two arguments of one type could be swapped (`~at ~h`), or to name a
  literal at a call site (`~tol:1e-6`). Otherwise arguments are positional.
- **`.mli` files.** Every library module has one, and it is the API documentation:
  - a header of 1 to 4 lines;
  - then, per item, what it means, its inputs and outputs, and how it fails: 1 to 3 lines.
  - A driver whose defaults and failure modes are its contract takes as many lines as those
    facts need, as `Adaptive.integrate` does.
- **`.ml` comments.**
  - A comment says why, usually in 1 to 3 lines: a numerical reason, an invariant, a non-obvious
    construct.
  - It never narrates the code or repeats the `.mli`.
  - Depth belongs in `docs/`. A module's comment lines should not outnumber its code lines.
- **Style.**
  - Match the surrounding code: indentation, naming, and lines near 120 columns.
  - Mathematics in Markdown looks as it does in a numerical analysis book, so that an applied
    mathematician reads it at once: every symbol, subscript, power, norm, absolute value, formula
    and number such as $10^{-6}$ is typeset in LaTeX, with upright operators ($\operatorname{diag}$,
    $\max_i$, $\mathrm{tol}$), never left as plain text or in a code span. The step ratio, for
    one, is $\omega = h / h_{\mathrm{prev}}$ and a norm is $\lVert x \rVert_\infty$.
  - The delimiters are standard LaTeX, which editors, GitHub and pandoc all render: inline math
    between single dollar signs, with no space just inside them, and displayed math between double
    dollar signs on lines of their own, with a blank line before and after. Never the GitHub-only
    form with backticks inside the dollars, and never a `math` fence. In a table cell write
    $\lvert x \rvert$, never a bare `|`, which splits the cell.
  - Write the TeX that GitHub passes on intact: `\lt` and `\gt` for the two inequality signs (a
    bare one reaches the renderer as `&lt;`), `\lbrace` and `\rbrace` for braces, and
    `\thinspace` for a thin space, because `\,`, `\;` and `\{` lose their backslash. A row break
    `\\` survives only at the end of a line of a display that spans several lines; anywhere else
    write `\cr`.
    `\lt` and `\gt` are MathJax and KaTeX macros, not LaTeX, so a PDF build through LaTeX must
    define them: `\newcommand{\lt}{<}\newcommand{\gt}{>}`.
  - Code stays code: identifiers and expressions (`Bdf2.coeffs`, `~tol:1e-6`), commands, file
    names, and program or expected-file lines quoted verbatim. A sentence about the mathematical
    quantity is math; a sentence that names the code is code. Math in code comments stays plain
    ASCII, since no renderer reads a comment.
- **Comment pitfalls.** OCaml comments nest, and string literals are lexed inside them, so:
  - keep `(*`, `*)` and double quotes out of comment text;
  - inside a function body, use `(* *)`: a documentation comment there is warning 50, which the
    dev profile treats as an error.

## Build and test

The toolchain, and how to install it, is listed in [docs/README.md](docs/README.md), "Before you
start".

| Command | Use |
|---|---|
| `dune build @check` | Type-check only: the fast loop |
| `dune build` | Compile everything, and run the test programs to record their output |
| `dune runtest` | Diff that output against the `.expected` files; silent when green |
| `dune exec ./test/corpus.exe` | Run the corpus and print its output unfiltered |
| `dune promote` | Record a new case, a re-pin or a reformatted line (H2); never a failure |

- **Warnings.** The dev profile turns warnings into errors. Fix the code, never the flags.
- **Runs that may not finish.** `dune build` and `dune runtest` run both test programs with no
  time limit, and stopping dune can leave a program running. After an edit to `Newton`, `Jac`,
  `Halving`, `Adaptive` or anything they call, build the programs alone, then run each under a
  limit:

  ```sh
  dune build ./test/corpus.exe ./test/soak.exe
  perl -e 'alarm 120; exec @ARGV' ./_build/default/test/corpus.exe
  perl -e 'alarm 120; exec @ARGV' ./_build/default/test/soak.exe
  ```
- **Trying things out.** Experiment outside the repository, in the probe project or in a scratch
  copy of the working tree ([docs/exercises.md](docs/exercises.md), "Setup"). Nothing stray then
  reaches a branch.

## Tests

[docs/testing.md](docs/testing.md) explains the test programs.
[docs/numerics/06-the-corpus.md](docs/numerics/06-the-corpus.md) explains what each corpus line
proves.

- **Expect tests.** `test/corpus.ml` and `test/soak.ml` print one line per case, and dune diffs
  each output against its `.expected` file.
  - A line says what it checks, as in `max error 3.68e-07 < 1e-06: true`.
  - It prints only what its check needs, so the same files pass on every platform (H1).
- **Adding a corpus case.**
  1. Write the comment above the new group first. In one sentence, it says which defect the case
     catches that no existing line catches. If you cannot write that sentence, do not add the
     case.
  2. Put a new problem in `test/problems.ml`, with `exact` when the solution is known. H7 says
     what a problem may carry. A new reference goes in a new module (H4).
  3. Append the case as a new group at the end of the list in `test/corpus.ml`: the documents
     cite line numbers.
  4. If the run might not finish, wrap the right-hand side in `Guard.budget` and the case in
     `Guard.run`.
  5. The diff must show your line and nothing else. Check the line by reasoning, then promote it.
  6. In the same commit, update every document that counts or cites corpus lines, such as
     [docs/numerics/06-the-corpus.md](docs/numerics/06-the-corpus.md) and [README.md](README.md).
- **Adding a soak case.** Write a module of type `Case` in `test/soak.ml`, and wrap an adaptive
  run with `on_budget`.
- **Adding a method or a controller.**
  1. Implement the `Ode` contract in `src/`, with an `.mli` that includes it.
     [src/bdf2.mli](src/bdf2.mli) and [src/halving.mli](src/halving.mli) are the models.
  2. List it in `src/vstiff.ml` and `src/vstiff.mli` to make it public.
  3. Add corpus cases. A controller must end every run of rejections (H5).

## Effects map

| Module | Effect |
|---|---|
| `Check` (library, internal) | Raises `Invalid_argument`, on purpose, for arguments that make no sense |
| `Instrument` (library) | A mutable counter of right-hand-side calls |
| `Guard` (tests) | Raises `Exhausted` past its call budget. Turns `Exhausted` and `Invalid_argument` into text (`run`), or `Exhausted` into `None` (`bounded`) |
| `Report` (tests) | Prints one line per case |
| Everything else | Pure. A standard-library function that raises on misuse, such as an out-of-bounds index, signals a bug in the caller |

## Commits

The format is `type(scope): imperative summary`.

- **The subject line** is under 72 characters, type and scope included. It starts with a
  lowercase letter, though names of modules and values keep their case, and it ends without a
  period.
- **The body** says why.
- **No credit.** No trailer, footer or comment credits a tool or an assistant.

| Type | For |
|---|---|
| `feat` | New behaviour or capability |
| `fix` | A bug fix, together with the line that pins it |
| `perf` | Faster, with every result still within its bounds |
| `refactor` | Structure only: no change of behaviour |
| `test` | Tests only |
| `docs` | Documents and comments only |
| `build` | dune and packaging |
| `ci` | Continuous integration |
| `chore` | Anything else: tooling, housekeeping |

- **Scope.** The scope names a module or an area: `newton`, `adaptive`, `corpus`, `soak`,
  `docs`, `agents`. For example, `fix(halving): stop when a halved step falls below the floor`.
- **One change per commit.** Do not mix a refactor with a behaviour change.
- **Branch names** are `<type>/<short-slug>`, for example `fix/halving-floor`.
- **Rewording a pushed commit.**
  1. Rebuild the commit from its original tree, author and dates.
  2. Check that the trees match.
  3. Push with `--force-with-lease`.

  Never do this on `main` (H14).

## Pull requests

The title is `Area: Sentence-case summary`. The area is one of these words, never a commit type:

| Area | For |
|---|---|
| `Core` | Solver and kernel behaviour |
| `Design` | Interfaces, architecture, refactors |
| `Spec` | Specifications and written proposals |
| `Proof` | Formal verification |
| `Harness` | Tests, properties, mutation testing, tripwires, corpus references |
| `Bench` | Benchmarks, their references, comparisons with other solvers |
| `Guide` | README and documentation |
| `Infra` | CI, build, tooling |

When a change fits more than one area, take the one that names what it does to the repository:

- A change to a contract in `ode.mli` is `Design`; a written proposal for one is `Spec`.
- A reference used by a corpus line is `Harness`; one used by a benchmark is `Bench`.

For example, `Core: Variable-order BDF from the point history`.

- **Body.**
  - **What**: 1 to 3 bullets.
  - **Why**: one line.
  - Leave out local test logs and pass counts. Mention testing only when a reviewer cannot
    repeat it with the commands above.
- **One logical change per pull request.** A larger change becomes a stack:
  - every pull request targets `main` and is built on the branch below it;
  - a problem in a lower pull request is fixed there, and the branches above it are rebased.

## Done means

Every hard rule holds. These are the checks most often missed:

- [ ] `dune build @check`, `dune build` and `dune runtest` pass, with no warnings and no diff.
- [ ] No check was weakened. A moved line was re-pinned for a deliberate change, or rewritten to
      print its check, with the old and the new line in the commit message (H2).
- [ ] Every behaviour change has a test that was seen to fail without the change (H3).
- [ ] Effects are still in their modules. Failures are `Fail.t` values, and every run ends (H5,
      H6).
- [ ] The public API is deliberate. There is no new dependency without agreement, and no copied
      code (H9, H10, H17).
- [ ] Every document and `.mli` comment that the change makes untrue is updated in the same
      change. Links are relative, numbers are sourced, snippets compile, and there is no process
      talk (H11 to H13).
- [ ] Commit messages and the pull request title follow the conventions above (H15, H16).
