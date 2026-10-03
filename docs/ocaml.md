# OCaml and tooling primer

This page covers the OCaml and the tools that vstiff and its tests use, and no more. It assumes you program in another language (Python or JavaScript, say) and have never used OCaml. Part 1 gets the toolchain running. Part 2 tours the language one feature at a time: a tiny example, then where the code uses it. [README.md](README.md) indexes the other pages, [architecture.md](architecture.md) maps the modules and [testing.md](testing.md) explains the tests. Official documentation: <https://ocaml.org/docs> (language and standard library), <https://dune.readthedocs.io> (build system), <https://opam.ocaml.org/doc/> (package manager).

**Reading the snippets.** Every `ocaml` snippet compiles on its own; `text` blocks hold excerpts and compiler messages. A comment `(* val f : ... *)` shows the type the compiler infers, and `(* = v *)` a value. You seldom write types yourself, but reading them is the most useful OCaml skill: a function's type says most of what it does, and editors show it on hover. To try a snippet, save it and run `ocaml snippet.ml` (the plain toplevel that comes with the compiler), or type `#use "snippet.ml";;` inside `ocaml` to see the type of every definition. `utop` or a probe ([testing.md](testing.md)) work too, but the main module of a probe rejects unused top-level definitions, so print what you define.

## Part 1. Tooling

| Piece | What it is | Analogy |
|---|---|---|
| `ocamlopt`, `ocamlc` | the compilers (native code, bytecode); dune runs them | `gcc` or `javac`: compiled ahead of time, where CPython and Node interpret |
| opam | the package manager, which also installs compilers; a **switch** is one self-contained installation | pip + venv + pyenv, or npm + nvm |
| dune | the build system: reads `dune-project` and `dune` files, orders the compilation, writes `_build/` | make or cmake, or npm scripts + a bundler |
| ocaml-lsp-server | editor integration (a language server) | pyright, tsserver |
| utop | an interactive prompt (optional) | the Python or Node REPL |

```sh
opam init                        # first time only: creates ~/.opam
opam switch create 5.5.0         # a switch holding OCaml 5.5.0
eval $(opam env)                 # point this shell at it; repeat in each new terminal
opam install dune                # the build system
opam install ocaml-lsp-server    # optional: editor support
opam install utop                # optional: interactive prompt
ocaml -version                   # The OCaml toplevel, version 5.5.0
dune --version                   # 3.x
```

These commands assume macOS, Linux or WSL (Windows Subsystem for Linux). The project is built and tested with OCaml 5.5.0; `Stepper.fixed` and `Adaptive.integrate` use modular explicits (below), which older releases reject. It uses only the standard library, so nothing else has to be installed. `opam switch list` shows your switches; `opam exec --switch=NAME -- dune runtest` runs one command in a switch without changing the shell (NAME may be the path of a switch kept in a directory).

### What dune reads

Three small files configure the build. They are lists in parentheses, and `;` starts a comment (the copies in the repository may carry comments, left out here).

```text
dune-project   (lang dune 3.0)           the dune language version, not the compiler's
               (name vstiff)
src/dune       (library
                (name vstiff)
                (modules_without_implementation ode))
test/dune      (tests
                (names corpus soak)
                (libraries vstiff))
```

`(lang dune 3.0)` has to be the first line of `dune-project`: dune rejects the file if a comment comes before it. Every `.ml` and `.mli` in `src/` becomes a module of one library called `vstiff`; there is no list of files, and dune orders the compilation from the module names each file mentions. `test/dune` declares two test programs, `corpus` and `soak` (their main modules are `test/corpus.ml` and `test/soak.ml`), linked against the library; the other `.ml` files in `test/` are ordinary modules both can use. Each program has an `.expected` file next to it.

### Everyday commands

Run these from the repository root. Dune builds in its **dev profile** by default, in which many compiler warnings are errors (see "Warnings are errors"). Everything lands under `_build/`, which git ignores.

| Command | What it does |
|---|---|
| `dune build` | Compiles everything. It also runs each test program to record its output under `_build/` (when that is out of date) but compares nothing: a crash fails it, a changed line does not. |
| `dune runtest` | Builds, then compares what each test program printed with its `.expected` file. Silent with status 0 when they match; otherwise a diff and status 1. Dune remembers a pass, so a repeat with nothing changed does nothing; a failing comparison is reported again every time. |
| `dune runtest --force` | Repeats the comparison although it passed. It does not run the programs again: dune keeps their recorded output until something they depend on changes. |
| `dune exec ./test/corpus.exe` | Runs one program and shows its output unfiltered. |
| `dune build @check` | Type-checks without linking or running: the quickest way to see compile errors. (`@check` is an alias, dune's name for a group of targets.) |
| `dune build --watch` | Rebuilds when a file changes (`dune runtest --watch` reruns the tests). |
| `dune promote` | Overwrites `.expected` files with the output of the last run. Read the next section first. |
| `dune clean` | Deletes `_build/`. |
| `dune build @doc-private` | Renders the `.mli` comments as HTML under `_build/default/_doc/_html/`; needs `opam install odoc`. (`@doc` alone builds nothing for a private library.) |

**Reading a failing test.** `dune runtest` prints a diff in which lines starting with `-` are in the expected file and lines starting with `+` were printed. Swapping two rows of the matrix that `Jac.forward` builds, for example, changes among others this line:

```diff
-bdf1 canary t=1 dt=2e-6: max error 3.68e-07 < 1e-06: true
+bdf1 canary t=1 dt=2e-6: max error 3.83e-07 < 1e-06: true
```

[testing.md](testing.md) shows the whole diff and how to read it.

**Why `dune promote` is rarely the right move.** It copies the program's actual output over the `.expected` file, and from then on `dune runtest` is silent whatever the output was, wrong answers included. A red test is information; promoting it away destroys that. The legitimate uses are recording the line of a newly added case and a deliberate, reviewed change of format ([testing.md](testing.md) has the rules).

**Editor and prompt.** Install `ocaml-lsp-server` and use an editor with Language Server support (VS Code with the OCaml Platform extension, Emacs, Vim or Neovim with an LSP client). You want type on hover, jump to definition and inline errors; open the editor at the repository root and run `dune build` once if they do not work at first. The repository has no `.ocamlformat`, so there is no formatter to run: follow the layout of the file you edit. `dune utop src` starts `utop` with the library loaded; type `open Vstiff;;` first (phrases end with `;;`), and `#show Newton;;` prints a module's interface.

## Part 2. The language as this code uses it

### How the code reads

- **Everything is an expression.** `if` and `match` have values and a function returns the value of its body; there is no `return`.
- **Application is juxtaposition.** Write `Vec.scale 2. v`, not `Vec.scale(2., v)`; parentheses only group. A negative literal needs them: `Vec.scale (-1.) fx` in [src/newton.ml](../src/newton.ml).
- **A name is usable only below its definition**, except that `let rec` may call itself. Helpers therefore come first, and modules cannot depend on each other in a cycle.
- `(* ... *)` is a comment; `(** ... *)` is a documentation comment attached to the item before or after it. `x'` is an ordinary name ("x prime", the next x). `a; b` runs `a` (of type `unit`), then `b`. `()` is the only value of `unit`, and a top-level `let () = e` runs `e` when the program starts.

### `let`, functions and shadowing

```ocaml
let area side = side *. side              (* val area : float -> float *)
let hyp2 a b =
  let a2 = a *. a and b2 = b *. b in      (* local names; and binds both at once *)
  a2 +. b2
let x = 1
let y = x + 1                             (* = 2 *)
let x = 10                                (* shadows the first x; y is still 2 *)
let twice f x = f (f x)                   (* val twice : ('a -> 'a) -> 'a -> 'a *)
let eight = twice (fun v -> v *. 2.) 2.   (* fun makes an anonymous function: = 8. *)
```

A `let` never assigns: a later `let` of the same name hides the earlier one. `'a` is a type variable: any type, the same one everywhere it appears. Where: `let ... and ...` in `Halving.rejected`; shadowing in `Stepper.fixed`, where `let* y, history = ... in` rebinds `history`.

### Partial application

```ocaml
let add x y = x +. y                      (* val add : float -> float -> float *)
let increment = add 1.                    (* val increment : float -> float *)
let three = increment 2.                  (* = 3. *)
let scale s = Array.map (( *. ) s)        (* an operator is a function: ( *. ) s multiplies by s *)
let doubled = scale 2. [| 1.; 2. |]       (* = [|2.; 4.|] *)
```

Every function takes one argument: `float -> float -> float` means `float -> (float -> float)`, so giving fewer arguments returns a function that waits for the rest. Write the product operator as `( *. )` with spaces, because the two characters `(*` open a comment. Where: `Vec.scale` and `Vec.norm_inf`; `let f = rhs t in` in `Stage.solve`; `Guard.run reached (fun () -> ...)` in [test/corpus.ml](../test/corpus.ml), a function still waiting for `()`; and `let bdf2_halving = Adaptive.integrate (module Bdf2) (module Halving)`, which fixes the first two arguments.

### Labelled and optional arguments

```ocaml
let ratio ~num ~den = num /. den             (* val ratio : num:float -> den:float -> float *)
let r1 = ratio ~num:1. ~den:4.               (* = 0.25 *)
let r2 = ratio ~den:4. ~num:1.               (* labelled arguments may come in any order *)
let num = 1. and den = 4.
let r3 = ratio ~num ~den                     (* punning: ~num means ~num:num *)

let scale_by ?(factor = 2.) x = factor *. x  (* val scale_by : ?factor:float -> float -> float *)
let six = scale_by 3.                        (* = 6. *)
let thirty = scale_by ~factor:10. 3.         (* = 30. *)
let capped ?limit x =                        (* with no default, limit arrives as a float option *)
  match limit with None -> x | Some l -> Float.min x l
```

**The labels rule.** A label appears only where two arguments of one type could be swapped, or to name a bare literal at a call site: `~y ~err` (two vectors) and `~at ~h` (two floats) in `Ode.Controller`, `~dt:2e-6` and `~tol:1e-6` at call sites. An optional parameter needs an unlabelled one after it, so that the call can settle it: `Adaptive.integrate ?dt0 ?dt_max ?max_rejects ~tol problem` ends with `problem`. Partial application keeps the optional parameters, which is how the corpus calls `bdf2_halving ~dt0:0.5 ~tol:1e-6 problem`. Without a default, an optional parameter is an option inside the function: `Option.value dt_max ~default:(span /. 10.)` in [src/adaptive.ml](../src/adaptive.ml).

### Floats, equality and `Printf`

```ocaml
let q = 7 / 2                            (* = 3: integer division *)
let x = 7. /. 2.                         (* = 3.5 *)
let y = float_of_int q +. x              (* = 6.5: ints and floats never mix *)
let k = Float.to_int 2.9                 (* = 2: truncates *)
let r = Float.round 2.5                  (* = 3.: ties go away from zero *)
let same = [| 1. |] = [| 1. |]           (* = true: = compares contents *)
let nan_eq = nan = nan                   (* = false: nan equals nothing, itself included *)
let line = Printf.sprintf "error %.3e, ok %b, steps %d" 3.116e-06 true 10
(* = "error 3.116e-06, ok true, steps 10" *)
```

`+ - * /` are for `int` and `+. -. *. /.` for `float`; `-. x` negates a float; a literal is a float if it has a point or an exponent (`1.`, `1e-8`). Useful: `exp`, `sqrt`, `Float.abs`, `Float.max`, `Float.epsilon` (about 2.2e-16), `Float.is_finite`. Float division by zero does not raise (`1. /. 0.` is `infinity`, `0. /. 0.` is `nan`), and `Float.max` returns `nan` if an argument is `nan`, so one `nan` entry makes `Vec.norm_inf` return `nan`. `=` is structural equality; `==` is physical (the same object) and rarely wanted. [test/soak.ml](../test/soak.ml) decides `identical` with `( = )`, so a result containing `nan` would print `identical: false`; `Linalg.solve` and `Newton.solve` compare with exactly `0.` on purpose.

`Printf.sprintf` formats a string and the compiler checks the arguments against the directives: `%s` string, `%d` int, `%b` bool, `%g` float, plain or scientific, up to 6 significant digits, `%.3e` scientific with 3 decimals, `%.12f` fixed with 12 decimals, `%c` a `char` (one character, written `'A'`). A flag after the `%` changes the layout: `%-8g` pads on the right (`-` left-justifies) and `% .6f` leaves room for a minus sign (a space). `%!` flushes the output. `^` concatenates strings. The tests build their printed lines this way.

### Arrays, tuples and lists

```ocaml
let v = [| 1.; 2.; 3. |]
let first = v.(0)                                          (* = 1.: indices start at 0 *)
let doubled = Array.map (fun a -> 2. *. a) v               (* = [|2.; 4.; 6.|] *)
let total = Array.fold_left ( +. ) 0. v                    (* = 6.: ((0. + 1.) + 2.) + 3. *)
let squares = Array.init 4 (fun i -> float_of_int (i * i)) (* = [|0.; 1.; 4.; 9.|] *)
let m = [| [| 3.; 1. |]; [| 1.; 2. |] |]                   (* a matrix: an array of rows *)
let entry = m.(0).(1)                                      (* row 0, column 1: = 1. *)
let (lo, hi) = (v.(0), v.(2))                              (* a tuple, taken apart by a pattern *)
```

Two traps. Elements are separated by `;`: a comma builds a tuple, so `[| 1., 2. |]` is an array holding one pair. And arrays are mutable in OCaml (`a.(i) <- x`) and often shared, which is safe here only because this code never writes into an array after creating it: treat every array you receive as read-only ([architecture.md](architecture.md)). An index out of range raises `Invalid_argument`, and so does `Array.map2` on arrays of different lengths. Where: `Vec`, `Jac.forward`, `Linalg` and `Stage.solve` build arrays with `Array.init`, `Array.map` and `Array.mapi`; `Vec.finite` is `Array.for_all`; `Linalg` uses `Array.sub` and `Array.append`. A tuple carries several results at once: `M.step` returns `(y, history)`, of type `Vec.t * history`, and `let* y, history = M.step ... in` takes the pair apart. Lists (`[1; 2]`) appear only in the tests, as the case tables.

### Records

```ocaml
type point = { t : float; y : float array }
let p0 = { t = 0.; y = [| 1. |] }          (* construction *)
let p1 = { p0 with t = 1. }                (* functional update: a copy with t replaced *)
let time { t; _ } = t                      (* a pattern: binds the name t *)
let y = [| 2. |]
let p2 = { t = 2.; y }                     (* punning: { y } means { y = y } *)
```

Records here are immutable. `Halving.rejected` returns `Ok { c with dt; streak = 0; failures; stats = { c.stats with ... } }`: the same controller with some fields changed, two of them punned. Several records have a field `t` (`Ode.point`, `Stage.equation`, `Adaptive.solution`), so the compiler must decide which one a field belongs to:

```ocaml
type point = { t : float; y : float array }
type equation = { t : float; gamma : float }
let a : point = { t = 1.; y = [||] }       (* the annotation selects point *)
let b = { t = 1.; gamma = 2. }             (* the other field, gamma, selects equation *)
let time_of (p : point) = p.t              (* the annotation selects point.t *)
```

The compiler uses an annotation, the expected type or the other fields, and otherwise takes the latest definition. That is why [src/bdf1.ml](../src/bdf1.ml) writes `{ Stage.t = at.t +. h; gamma = h; psi = at.y }`: naming the module on one field names the type, and `open Ode` at the top of the file makes `at.t` mean `Ode.point`'s `t`. [test/problems.ml](../test/problems.ml) builds `{ Ode.rhs; t0 = 0.; t_end = 1.; y0 }` the same way.

### Variants, pattern matching and inline records

```ocaml
type failure = Diverged | StepRejected of int | Nan     (* Fail.t, renamed *)

let to_string = function
  | Diverged -> "Diverged"
  | StepRejected n -> Printf.sprintf "StepRejected %d" n
  | Nan -> "Nan"

let classify = function
  | Some x when x < 0. -> "negative"       (* a guard: the clause applies only if it holds *)
  | Some _ -> "non-negative"               (* _ matches anything and ignores it *)
  | None -> "missing"

let first_error a b =
  match (a, b) with
  | Ok x, Ok y -> Ok (x, y)
  | Error e, _ | _, Error e -> Error e     (* an or-pattern: either side binds e *)

type history = Start | After of { h_prev : float; y_prev : float array }
let h_prev_of = function Start -> None | After { h_prev; _ } -> Some h_prev
let next = After { h_prev = 0.5; y_prev = [| 1. |] }
```

A **variant** is a type whose values are one of several named constructors, each carrying data or not. `match` takes a value apart: the first pattern that fits wins and binds the names in it, and `function | p -> a | q -> b` is short for `fun x -> match x with | p -> a | q -> b`. The compiler checks that a match covers every case; a missing case is warning 8, an error here, so adding a constructor to `Fail.t` points at every `match` that must change. `option` (`None | Some x`), `result` (`Ok x | Error e`) and `bool` are variants too. A constructor can carry an **inline record**: named fields without a separate type, built and matched as above. `Bdf2.history` in [src/bdf2.ml](../src/bdf2.ml) has this shape, with `Vec.t` (that is, `float array`) for the vector. `Adaptive.integrate` matches a step's outcome with a guard:

```text
match M.step_with_error p.rhs h history at with
| Ok (y, err, next) when C.acceptable c ~y ~err -> go (C.accepted c) next { Ode.t = ...; y }
| Ok _ -> retry Ode.Too_large
| Error e -> retry (Ode.Solver e)
```

### `result`, `option` and `Fail.Syntax`

```ocaml
let ( let* ) = Result.bind                  (* the definition in Fail.Syntax *)
let ( let+ ) r f = Result.map f r           (* and this one *)

let parse s =
  match float_of_string_opt s with
  | Some x -> Ok x
  | None -> Error ("not a number: " ^ s)

let sum s1 s2 =
  let* x = parse s1 in                      (* Result.bind (parse s1) (fun x -> ...) *)
  let+ y = parse s2 in                      (* Result.map (fun y -> ...) (parse s2) *)
  x +. y
(* val sum : string -> string -> (float, string) result *)
```

`sum "1" "2"` is `Ok 3.` and `sum "1" "x"` is `Error "not a number: x"`: after the first `Error` the rest is skipped. `let*` continues with a function that may fail again, `let+` ends with a plain value. A type with arguments is written with the arguments first: `float option`, `(float, string) result`, `Halving.stats Adaptive.solution`. `Fail.Syntax` ([src/fail.ml](../src/fail.ml)) defines both, and `open Fail.Syntax` brings them in. Where: `Bdf2.step_with_error` (`let* be = ...`, then `let+ y = ...`) and `Stepper.fixed` (`let* y, history = M.step ...`). `Linalg.solve` returns an `option` (`None` for a zero pivot), which `Newton.solve` turns into `Error Diverged`. Numerical failures are `Fail.t` values in a `result`, never exceptions.

### Modules, `.mli` files and abstraction

```ocaml
module Canary = struct
  let lambda = [| 1.; 100.; 1e4 |]
  let y0 = [| 1.; 1.; 1. |]
end
let n = Array.length Canary.lambda            (* qualified access: = 3 *)
let m = let open Canary in Array.length y0    (* local open: Canary's names, inside this expression only *)

module Counter : sig
  type t                                      (* abstract: the representation is hidden *)
  val zero : t
  val incr : t -> t
  val to_int : t -> int
end = struct
  type t = int
  let zero = 0
  let helper n = n + 1                        (* not in the signature: invisible outside *)
  let incr n = helper n
  let to_int n = n
end
let two = Counter.(to_int (incr (incr zero))) (* = 2 *)
```

Every `.ml` file is a module named after it (`vec.ml` is `Vec`), and `module Name = struct ... end` nests one (`Problems.Canary`). `open M` exposes M's names for the rest of the file, `let open M in e` only inside `e`; `Counter.(...)` is the same for one expression, and [test/corpus.ml](../test/corpus.ml) starts its groups with `let open Problems.Canary in`. Dune wraps the library: from outside, the modules are `Vstiff.Vec`, `Vstiff.Newton` and so on, so the tests start with `open Vstiff`; inside `src/` the short names work. `Stdlib` is always open, which is why `float_of_int` and `invalid_arg` need no prefix.

Each module in `src/` has an **interface**, `foo.mli`: the list of what other modules may use. Everything else in `foo.ml` is private, and a type declared without a definition is **abstract**: outsiders can pass its values around but cannot build or inspect them. `Counter.incr 3` is rejected:

```text
Error: The constant 3 has type int but an expression was expected of type
         Counter.t
```

A name missing from the `.mli` does not exist outside (`Vstiff.Newton.tol` is `Unbound value`). That is how Newton's tuning constants and `Linalg.pivot` stay internal, how `Halving.t` stays opaque, and how `Bdf2.history` hides `Start` and `After`. The documentation comments in the `.mli` files are the API documentation.

### Module types, `include` and interface-only modules

```ocaml
module type Policy = sig
  type t
  type stats
  val init : t
  val stats : t -> stats
end

module type Counting = sig
  type stats = { accepted : int }
  include Policy with type stats := stats
end

module Count : Counting = struct
  type stats = { accepted : int }
  type t = int
  let init = 0
  let stats n = { accepted = n }
end

let n = (Count.stats Count.init).accepted
```

A **module type** lists the types and values a module must provide; a module satisfies it when it provides at least those, with compatible types. There is no `implements` keyword. `include S` copies the items of `S` into another signature: `Ode.Embedded` is `Ode.Method` plus `step_with_error`. `with type stats := stats` substitutes: it removes `Policy`'s own `stats` and uses ours, without which the signature would define `stats` twice (`Multiple definition of the type name stats`). `halving.mli` does exactly this with `Ode.Controller`, so `Halving.stats` is a record whose fields (`accepted_steps`) the tests can read. `bdf1.mli` is `include Ode.Method` and `bdf2.mli` is `include Ode.Embedded` plus `coeffs`.

[src/ode.mli](../src/ode.mli) has no `ode.ml`: it holds types and module types, so there is nothing to run. Dune has to be told, by `(modules_without_implementation ode)` in `src/dune`; without it the build stops with `Some modules don't have an implementation`.

### First-class modules and modular explicits

```ocaml
module type Shape = sig
  type t
  val make : float -> t
  val area : t -> float
end

module Square = struct
  type t = float
  let make side = side
  let area s = s *. s
end

(* modular explicit: a module parameter whose types the result type may name *)
let make_twice (module S : Shape) size : S.t * S.t = (S.make size, S.make size)
(* val make_twice : (module S : Shape) -> float -> S.t * S.t *)
let (s1, s2) = make_twice (module Square) 3.      (* both have type Square.t *)

(* first-class module: an ordinary value of type (module Shape) *)
let packed = (module Square : Shape)
let area_of (m : (module Shape)) size = let module S = (val m) in S.area (S.make size)
let nine = area_of packed 3.                      (* = 9.; S.t may not escape area_of *)
```

A **first-class module** is a module packed into a value, so it can sit in a list or a record. Inside the function that unpacks it, its types are abstract: returning an `S.t` is rejected (`The type constructor S.t would escape its scope`), unless the type of the packed module says what `t` is, as in `(module Shape with type t = a)`. A **modular explicit** (OCaml 5.5) is a module parameter written `(module S : Shape)`; the result type may name `S.t`, and the caller receives the types of the module it passed. The argument must be a static module, `(module Square)`; a packed value is refused:

```text
Error: This expression has type (module S : Shape) -> float -> S.t * S.t
       but an expression was expected of type (module Shape) -> 'a
       The module S would escape its scope
This function is module-dependent. The dependency is preserved
when the function is passed a static module argument (module M : S)
or (module M). Its argument here is not static, so the type-checker
tried instead to change the function type to be non-dependent.
```

Where: [src/adaptive.mli](../src/adaptive.mli) declares

```text
val integrate :
  (module M : Ode.Embedded) -> (module C : Ode.Controller) ->
  ?dt0:float -> ?dt_max:float -> ?max_rejects:int -> tol:float -> Ode.problem ->
  (C.stats solution, Fail.t) result
```

The result type names `C.stats`: `Adaptive.integrate (module Bdf2) (module Halving) ~tol problem` returns `(Halving.stats Adaptive.solution, Fail.t) result`, so the tests read `s.stats.rejected_steps`, and another controller with another stats type gives another result type. `repeat (module C : Case) : (C.t, Fail.t) result list` in [test/soak.ml](../test/soak.ml) has the same shape. `Stepper.fixed (module M : Ode.Method)` takes a module the same way.

**Why this code prefers modular explicits** to functors and to first-class modules packed into values. The result type names the module's own types (`C.stats`), without the `with type` plumbing that a packed module needs. A functor (a function from modules to a module) would turn every combination into a named module (`module Run = Adaptive.Make (Bdf2) (Halving)`) before it could be called; here the call site just writes `(module Bdf2) (module Halving)`, and partial application works. The cost is that the argument must be written out in the source, so the choice of method cannot come from a runtime list; a first-class module is the tool for that.

### Exceptions and mutable cells

```ocaml
exception Exhausted

let positive dt = if dt > 0. then dt else invalid_arg "dt must be positive"
let outcome =
  match positive (-1.) with
  | v -> Ok v
  | exception Invalid_argument msg -> Error msg      (* = Error "dt must be positive" *)

let counter () =
  let calls = ref 0 in                               (* a mutable cell *)
  ((fun () -> incr calls), fun () -> !calls)         (* incr adds one; ! reads the cell *)
```

`invalid_arg msg` raises `Invalid_argument msg`; a `match` case written `| exception ...` catches an exception raised while the matched expression is evaluated. These tools are quarantined. Invalid arguments are programming errors, and in the library only `Check` ([src/check.ml](../src/check.ml)) raises on purpose; in the tests only [test/guard.ml](../test/guard.ml) raises (`Exhausted`, when a call budget runs out) or catches. The only mutable state in the library is the `ref` in `Instrument.count`, whose shape `counter` copies. A `ref` is read with `!` (not a negation) and written with `:=`. [architecture.md](architecture.md) gives the reasons.

### Tail recursion

```ocaml
let rec sum_to n = if n = 0 then 0 else n + sum_to (n - 1)   (* the + waits for the result: not a tail call *)

let sum_to_acc n =
  let rec go acc k = if k = 0 then acc else go (acc + k) (k - 1) in
  go 0 n                                                      (* go is the last thing go does: a tail call *)
```

This code has no `for` or `while` loops (OCaml has them): repetition is recursion. A call in tail position reuses the stack frame, so a loop written that way runs in constant stack however long it is. `go` in `Stepper.fixed` and `Adaptive.integrate`, and `iterate` and `damp` in `Newton.solve`, are tail recursive (`let*` is `Result.bind`, which calls its function last), so the number of steps is not limited by the stack. The recursion in `Linalg` is not a tail call; its depth is the number of unknowns.

### Warnings are errors

In the dev profile most compiler warnings stop the build. The ones you will meet:

- **26, 27: unused variable.** `Error (warning 27 [unused-var-strict]): unused variable y.` Remove it, or start its name with an underscore (`_t` in `let rhs _t y = ...`, [test/problems.ml](../test/problems.ml)).
- **8: a `match` that misses a case**: `Error (warning 8 [partial-match]): this pattern-matching is not exhaustive.` followed by an example of the missing case.
- **32: unused top-level value**, in a module that has an `.mli` and in the main module of a test or probe program (dune gives those an empty interface).
- **33, 39, 11, 16**: unused `open`, unused `rec`, an unused match case, an optional argument that can never be defaulted.
- **50: a documentation comment in the wrong place**: `Error (warning 50 [unexpected-docstring]): unattached documentation comment (ignored)`.

Fix the code; never loosen the flags. Three habits break the build when you write comments: comments nest, so the two characters `(*` inside a comment open another one (write `( *. )` with spaces); string literals are lexed inside comments, so a double quote must belong to a closed string (avoid it; the same goes for the two characters `{|`); and a `(** ... *)` comment must sit directly before or after the item it documents (a definition, a constructor or a record field), while inside a function body only `(* ... *)` is allowed. Documentation comments use the odoc markup: `[code]`, `{[ ... ]}` for a block, `{b bold}`, `{!Module.name}` for a link and `-` for list items. Dune files use `;` for comments.
