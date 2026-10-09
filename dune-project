(lang dune 3.0)
; The dune language version, not the compiler's. It must be the first line: dune rejects a comment before it.

(name vstiff)
; The project name. The library's own name is set in src/dune.

; Nothing here pins the OCaml version: the project is built and tested with OCaml 5.5.0 and dune 3.x, and
; older OCaml releases reject the modular explicits in the sources. docs/ocaml.md has the setup.
