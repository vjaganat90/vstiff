# verif: proofs about vstiff in Rocq and MathComp

This directory holds machine-checked theorems about a Rocq mirror of vstiff's numerics, proved with
the Mathematical Components libraries. The main build ignores it (the root `dune` file says
`data_only_dirs verif`), so the library and its tests never need Rocq.

## The proof switch

The proofs build in their own opam switch, on the same OCaml as the project. MathComp's
dependency `rocq-elpi` needs dune below 3.24, so the development switch cannot host it.

```sh
opam switch create vstiff-proofs ocaml-base-compiler.5.5.0 --no-switch \
  --repos=default,rocq-released=https://rocq-prover.org/opam/released
opam install --switch=vstiff-proofs dune.3.23.1 rocq-core.9.2.0 rocq-stdlib.9.2.0 \
  rocq-mathcomp-ssreflect.2.6.0 rocq-mathcomp-algebra.2.6.0
```

That gives OCaml 5.5.0, dune 3.23.1, Rocq 9.2.0 and MathComp 2.6.0. Move to Rocq 9.3 once
MathComp lists it.

## Building

From the repository root:

```sh
opam exec --switch=vstiff-proofs -- dune build --root verif
```

A green build means every proof in `theories/` checks. `theories/dune` lists the whole theory
closure MathComp 2.6 needs (the standard library, micromega, Elpi, Hierarchy Builder, then the
MathComp parts); MathComp 2.6 also renamed `all_ssreflect` to `boot` and `order`.

## What is proved

| File | Statement | Mirrors |
|---|---|---|
| [theories/Bdf2_weights.v](theories/Bdf2_weights.v) | `weights_sum_to_one`: over any field, `a1 w + a0 w = 1` whenever `1 + 2 w <> 0` | `Bdf2.coeffs` in [src/bdf2.ml](../src/bdf2.ml) |

Next come the proofs that the formula is exact on quadratics and that variable-step BDF2 is
zero-stable.
