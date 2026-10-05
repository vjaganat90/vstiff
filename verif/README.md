# verif: proofs about vstiff in Rocq and MathComp

This directory holds machine-checked theorems about a Rocq mirror of vstiff's numerics, proved with
the Mathematical Components libraries. The main build ignores it (the root `dune` file says
`data_only_dirs verif`), so the library and its tests never need Rocq. The tripwire that ties the
mirror to the code is in the main build, in [sync/](../sync).

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

A green build means every proof in `theories/` checks. None uses `Admitted`, `admit` or an axiom:
`Print Assumptions` prints `Closed under the global context` for each theorem below.
`theories/dune` lists the whole theory closure MathComp 2.6 needs (the standard library,
micromega, Elpi, Hierarchy Builder, then the MathComp parts); MathComp 2.6 also renamed
`all_ssreflect` to `boot` and `order`. Its `ring`, `field`, `lra` and `nra` tactics come with
`algebra`.

## The mirror

Each mirror keeps the operation order of the OCaml expression it transcribes, and a comment names
the OCaml it mirrors.

| Rocq, in [theories/Bdf2_weights.v](theories/Bdf2_weights.v) | Mirrors, in [src/bdf2.ml](../src/bdf2.ml) |
|---|---|
| `den w` | the local `d` of `Bdf2.coeffs omega` |
| `a1 w`, `a0 w`, `beta w` | the fields `a1`, `a0` and `beta` of `Bdf2.coeffs omega` |

The other definitions are built on the mirror: `residual w hp tn p`, the residual of one step on a
polynomial `p`, and `q w = - a0 w`, the parasitic root. [`Clock.resolution`](../src/clock.ml) has no
mirror yet; the tripwire watches it already, for a future proof about the controller.

## What is proved

`w` is the step ratio `h / h_prev`, `hp` the previous step and `h = w hp` the new one.

**T1**, in [theories/Bdf2_weights.v](theories/Bdf2_weights.v), over any field, whenever
`den w != 0`:

- `weights_sum_to_one`: `a1 w + a0 w = 1`.
- `bdf2_exact_on_quadratics`: every polynomial `p` of degree at most 2 satisfies the step exactly,
  `p(tn + h) = a1 w p(tn) + a0 w p(tn - hp) + beta w h p'(tn + h)`.
- `bdf2_residual_on_cubics`: on a polynomial of degree at most 3, the residual (left side minus
  right side) is `- p3 beta w h^2 (h + hp)`, where `p3` is its coefficient of `t^3`.
- `bdf2_cubic_residual`: the residual on `t^3` is `- beta w h^2 (h + hp)`.
- `bdf2_cubic_residual_shifted`: `a1 w h^3 + a0 w (h + hp)^3 = - beta w h^2 (h + hp)`, the
  residual on $`(t - (t_n + h))^3`$, the cubic with its origin at the new time point.

**T2**, in [theories/Bdf2_zero_stability.v](theories/Bdf2_zero_stability.v), for the recurrence
`y (n+2) = a1 (w n) y (n+1) + a0 (w n) y n + d n` with ratios `w n` and perturbations `d n`
(rounding, the stage residual, or the f-term); `d = 0` is the homogeneous recurrence:

- `qE`, over any field: `q w = w^2 / (1 + 2 w)`.
- `bdf2_increment`, over any field, when `den w != 0`: `a1 w y1 + a0 w y0 - y1 = q w (y1 - y0)`,
  so a step multiplies the step difference by `q w`.
- `q_le`, over an ordered field: `0 <= w <= ws` implies `q w <= q ws`.
- `bdf2_contraction`, over an ordered field: if every `w n` lies in `[0, ws]`, the homogeneous
  recurrence has `|y (n+1) - y n| <= q ws ^ n |y 1 - y 0|`.
- `bdf2_bounded`, over an ordered field: if moreover `q ws < 1`, it has
  `|y n| <= |y 0| + |y 1 - y 0| / (1 - q ws)`.
- `bdf2_stable`, over an ordered field: under the same hypotheses the perturbed recurrence has
  `|y n| <= |y 0| + (|y 1 - y 0| + sum over k < n of |d k|) / (1 - q ws)`.
- `q_lt1`, over a real closed field: for `w >= 0`, `q w < 1` exactly when `w < 1 + sqrt 2`.
- `q_threshold`, over a real closed field: `q (1 + sqrt 2) = 1`.
- `bdf2_zero_stable`, over a real closed field: if `ws < 1 + sqrt 2` and every `w n` lies in
  `[0, ws]`, then `q ws < 1` and the bound of `bdf2_stable` holds.
- `bdf2_sharp`, over a real closed field: for a constant ratio `w >= 1 + sqrt 2`, `q w >= 1` and
  every homogeneous solution has `|y n - y 0| >= n |y 1 - y 0|`, unbounded unless `y 1 = y 0`.

The bounds hold componentwise, hence in the max norm. Zero-stability is the $`h \to 0`$ notion:
it says nothing about stiff stability at a step ratio of 2.

## The tripwire

The theorems are about the mirror, not the compiled OCaml. The tripwire stops the code from
drifting away from the mirror unnoticed. It is in the main build because a separate dune root
cannot read `../src`:

- [sync/canon.ml](../sync/canon.ml), an executable on compiler-libs: `canon.exe file.ml` prints the
  canonical form of `file.ml`, its parse tree without comments or doc comments, laid out by the
  compiler's own printer, so that layout does not count either.
- [sync/dune](../sync/dune): `dune runtest` diffs the canonical forms of
  [src/bdf2.ml](../src/bdf2.ml) and [src/clock.ml](../src/clock.ml) with the committed
  [sync/bdf2.ml.canon](../sync/bdf2.ml.canon) and [sync/clock.ml.canon](../sync/clock.ml.canon).

A code change anywhere in either module, even a renamed local, fails `dune runtest` with a diff
of the canonical form; an edit to comments or layout passes. When it fails:

1. Re-check the mirror: compare the diff with the Rocq definitions of the table above, update them
   and the proofs if the mirrored code moved, and rebuild the proofs.
2. Promote that file alone, for example `dune promote sync/bdf2.ml.canon`: a bare `dune promote`
   would also accept a pending diff of the expect files.
3. Commit the new `.canon` with the change to the mirror, if there is one.

The tool reads `.ml` files only, so a contract change in an `.mli` shows only when the `.ml`
changes too. The compiler's printer can change its layout between releases, so an OCaml upgrade
can make every `.canon` differ: promote them once the diffs show layout only.

## Next

A proof that the controller never proposes a step ratio at or above the bound of `q_lt1` waits for
a controller that enforces one. No CI job builds the proofs yet.
