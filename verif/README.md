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
polynomial $p$, and `q w`, the parasitic root $q(\omega) = -a_0(\omega)$.
[`Clock.resolution`](../src/clock.ml) has no mirror yet; the tripwire watches it already, for a future proof about the controller.

## What is proved

$\omega = h / h_{\mathrm{prev}}$ is the step ratio, $h_{\mathrm{prev}}$ the previous step and
$h = \omega\thinspace h_{\mathrm{prev}}$ the new one. $a_1(\omega)$, $a_0(\omega)$ and $\beta(\omega)$ are the
coefficients that `Bdf2.coeffs` returns and the mirror (`a1`, `a0`, `beta`) transcribes;
$\operatorname{den}(\omega)$ (`den`) is their shared denominator.

**T1**, in [theories/Bdf2_weights.v](theories/Bdf2_weights.v), over any field, whenever
$\operatorname{den}(\omega) \ne 0$:

- `weights_sum_to_one`: $a_1(\omega) + a_0(\omega) = 1$.
- `bdf2_exact_on_quadratics`: every polynomial $p$ of degree at most 2 satisfies the step
  exactly, $p(t_n + h) = a_1(\omega)\thinspace p(t_n) + a_0(\omega)\thinspace p(t_n - h_{\mathrm{prev}}) + \beta(\omega)\thinspace h\thinspace p'(t_n + h)$.
- `bdf2_residual_on_cubics`: on a polynomial of degree at most 3, the residual (left side minus
  right side) is $-p_3\thinspace \beta(\omega)\thinspace h^2\thinspace (h + h_{\mathrm{prev}})$, where $p_3$ is its coefficient of $t^3$.
- `bdf2_cubic_residual`: the residual on $t^3$ is $-\beta(\omega)\thinspace h^2\thinspace (h + h_{\mathrm{prev}})$.
- `bdf2_cubic_residual_shifted`: $a_1(\omega)\thinspace h^3 + a_0(\omega)\thinspace (h + h_{\mathrm{prev}})^3 = -\beta(\omega)\thinspace h^2\thinspace (h + h_{\mathrm{prev}})$,
  the residual on $(t - (t_n + h))^3$, the cubic with its origin at the new time point.

**T2**, in [theories/Bdf2_zero_stability.v](theories/Bdf2_zero_stability.v), for the recurrence
$y_{n+2} = a_1(\omega_n)\thinspace y_{n+1} + a_0(\omega_n)\thinspace y_n + d_n$ with ratios $\omega_n$ and perturbations
$d_n$ (rounding, the stage residual, or the $f$-term); $d = 0$ is the homogeneous recurrence.
$q(\omega) = -a_0(\omega)$ is the factor a step applies to the step difference:

- `qE`, over any field: $q$ in the closed form the other lemmas use.
- `bdf2_increment`, over any field, when $\operatorname{den}(\omega) \ne 0$:
  $a_1(\omega)\thinspace y_1 + a_0(\omega)\thinspace y_0 - y_1 = q(\omega)\thinspace (y_1 - y_0)$.
- `q_le`, over an ordered field: $0 \le \omega \le \omega_{\mathrm{s}}$ implies $q(\omega) \le q(\omega_{\mathrm{s}})$.
- `bdf2_contraction`, over an ordered field: if every $\omega_n$ lies in $[0, \omega_{\mathrm{s}}]$, the
  homogeneous recurrence has $\lvert y_{n+1} - y_n \rvert \le q(\omega_{\mathrm{s}})^n\thinspace \lvert y_1 - y_0 \rvert$.
- `bdf2_bounded`, over an ordered field: if moreover $q(\omega_{\mathrm{s}}) \lt  1$, it has
  $\lvert y_n \rvert \le \lvert y_0 \rvert + \lvert y_1 - y_0 \rvert / (1 - q(\omega_{\mathrm{s}}))$.
- `bdf2_stable`, over an ordered field: under the same hypotheses the perturbed recurrence has
  $\lvert y_n \rvert \le \lvert y_0 \rvert + \bigl(\lvert y_1 - y_0 \rvert + \sum_{k\lt n} \lvert d_k \rvert\bigr) / (1 - q(\omega_{\mathrm{s}}))$.
- `q_lt1`, over a real closed field: for $\omega \ge 0$, $q(\omega) \lt  1$ exactly when $\omega \lt  1 + \sqrt{2}$.
- `q_threshold`, over a real closed field: $q(1 + \sqrt{2}) = 1$.
- `bdf2_zero_stable`, over a real closed field: if $\omega_{\mathrm{s}} \lt  1 + \sqrt{2}$ and every $\omega_n$ lies
  in $[0, \omega_{\mathrm{s}}]$, then $q(\omega_{\mathrm{s}}) \lt  1$ and the bound of `bdf2_stable` holds.
- `bdf2_sharp`, over a real closed field: for a constant ratio $\omega \ge 1 + \sqrt{2}$,
  $q(\omega) \ge 1$ and every homogeneous solution has $\lvert y_n - y_0 \rvert \ge n\thinspace \lvert y_1 - y_0 \rvert$,
  unbounded unless $y_1 = y_0$.

The bounds hold componentwise, hence in the max norm. Zero-stability is the $h \to 0$ notion:
it says nothing about stiff stability at a step ratio of $2$.

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

`Halving` at most doubles a step, and the driver's snapping moves a step by at most $1/32$ of its
length, so the step ratios of a run stay at or below
$2 \cdot \dfrac{33/32}{31/32} = \dfrac{66}{31}$, inside the bound of `q_lt1`. A proof of
that, on a mirror of `Halving` and the driver, comes next. No CI job builds the proofs yet.
