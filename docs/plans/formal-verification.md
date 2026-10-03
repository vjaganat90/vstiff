# Formal verification of vstiff with Rocq and MathComp

Landscape, bridge to OCaml, ranked candidate theorems, and a pilot plan. Companion to the [roadmap](roadmap.md).

- Date: 2026-10-03. Code read at `f6d9b4e`. Two fixes landed afterwards, `052fdc7` (`Newton.solve` reports a converged step that overflows as `Error Nan`) and `dc3bb67` (`Adaptive.integrate` snaps steps to the clock and rejects steps that cannot move `t`; the floor `16 eps |t|` is now held by the new module `Clock`), and the text states the current behaviour where they matter (findings 1 and 2, T3, T4, T10). The comment-only rewrite committed in between (`10ff44d`) leaves the comment-free canonical sources of all 19 `.ml` files in `src/` and `test/` identical to `f6d9b4e` (section 2.1), and the dune files gained only comments.
- No Rocq toolchain was installed. **Every Rocq, Gappa and dune snippet below is an unchecked sketch.** Library versions and the solver results in section 1.2 were checked: the solver results come from opam `--dry-run` runs in a throwaway opam root.
- Every OCaml claim about runtime behaviour was checked by running the committed code in throwaway programs (Appendix A).
- The main solver rows of section 1.2 were run twice (dry runs, nothing installed) with the same results, and the opam constraints behind the other rows were re-read.

---

## 0. Summary

**What is realistic.** Real-arithmetic theorems about the method and the step-control logic can be proved with MathComp. The proofs would run against a hand-written Rocq "mirror" of a few OCaml functions. Tests and a drift tripwire tie the mirror to the code. A few binary64 theorems are also within reach: special-value safety, and rounding bounds on small closed formulas. Proving the whole adaptive integrator correct in binary64 is not realistic.

**What OCaml does and does not buy.** Rocq is written in OCaml and extracts to OCaml, so the Rocq-to-OCaml direction is mature. The OCaml-to-Rocq direction is weak for numerical code:

- `rocq-of-ocaml` models `float` as `Z`, and turns float literals into integers.
- CFML needs OCaml < 5.
- None of the translators can read OCaml 5.5 source such as modular explicits; `rocq-of-ocaml` is limited to OCaml 5.4 by its Rocq 9.0 dependency.

Extraction has its own catch. Rocq's primitive floats have no fused multiply-add, but `ocamlopt` fuses five vstiff expressions on arm64. Extracted kernels would therefore not be bit-identical to the current build.

**Toolchain fact that shapes everything.** MathComp cannot go into the current development switch:

- dune 3.24 removed the old `(using coq ...)` extension.
- `rocq-elpi`, which Hierarchy Builder and therefore MathComp 2.6 depend on, is bounded to `dune < 3.24`.
- The development switch has dune 3.24.2.

The solver accepts this separate switch: **dune 3.23.1 + Rocq 9.2.0 + MathComp 2.6.0** on OCaml 5.5.0, and also on OCaml 4.14.2, the lower-risk choice (4.3). It also accepts adding MathComp-Analysis 1.18.0, Flocq 4.2.2, Coquelicot 3.4.5 and CoqInterval 4.11.5, and even VCFloat 2.4.2 and LAProof 2.0.1.

**Findings from stating theorems precisely.** All five were confirmed by running the code at `f6d9b4e`; the first two have since been fixed:

1. **[`Adaptive.integrate`](../../src/adaptive.ml) could loop forever, and could return a wrong state as `Ok`.** Fixed in `dc3bb67`. At `f6d9b4e` the loop did not end when every proposed step was below half an ulp of `t`. Example: `t0 = 1.`, `t_end = Float.succ 1.`, all options at their defaults, `y' = -y`. After 3e6 right-hand-side calls, the largest `t` ever passed to `rhs` was still exactly `1.0`. The same happened with `t0 = 1e10`, `dt_max = 5e-7`. The 16 eps |t| floor applied only on rejection, and here every step was accepted. A quieter variant: when `h` was only a few ulps of `t`, the state advanced by `h` but the clock by `RN(t + h) - t`, and the run ended `Ok` at `t = t_end` with a wrong state. For `y' = 1` from `t0 = 1e15` over a span of 100 (ulp(t) = 0.125), `dt_max = 0.19` returned `y = 76.36` instead of 100, and the defaults (`dt_max = 10`) returned `y = 100.32`. Since the fix every non-final step is snapped to the floats and a step that cannot move `t` is rejected as `Too_small` without calling the method (T4); corpus lines 25 to 27 pin the three cases.
2. **[`Newton.solve`](../../src/newton.ml) could return `Ok [|infinity|]`.** Fixed in `052fdc7`. The converged result `x + dx` was never checked for finiteness. Example: `f _ = [|-1e298|]`, Jacobian `[[1]]`, `x0 = [|max_float|]`. A converged step that overflows is now reported as `Error Nan`; corpus line 24 pins it (it printed `Ok [inf]` before the fix).
3. **[`Stepper.fixed`](../../src/stepper.ml) does not land exactly on `t_end`.** The `.mli` at `f6d9b4e` said the last step "lands exactly on `t_end`", which holds only in real arithmetic; the current `.mli` no longer says so, and the code is unchanged. With `t_end = 1`, `dt = 0.1`, the last stage equation is solved at `t = 0.99999999999999989`.
4. **BDF2 never sees a step ratio above 2 under [`Halving`](../../src/halving.ml) in real arithmetic.** Instrumented runs at `f6d9b4e` on the three adaptive corpus problems reached exactly ω = 2: 358, 1362 and 842 times, and never above. That is safely inside the sharp zero-stability bound 1 + √2 ≈ 2.414. It is above the classical energy-stability bounds for stiff dissipative problems (1.87 to 1.94). Only newer results cover ω = 2 there (section 1.6). Since `dc3bb67` non-final steps are snapped, so a ratio can exceed 2 by rounding (T3).
5. **`ocamlopt` 5.5.0 on arm64 emits `fmadd`/`fmsub` for five expressions.** They are in [`Vec.axpy`](../../src/vec.ml), [`Jac.forward`](../../src/jac.ml), [`Stage.solve`](../../src/stage.ml), the Newton Armijo threshold and [`Bdf2.coeffs`](../../src/bdf2.ml). Blocking all five with `Sys.opaque_identity` changed nothing:
   - no corpus or soak line;
   - no step count;
   - no bit of the final Robertson, logistic and van der Pol states.

**Recommended pilot.** Three theorems, all in real arithmetic, needing only MathComp:

- **T1:** the variable-step BDF2 coefficients are exact on quadratics.
- **T2:** variable-step BDF2 is zero-stable, uniformly in the grid, when every step ratio is ≤ ω* < 1 + √2.
- **T3:** `Halving`, `Adaptive` and `Bdf2` together never attempt ω > 2 in real arithmetic, except for a final step within the clock's resolution. This is a cross-module invariant.

Together they say every BDF2 run vstiff can produce is zero-stable, with a step-difference contraction factor q(2) below 1 (the parasitic root at the largest ratio, from `Bdf2.coeffs`; T2). Zero-stability is the h → 0 property of the homogeneous recurrence; it says nothing about stiff stability at ratio 2 (that is T12). Estimate: about 1 to 2 weeks for someone fluent in MathComp, or 3 to 6 weeks for a newcomer who already writes Rocq proofs, including a one-time ramp-up; add more if Rocq itself is new (OCaml fluency helps with the mirror, not with the proofs). The pilot's one code prerequisite, a fix for finding 1 landed together with its pins, is in place (`dc3bb67`, corpus lines 25 to 27) and discharges T3's hypothesis h_prev > 0 (section 4.2).

---

## 1. Landscape (checked 2026-10-03)

### 1.1 Rocq versions and OCaml compatibility

Rocq is the Coq prover renamed from version 9.0 on. The opam packages are `rocq-runtime`, `rocq-core`, `rocq-stdlib`, the meta package `rocq-prover`, and the compatibility packages `coq`/`coq-core`/`coq-stdlib`.

| Rocq | Released | OCaml bound (opam) | dune bound (opam) | Notes |
|---|---|---|---|---|
| 9.0.0 / 9.0.1 | 2025-03-12 / 2025-09-27 | `< 5.5` | `< 3.24` | Not installable with OCaml 5.5 |
| 9.1.0 | 2025-09-15 | `>= 4.14 & < 5.5` | `< 3.24` | Not installable with OCaml 5.5 |
| 9.1.1 | 2026-02-09 | `>= 4.14` | `< 3.24` | OCaml 5.5 allowed |
| 9.2.0 | 2026-03-27 | `>= 4.14` | `>= 3.8 & < 3.24` | Recommended (widest library support) |
| 9.3.0 | 2026-09-19 | `>= 4.14` | `>= 3.21` | Newest; VCFloat, LAProof and taylor_rocqs not yet compatible |

- opam-repository bounded `rocq-runtime`/`coq-core` versions from 8.19.2 up to (not including) 9.1.1 below OCaml 5.5 on 2026-03-07 (commit "rocq-runtime/coq-core >= 8.19.2 < 9.1.1 is not compatible with OCaml 5.5"). Later releases have no upper OCaml bound.
- Rocq's own `INSTALL.md` at V9.3.0 still says "Support for OCaml 5.x remains experimental". On master it says "tested up to OCaml 5.4.1".
- Rocq 9.3 "Restored support for `native_compute` with OCaml 5 on ARM64".
- The Rocq opam documentation warns that OCaml 5 with a system ocamlfind breaks installation, and recommends an independent switch.

**Can Rocq sit next to OCaml 5.5?** Yes for Rocq ≥ 9.1.1 as far as opam metadata goes, and the solver agrees (1.2). Metadata is not testing: upstream calls OCaml 5.x experimental, lists no test of 5.5, and on arm64 `native_compute` with OCaml 5 exists only from Rocq 9.3 (the 9.2 changelog: "only some x86 setups are supported"). With MathComp it needs dune ≤ 3.23.1, so not in the current development switch, which has dune 3.24.2 and should be left alone anyway.

### 1.2 Switch layouts tried with the opam solver (dry-run only)

Setup: opam 2.5.2, repositories `default` (opam.ocaml.org) and `rocq-released` (https://rocq-prover.org/opam/released), both fetched 2026-10-03. Each layout was solved with `opam switch create probe --dry-run` in a throwaway `OPAMROOT`. Nothing was installed.

| Layout | Solver result |
|---|---|
| OCaml 5.5.0 + **dune 3.24.2** + any Rocq + MathComp 2.6.0 | **No solution**: "rocq-hierarchy-builder → … → dune < 3.24" |
| OCaml 5.5.0 + dune 3.24.2 + Rocq 9.3.0 + Stdlib + Flocq 4.2.2 (no MathComp) | Solution, 22 packages |
| **OCaml 5.5.0 + dune 3.23.1 + Rocq 9.2.0 + Stdlib 9.2.0 + MathComp 2.6.0 + MathComp-Analysis 1.18.0 + Flocq 4.2.2 + Coquelicot 3.4.5 + CoqInterval 4.11.5** | **Solution, 72 packages** (includes `rocq-elpi 3.5.1`, `rocq-hierarchy-builder 1.10.3`) |
| The same with Rocq 9.2.0 + VCFloat 2.4.2 + LAProof 2.0.1 + gappalib-coq 1.11.0 | Solution, 86 packages (pulls in `coq-compcert 3.18` and `coq-vst 2.17`) |
| OCaml 5.4.1 + dune 3.23.1 + Rocq 9.2.0 + MathComp 2.6.0 + Flocq | Solution, 58 packages |
| OCaml 5.5.0 + dune 3.23.1 + Rocq 9.3.0 + MathComp 2.6.0 + Analysis 1.18.0 + Flocq | Solution, 66 packages |
| Rocq 9.3.0 + VCFloat (any version) | **No solution**: VCFloat needs Rocq < 9.3 |

A second run (throwaway opam root, `opam switch create --dry-run`, nothing installed): the 72-package, 66-package and 22-package solutions and the dune 3.24.2 + MathComp no-solution rows reproduce; the VCFloat row reproduces. The same sets also solve with OCaml 4.14.2 (pilot set 51 packages, full set with VCFloat and LAProof 86), and `coq-mathcomp-algebra-tactics` together with MathComp 2.6.0 has no solution. Package counts shift by a few with the exact set and OCaml version. No `rocq-stdlib` 9.3.0 exists yet, so the Rocq 9.3.0 rows use `rocq-stdlib` 9.2.0, which opam allows up to `rocq-core < 9.4~`.

### 1.3 dune's Rocq support

dune's Rocq support changed a lot in 2026 ([dune CHANGES](https://github.com/ocaml/dune/blob/main/CHANGES.md), [dune Rocq docs](https://dune.readthedocs.io/en/latest/rocq.html)):

| dune | Date | Rocq build language |
|---|---|---|
| 3.21.0 | 2026-01-12 | New `(using rocq 0.11)` build language with `rocq.theory`. `(lang coq)` deprecated. |
| 3.22.0 | 2026-03-18 | Rocq lang 0.12: expected-output tests (a `.v` file with a `.expected` file, diffed under `runtest`). |
| 3.23.0 | 2026-05-04 | Rocq lang 0.13: `rocq.extraction` takes `extracted_files`. |
| 3.24.0 | 2026-06-21 | Removes `(using coq ...)` entirely. Rocq lang 0.14 replaces `(stdlib no)` with `(no_corelib)`. |

Rocq 9.3 itself requires dune ≥ 3.21 for `rocq.theory`.

Because MathComp's dependencies still build through the removed extension, the `rocq-released` repository (github.com/rocq-prover/opam, not opam-repository) bounded `rocq-elpi` below dune 3.24 on 2026-08-06 and marked it "incompatible with dune 3.24" on 2026-09-16 (`rocq-elpi` 3.5.1 has `dune < 3.24`). The default opam-repository had already bounded `rocq-runtime` 9.0 to 9.2 below dune 3.24 on 2026-06-03, so only Rocq 9.3 can live in a dune 3.24 switch. The elpi constraint is temporary and will probably lift. Until then a MathComp project uses dune 3.23.x and `(using rocq 0.13)`.

### 1.4 Libraries

| Library | Latest (opam) | Date | Rocq | What it is, and what it is for in vstiff |
|---|---|---|---|---|
| [MathComp](https://github.com/math-comp/math-comp) (`rocq-mathcomp-*`) | 2.6.0 | 2026-07-10 | 9.0–9.2 per its changelog; opam also allows 9.3, and a `2.6.0-rocq-prover-9.3` Docker image exists | Algebra hierarchy, `'M[R]_(m,n)` matrices, `{poly R}`, ordered fields. 2.6 ships `ring`, `field`, `lra`, `nra` and `psatz` for its own structures (files `ring_tactic.v`, `field_tactic.v`, `arithmetic_tactic.v`; the old `ring.v` and `lra.v` are deprecated shims, the standalone `coq-mathcomp-algebra-tactics` 1.2.7 requires MathComp < 2.6, and `rocq-micromega-plugin` is a new dependency). `field` closes denominator side conditions with `done` (`field by tac`, or `field?` to inspect them). Logical paths: `mathcomp.boot`, `.order`, `.ssreflect`, `.algebra`, and `mathcomp.finite_group` (formerly `fingroup`). `algebra/matrix.v` contains `cormen_lup` with `cormen_lup_correct : P * A = L * U`, a recursive LUP with the same shape as [`Linalg.solve`](../../src/linalg.ml). **T1–T3, T6a, T9.** |
| [MathComp-Analysis](https://github.com/math-comp/analysis) | 1.18.0 | 2026-09-02 | 9.0–9.3 | Real analysis over MathComp: `derive`, `MVT`, `Rolle`, Banach contraction (`contraction_cvg_fixed`), Lebesgue integral, FTC, normed modules over matrices. **T7, T8, T11.** Multivariate calculus is the thin part. |
| [Flocq](https://flocq.gitlabpages.inria.fr/) | 4.2.2 | 2026-02-19 | ≥ 8.17 (any 9.x) | IEEE-754 formalization: `Binary.v` and `BinarySingleNaN.v` with `Bplus_correct`, `Bmult_correct`, `Bdiv_correct` and `Bfma_correct`. `IEEE754/PrimFloat.v` relates Rocq primitive floats to Flocq binary64 (`add_equiv`, `mul_equiv`, `next_up_equiv`, …). **All binary64 work.** |
| [Coquelicot](https://coquelicot.gitlabpages.inria.fr/) | 3.4.5 | 2026-08-27 | any | Real analysis on Stdlib's `R`, including Taylor–Lagrange. Alternative to MathComp-Analysis for 1-D calculus (T7, T11). |
| [CoqInterval](https://coqinterval.gitlabpages.inria.fr/) | 4.11.5 | 2026-09-16 | any 9.x | `interval` tactic for bounds on real expressions. Concrete constants (values of q at a given ratio, Peano-kernel constants). |
| [Gappa](https://gappa.gitlabpages.inria.fr/) | tool 1.6.0 in opam (upstream tags reach 1.8.3, 2026-09-09); `coq-gappa` (gappalib-coq) 1.11.0 | opam package 2025-11; 2026-06-02 | any | Automatic rounding-error bounds with Rocq certificates. **T1f** (coefficient rounding). Needs gmp, mpfr and boost on the system. `coq-gappa` declares no opam dependency on the tool, and it was not checked that 1.11.0 accepts tool 1.6.0 (its NEWS lists support for tool versions up to 1.8). |
| [VCFloat](https://verinum.org/vcfloat/) | 2.4.2 | 2026-09-22 | 9.0–9.2 | Automated round-off analysis of float expressions, as used by VST proofs ([VCFloat2, CPP 2024](https://dl.acm.org/doi/10.1145/3636501.3636953)). Pulls in CompCert 3.18. Optional for T7 and T6b. |
| [LAProof](https://github.com/VeriNum/LAProof) | 2.0.1 | 2026-09-25 | 9.1–9.2 | Proven accuracy of dot products (plain and FMA), gemv, gemm, sums, Cholesky and triangular solves, on MathComp matrices ([ARITH 2023](https://doi.org/10.1109/ARITH58626.2023.00021)). No LU/GEPP. Building blocks for T6b. Depends on VCFloat, libValidSDP and VST. |
| [ValidSDP / libValidSDP](https://github.com/validsdp/validsdp) | 1.1.2 | 2026-09-23 | stdlib < 9.3 | Floating-point Cholesky error bounds (Roux). `coq-libvalidsdp` 1.1.2 accepts MathComp < 2.7 and is a dependency of LAProof 2.0.1; only the full `coq-validsdp` (needs `osdp` and `coqeal`) is capped at MathComp < 2.6. Not a direct dependency of the plan. |
| CoqEAL | 2.1.2 | — | 9.0–9.3 | Refinements from proof-friendly to efficient data. Not needed for the pilot. |
| [rocq-num-analysis](https://lipn.univ-paris13.fr/rocq-num-analysis/) | 2.2.0 | 2026-03-11 | — | Lebesgue, Lax–Milgram and FEM (Boldo, Clément, Mayero et al.). Background, not a dependency. |
| [taylor_rocqs](https://github.com/holgerthies/taylor_rocqs) | 1.1 | 2026-09-23 | 9.2 only | Verified Taylor-series ODE solver, exact reals or CoqInterval back ends ([CPP 2026](https://doi.org/10.1145/3779031.3779097)). Comparison point, not a dependency. |
| [rocq-verified-extraction](https://github.com/MetaRocq/rocq-verified-extraction) | 1.0.0+9.1 | — | ≥ 9.1 | MetaRocq verified extraction to OCaml through Malfunction ([PLDI 2024](https://doi.org/10.1145/3656379)). Supports primitive floats and arrays. Research-grade install. |

### 1.5 Prior work

**In Rocq/Coq:**

- **Wave equation, PDE to C program.** Boldo, Clément, Filliâtre, Mayero, Melquiond and Weis. They proved convergence of the finite-difference scheme in Coq, then the round-off of the actual C code with Frama-C, Why3, Flocq and Gappa ([method error, arXiv 1005.0824](https://arxiv.org/abs/1005.0824); JAR 50 (2013); ["Trusting computations", CAMWA 2014, arXiv 1212.6641](https://arxiv.org/abs/1212.6641)). This is the reference "full stack" example. It took a team of experts years for one explicit scheme on one PDE.
- **Round-off of explicit Runge–Kutta.** Boldo, Faissole and Chapoutot, IEEE TC 69(12), 2020 (online 2019), including underflow and overflow. Faissole formalized it in Coq for Euler and RK2 on linear systems ([JAR 2023](https://link.springer.com/article/10.1007/s10817-023-09686-y)).
- **Verified leapfrog ODE solver.** Kellison and Appel, NSV 2022 ([paper](https://link.springer.com/chapter/10.1007/978-3-031-21222-2_9)). A C program verified with VST, a float functional model with VCFloat, and a real functional model with a global error bound. Their two-model method (float model ↔ real model) is the template recommended here.
- **Jacobi iteration.** Tekriwal, Appel, Kellison, Bindel and Jeannin, CICM 2023 ([paper](https://link.springer.com/chapter/10.1007/978-3-031-42753-4_14)): correctness, accuracy and convergence.
- **Lax equivalence theorem for finite differences.** Tekriwal, Duraisamy and Jeannin, NFM 2021 ([arXiv 2103.13534](https://arxiv.org/abs/2103.13534)).
- **Newton's method.** Paşca formalized Kantorovich's theorem and Newton with rounding at each step in Coq/SSReflect, around 2010 ([code](http://www-sop.inria.fr/marelle/Ioana.Pasca/code/)). It predates MathComp 2 and is not maintained.
- **Picard iteration.** Makarov and Spitters, ITP 2013 ([paper](https://doi.org/10.1007/978-3-642-39634-2_34)), on constructive reals.
- **Taylor models and power series for ODEs.** Park and Thies, ITP 2024 ([paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ITP.2024.30)), and Thies, CPP 2026 (taylor_rocqs).
- **Approximation models.** Bréhard, Mahboubi and Pous (`coq-approx-models`).

**In other provers:**

- **Isabelle/HOL.** Immler's verified ODE solver certified Tucker's Lorenz-attractor computations ([JAR 61 (2018)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6044317/)). It uses Runge–Kutta with affine arithmetic and generates SML code from the proofs. "Verifying Numerical Methods with Isabelle/HOL" (Bryant, Huerta y Munive, Foster; [arXiv 2511.20550](https://arxiv.org/abs/2511.20550)) covers bisection, fixed-point iteration and gradient descent.
- **Lean/Mathlib.** ODE existence and uniqueness (Picard–Lindelöf) and Gronwall-type estimates. No Flocq-class float library.

**No formal verification of a BDF method or of any stiff solver** was found in any prover (a second search found none either; this is a negative search result, not a proof of absence). Rigorous ODE work targets enclosures (interval, affine or Taylor models) with explicit methods. vstiff would be breaking new, if modest, ground.

### 1.6 Variable-step BDF2 stability results that matter here

- **Zero-stability.** The sufficient condition r = h_n/h_{n-1} < 1 + √2 ≈ 2.414 is sharp: for any constant ratio of 1 + √2 or more, solutions of the homogeneous recurrence do not stay bounded. A single larger ratio, or ratios that alternate (3, 1/3, 3, …), do not by themselves destroy stability, so for an individual grid the condition is sufficient, and sharp only as a uniform bound. See Grigorieff, Numer. Math. 42 (1983); Crouzeix and Lisbona, SINUM 21 (1984); Hairer–Nørsett–Wanner I, p. 405. All three are cited as such in [Akrivis et al., BIT 2024](https://www.cs.uoi.gr/~akrivis/ACHYZ.pdf) (checked; the paper is "The variable two-step BDF method for parabolic equations" and its own bound is r ≤ 1.9398).
- **Classical energy bounds for stiff dissipative (self-adjoint parabolic) problems:**
  - Grigorieff: (1+√3)/2 ≈ 1.366;
  - Becker: (2+√13)/3 ≈ 1.8685;
  - Emmrich: 1.9104;
  - Akrivis, Chen, Han, Yu and Zhang (BIT 2024): ≈ 1.94, using two multipliers.

  **None covers vstiff's ratio 2.**
- **Estimates "of a different kind".** These put (u₁ − u₀)/k₁ on the right-hand side:
  - Liao and Zhang ([arXiv 1912.11182](https://arxiv.org/abs/1912.11182)): unconditional stability for r ≤ (3+√17)/2 ≈ 3.561;
  - Li and Liao ([SINUM 2022](https://doi.org/10.1137/21M1462398)): stable "on arbitrary time grids", via discrete orthogonal convolution kernels.

  These do cover ratio 2, but are much harder to formalize. Their bounds carry mesh-dependent factors, so they do not contradict the mesh-independent 1 + √2 limit; the Li–Liao preprint ([arXiv 2201.00527](https://arxiv.org/abs/2201.00527)) itself quotes 1 + √2 as the zero-stability limit of variable-step BDF2.

---

## 2. The bridge to OCaml

### 2.0 What the compiled code computes: FMA contraction

OCaml's manual says so explicitly (native-code chapter, "Compatibility with the bytecode compiler"):

> On ARM and PowerPC processors … fused multiply-add (FMA) instructions can be generated for a floating-point multiplication followed by a floating-point addition or subtraction, as in `x *. y +. z`.

The OCaml 5.5.0 backends `asmcomp/{arm64,power,riscv,s390x}/selection.ml` all contain the pattern. `amd64` does not.

What `ocamlopt` 5.5.0 produced for vstiff on arm64 (checked with `ocamlopt -S` and `otool -tv` on the dune-built objects):

| Module | Expression | arm64 | Numerical effect |
|---|---|---|---|
| [`Vec.axpy`](../../src/vec.ml) | `(a *. x.(i)) +. y.(i)` | `fmadd` | Used in three places: (1) the Newton line search with λ = 2⁻ᵏ, where the product is exact; (2) the BDF2 predictor: ω is a power of two (product exact) except on the step after the `dt_max` cap engages, on a cut final step and, since `dc3bb67`, wherever `t + dt` rounds (product inexact; it only moves Newton's starting guess); (3) the first-step error predictor `y + h f`, where the product is inexact and feeds the accept/reject decision. |
| [`Jac.forward`](../../src/jac.ml) | `yk +. step yk`, with `step` inlined to `1e-8 *. (1. +. abs yk)` | `fmadd` | The perturbed coordinate. The code divides by the representable difference anyway. |
| [`Stage.solve`](../../src/stage.ml) | `(if i = j then 1. else 0.) -. (gamma *. v)` | `fmsub` | Entries of I − γJ. Off-diagonal entries are unaffected (0 − p is exact). |
| [`Newton.solve`](../../src/newton.ml) | `(1. -. (armijo *. lambda)) *. r` | `fmsub` | Armijo threshold. |
| [`Bdf2.coeffs`](../../src/bdf2.ml) | the denominator shared by the coefficients (a product by 2 added to 1) | `fmadd` | None: the product by 2 is exact. |
| [`Linalg`](../../src/linalg.ml) | `row.(j+1) -. (m *. top.(j+1))` | `fmul` + `fsub` | Not fused. Likely cause: both operands are bounds-checked loads, so instruction selection evaluates the right operand (the product) first into a register and the multiply is no longer visible to the pattern. A small example reproduces it: `a.(i) -. 2. *. b.(i)` and `a.(i) +. 2. *. b.(i)` do not fuse, `2. *. b.(i) +. a.(i)` does, and all of them fuse with `-unsafe`; `Linalg` fuses with `-unsafe` too (six sites in vstiff). `-inline 200` adds fused copies of `Vec.axpy` inside `Bdf2`. So which expressions fuse depends on flags, inlining and operand order as well as on the architecture. |

Two experiments:

- **Removing all fusions changes nothing observable.** vstiff was rebuilt with all five products wrapped in `Sys.opaque_identity`, which removes every fused instruction. `dune test` still passed. The Robertson (tol 1e-6), logistic (dt0 0.5) and van der Pol runs gave **bit-identical final states and identical step counts** (2531/831, 1109/374, 4145/1349). So FMA is semantically present but has been numerically invisible on the corpus. A second run agrees, also bit-for-bit for Robertson at tol 1e-9 (79869 accepted steps). There is a structural reason: the fused sites are the Jacobian entries and perturbation, Newton's starting point, and two discrete tests (the Armijo threshold and the first-step error estimate). None touches the residual or ψ, which fix Newton's root, so output changes only if a test is borderline.
- **How the operations are written decides fusion.** `a +. b *. c` written inline fuses. A let-bound product does not. The same computation through `[@@ocaml.inline always]` helper functions, the shape of Rocq's extracted `Float64.add`/`Float64.mul`, does not fuse either. Neither does `Sys.opaque_identity`. This is compiler behaviour, not a language guarantee.

Consequences:

- Any binary64 theorem about vstiff is **per architecture and per compiler flags** (`-unsafe` and inlining change which expressions fuse) unless the code pins the operation shape.
- Flocq can model both semantics: `Bfma` versus `Bplus ∘ Bmult`.
- Real-arithmetic theorems are immune to this. That is one reason to start there.

### 2.1 Option A, recommended: model in Rocq, test the OCaml against it

Write a small Rocq **mirror** of the OCaml functions under study: [`Bdf2.coeffs`](../../src/bdf2.ml), [`Halving.accepted`](../../src/halving.ml) and `Halving.rejected` (with its floor, [`Clock.resolution`](../../src/clock.ml)), and the attempt logic of [`Adaptive.integrate`](../../src/adaptive.ml). Keep it syntactically close to the OCaml and generic over a numeric structure, so one definition can be read three ways:

| Instance | Purpose | Example |
|---|---|---|
| `R : realFieldType` (MathComp) | Mathematical theorems | T1–T3, T9 |
| Flocq `binary_float 53 1024`, with `Bfma` where arm64 fuses | IEEE theorems | T1f, T4, T10 |
| `PrimFloat.float` | Execute in Rocq (`vm_compute`) and pin outputs as Rocq expected-output tests (dune Rocq lang ≥ 0.12) | Differential checks against the OCaml corpus |

This is the two-model method of Kellison and Appel, minus the C verification layer.

**Keeping mirror and code in sync.** A comment-free canonical source of each mirrored OCaml module acts as a tripwire. Prototype, tested on the trees before and after the comment-only rewrite:

```ocaml
(* sync/canon.ml: print an implementation without comments or doc
   comments, so that only code changes alter the output. Needs only
   compiler-libs from an OCaml >= 5.5 switch (modular explicits). *)
let is_doc (a : Parsetree.attribute) =
  match a.attr_name.txt with "ocaml.doc" | "ocaml.text" -> true | _ -> false

let strip =
  let open Ast_mapper in
  { default_mapper with
    attributes = (fun m l -> default_mapper.attributes m (List.filter (fun a -> not (is_doc a)) l));
    structure = (fun m items ->
      default_mapper.structure m
        (List.filter (fun (i : Parsetree.structure_item) ->
             match i.pstr_desc with Pstr_attribute a -> not (is_doc a) | _ -> true) items)) }

let () =
  let ast = Pparse.parse_implementation ~tool_name:"canon" Sys.argv.(1) in
  Format.printf "%a@." Pprintast.structure (strip.structure strip ast)
```

A dune rule in the **main workspace** (a separate `--root` project cannot read `../../src`, see 4.4) diffs `canon %{dep:../src/halving.ml}` against a promoted `halving.ml.canon`:

- A **code** change shows up as a diff, which means re-checking the Rocq mirror and then promoting.
- A comment change shows up as nothing: the comment-only rewrite of commit `10ff44d` changed the comments across `src/` and `test/`, and the canonical digests of all 19 `.ml` files are identical before and after.
- Caveats: `Pprintast`'s output format is not guaranteed stable across compiler releases, so a compiler upgrade can force re-promoting every snapshot (pin the compiler in CI). The tool reads `.ml` files only, so an `.mli` contract change is seen only if the `.ml` changes too.

**Strengths:**

- No change to the OCaml code, the build or the dependencies.
- Proofs survive bit-level refactors.
- It works today, with every OCaml 5.5 feature vstiff uses.

**Weakness:** the theorem is about the mirror. The OCaml can differ in ways the tripwire does not catch, such as a slip made while editing both. Differential tests (PrimFloat instance against OCaml) and theorem-derived property tests (section 4.2) shrink that gap. Keeping one definition generic over the real, Flocq and PrimFloat instances is the main design risk of this option; the pilot needs only the real instance.

### 2.2 Option B: write kernels in Rocq and extract them to OCaml

- **Floats.** `PrimFloat.add`, `mul` and the rest extract to `Float64.add` and friends. These are `[@@ocaml.inline always]` wrappers over `+.`. They live in Rocq's kernel and are packaged for extracted code as `rocq-primitive` 9.0.0 (LGPL-2.1-only, released from the peregrine-project GitHub organization). That would be a new runtime dependency. The alternative is to remap them with `Extract Inlined Constant PrimFloat.add => "( +. )"`, which is a small unverified mapping. Rocq 9.3's `PrimFloat` has **no `fma`**, and no `min` or `max` (define them from comparisons).
- **Arrays and integers.** `PArray` extracts to Rocq's persistent arrays; lists would make the O(n³) solve O(n⁴). `nat` must not be extracted as unary numbers: use `PrimInt63`, or the unverified `nat → int` mapping.
- **Placement in the module structure.** Extraction would produce something like `halving_core.ml`/`.mli`. The hand-written [`halving.ml`](../../src/halving.ml) would keep its labels and the [`Ode.Controller`](../../src/ode.mli) signature and delegate:

  ```ocaml
  let init ~tol ~dt0 ~dt_max ~max_rejects = Halving_core.init tol dt0 dt_max max_rejects
  ```

  Modular explicits, labels and module types stay in hand-written OCaml. Building through dune's `rocq.extraction` needs Rocq in the main build switch, which rules out MathComp there (1.2). So either commit the generated files and regenerate them in CI, or limit extracted theories to Rocq + Stdlib + Flocq. Rocq 9.3 + Flocq does solve with dune 3.24.2.
- **Bit-identity.** It would hold for [`Halving`](../../src/halving.ml) and for [`Bdf2.coeffs`](../../src/bdf2.ml), where the one fused product (by 2) is exact. It would **not** be guaranteed for [`Vec.axpy`](../../src/vec.ml), [`Stage`](../../src/stage.ml), [`Newton`](../../src/newton.ml) or [`Jac`](../../src/jac.ml). PrimFloat computes two roundings where the current arm64 binary fuses, which conflicts with the "refactors bit-identical" rule. In the experiment nothing visible changed, and there is a structural reason (section 2.0: the fused sites do not touch the residual), but not a guarantee.
- **Trust.** Classic extraction is unverified. MetaRocq verified extraction (PLDI 2024) closes the extraction gap but goes through the same OCaml backend, so FMA selection still applies. Its documented install path is an OCaml 4.14 flambda switch with pinned MetaRocq and Malfunction.
- **Verdict.** This is reasonable only for discrete control logic such as `Halving`, and even there it adds generated code and a dependency for little gain over Option A. Not recommended now.

### 2.3 Option C: translate the OCaml into Rocq

| Tool | Status (2026-10) | Fit for vstiff |
|---|---|---|
| [rocq-of-ocaml](https://github.com/formal-land/rocq-of-ocaml) (formerly coq-of-ocaml) | Master May 2026 requires OCaml 5.4–5.6 **and** Rocq 9.0.x, but Rocq 9.0 excludes OCaml 5.5, so it only works on 5.4. opam has `coq-of-ocaml` 2.5.3 only as +4.12, +4.13 and +4.14 builds (the last needs `ocaml >= 4.14 & < 4.15`). | **Unusable for numerics.** `proofs/Basics.v` defines `float := Z` and `array A := list A`, and `src/constant.ml` "approximates" float literals by integers (`1e-8` becomes 0). It is limited to OCaml 5.4 by its Rocq 9.0 dependency, so it cannot read OCaml 5.5 source such as modular explicits. |
| [CFML](https://github.com/charguer/cfml) | Needs OCaml < 5 and Coq 8.20. Last push March 2025. | Separation logic for imperative code. No real float support. Wrong tool for pure numeric code. |
| Osiris/Horus ([ICFP 2025](https://iris-project.org/pdfs/2025-icfp-osiris.pdf)) | Research prototype: a deep embedding of an OCaml fragment in Rocq/Iris, with Horus, a Hoare logic for pure fragments | Worth watching. Floats and arrays are unclear. Not usable now. |
| Hand-written mirror (Option A) | — | **What remains.** |

The subset vstiff uses:

- **Fits any functional translator:** records, variants, `result` with binding operators, local recursion, `Array.init`/`map`/`fold_left`.
- **Fits none of the current ones:** `float` semantics, the `Float.*` functions (`max`/`min` NaN propagation, `round`, `to_int`), modular explicits, first-class module packing, `Printf`, and `invalid_arg` in [`Check`](../../src/check.ml).

### 2.4 Option D (non-Rocq): Gospel, Ortac, Cameleer and Why3

- **[Gospel](https://github.com/ocaml-gospel/gospel).** Specifications inside `.mli` comments. Release 0.3.1 (2026-03-18) is bounded to OCaml ≤ 5.3.0 in opam; the main branch relaxes this to ≥ 4.14. Whether it parses [`stepper.mli`](../../src/stepper.mli) and [`adaptive.mli`](../../src/adaptive.mli) (modular explicits) was not verified. Gospel has a built-in `float` type but no real-number theory in its standard library, so accuracy specifications would need user-declared logic.
- **[Ortac](https://github.com/ocaml-gospel/ortac).** 0.8.0, 2026-03-23, pinned to Gospel 0.3.1, so OCaml ≤ 5.3. Turns Gospel specifications into runtime checks and QCheck-STM tests. It would be a cheap middle level between expect tests and proofs once it supports 5.5.
- **[Cameleer](https://github.com/ocaml-gospel/cameleer).** Install by pinning from GitHub. Generates Why3 verification conditions for OCaml with Gospel specifications and discharges them with SMT solvers. Good for integer and structural contracts: [`Halving`](../../src/halving.ml) counters, the termination measure, `streak mod 3`. Weak for real and float numerics. Why3 1.8.2's Rocq realizations stop at Rocq 9.0 (`why3-coq`: `coq-stdlib < 9.1~`).
- **Verdict.** Interesting for contracts on the `.mli` files, but blocked on OCaml 5.5 support today, and it does not reach the numerical theorems that matter.

---

## 3. Candidate theorems, ranked by value for effort

Effort is in focused person-days for someone already comfortable with MathComp; add a one-time ramp-up of 1 to 2 weeks for a newcomer who already writes Rocq proofs (a guess that could not be checked; several weeks if Rocq itself is new). "Real" means a theorem over an ordered or real field; "IEEE" means over binary64.

| Rank | Id | Theorem (one line) | Arithmetic | Libraries | Effort | Value |
|---|---|---|---|---|---|---|
| 1 | T1 | Variable-step BDF2 coefficients are exact on quadratics; residual on cubics is −β h²(h + h_prev) | Real | MathComp algebra | 0.5–1 d | Medium: pins the formula; template for variable order |
| 2 | T2 | Variable-step BDF2 is zero-stable on every grid whose ratios are all ≤ ω* < 1 + √2 (sharp for constant ratios): difference contraction q(ω) = −a₀(ω); explicit perturbation bound | Real | MathComp algebra | 2–4 d | High, with T3 |
| 3 | T3 | [`Halving`](../../src/halving.ml) + [`Adaptive`](../../src/adaptive.ml) + [`Bdf2`](../../src/bdf2.ml) never attempt ω > 2 (cross-module; real arithmetic, and binary64 up to the half-ulp snapping of steps) | Real, then IEEE | MathComp (+ Flocq) | 2–4 d real, +2–4 d binary64 | High |
| 4 | T4 | Termination of `Adaptive.integrate`. **False at `f6d9b4e`** (finding 1); true since `dc3bb67`, with at most (max_rejects + 1)(N + 1) attempts, N ≈ 2⁶² on (0, 1] | IEEE | Flocq | 1–3 wk | High |
| 5 | T6a | [`Linalg.solve`](../../src/linalg.ml) correct over an ordered field: `Some x ⇒ A x = b`; `None ⇒ det A = 0`; `det A ≠ 0 ⇒ Some` | Real | MathComp (`matrix.v`, `cormen_lup` as guide) | 1–3 wk | Medium-high |
| 6 | T9 | Linear invariants: if c·f ≡ 0 then c·yₙ = c·y₀ exactly for BDF1/BDF2 with *inexact* damped Newton and forward-difference Jacobians | Real | MathComp; needs only the soundness half of T6a | 1–2 wk | Medium-high |
| 7 | T5 | Linear stability: backward Euler contractive for Re z ≤ 0; constant-step BDF2 A-stable; for linear f the code path is exact in one Newton step | Real/complex | MathComp (`algC`, or `real_closed` complex) | 1–2 wk | Medium |
| 8 | T10 | IEEE safety: every accepted state is finite. Counterexamples on record: [`Newton.solve`](../../src/newton.ml) could return `Ok` with an infinite entry (fixed in `052fdc7`); [`Stepper.fixed`](../../src/stepper.ml) "lands exactly" is false in binary64 (its `.mli` no longer claims it) | IEEE | Flocq | 3–7 d | Medium |
| 9 | T11 | Error-estimator asymptotics: the estimate is first order, −(h²/2)y'' + O(h³) (BDF2 minus backward Euler), while the BDF2 local error is −(β/6)h²(h + h_prev)y''' + O(h⁴) | Real | Coquelicot or MathComp-Analysis | 1–2 wk | Medium: quantifies the "conservative estimate" gap |
| 10 | T1f | Binary64 coefficient defect: computed a1 + a0 is within a few ulps of 1 for ω ∈ (0, 2]; with T2, a per-step drift bound for y' = 0 | IEEE | Flocq + Gappa | 1–3 d after float setup | Low-medium |
| 11 | T7 | Forward-difference Jacobian error in binary64, FMA-aware: truncation (M/2)h̃ⱼ + amplified evaluation error (ε+ε′)/h̃ⱼ + O(u)·\|J̃ᵢⱼ\|; exact step when \|yⱼ\| ≳ 2e-8 | IEEE + real | Flocq + Coquelicot or MathComp-Analysis (+ VCFloat) | 2–4 wk | Medium |
| 12 | T8 | Damped Newton local convergence with an inexact Jacobian: full steps pass the Armijo test near a regular root; linear rate ≈ βη | Real | MathComp-Analysis | 3–8 wk | Medium |
| 13 | T6b | GEPP backward error in binary64: (A + ΔA)x̂ = b, \|ΔA\| ≤ γ₃ₙ\|L̂\|\|Û\|, multipliers \|m\| ≤ 1 exactly | IEEE | Flocq + LAProof + VCFloat | 2–4 mo | Medium (textbook) |
| 14 | T12 | Stability of variable-step BDF2 at ω ≤ 2 for dissipative linear problems (λ ≤ 0), via discrete orthogonal convolution kernels | Real | MathComp-Analysis | Months, research-level | Low-medium |
| 15 | — | End-to-end binary64 correctness of the adaptive integrator (global error ≤ C·tol) | IEEE | — | Impractical | — |

All sketches below are **not machine-checked**. They use MathComp 2.6 conventions. Tactic and module names follow MathComp 2.6's source (`field_tactic`, `arithmetic_tactic`) but have not been run.

### T1. BDF2 coefficients are exact on quadratics (consistency of order 2)

**Math.** Let ω > 0, h_prev arbitrary, h = ω h_prev, and let a₁, a₀, β be the three coefficients that [`Bdf2.coeffs`](../../src/bdf2.ml) returns for ω. No closed form is repeated here: `Bdf2.coeffs` is their single source, and [numerics/04-bdf.md](../numerics/04-bdf.md) sets up the equations that fix them. Then for every polynomial p with deg p ≤ 2:

p(tₙ + h) = a₁ p(tₙ) + a₀ p(tₙ − h_prev) + β h p′(tₙ + h).

The three exactness conditions (for p = 1, t − tₙ₊₁ and (t − tₙ₊₁)²) fix the coefficients; the first is a₁ + a₀ = 1.

On the cubic (t − tₙ₊₁)³ the residual is exactly a₁h³ + a₀(h + h_prev)³ = −β h² (h + h_prev). By Taylor's theorem this gives a local truncation error of order −(β/6) h²(h + h_prev) y‴. No condition on h_prev is needed; ω > 0 is enough for the coefficients to be defined.

The identities were checked with exact rationals at 2000 random points, and the cubic residual at ω ∈ {1/4, 1/2, 1, 2}.

```coq
(* Sketch, not machine-checked. *)
From mathcomp Require Import all_ssreflect all_algebra field_tactic arithmetic_tactic.
Import GRing.Theory Num.Theory.
Local Open Scope ring_scope.

Section Bdf2Coeffs.
Variable R : realFieldType.
Implicit Types w hp tn p0 p1 p2 : R.

(* Mirror of Bdf2.coeffs (src/bdf2.ml): den, a1, a0 and beta : R -> R are transcribed from it, with its
   operation order. Their expressions are deliberately not repeated in this document. *)

Lemma den_neq0 w : 0 < w -> den w != 0.
Proof. by move=> w0; apply/lt0r_neq0; rewrite /den; lra. Qed.

Theorem bdf2_exact_on_quadratics w hp tn p0 p1 p2 : 0 < w ->
  let h := w * hp in
  let p t := p0 + p1 * t + p2 * t ^+ 2 in
  let dp t := p1 + 2 * p2 * t in
  p (tn + h) = a1 w * p tn + a0 w * p (tn - hp) + beta w * h * dp (tn + h).
Proof.
move=> w0 /=; have := den_neq0 w0.   (* `field` closes denominators with `done`: keep the fact in context *)
by rewrite /a1 /a0 /beta /den; field.
Qed.

Theorem bdf2_cubic_residual w hp : 0 < w ->
  let h := w * hp in
  a1 w * h ^+ 3 + a0 w * (h + hp) ^+ 3 = - (beta w * h ^+ 2 * (h + hp)).
Proof.
move=> w0 /=; have := den_neq0 w0.
by rewrite /a1 /a0 /beta /den; field.
Qed.
End Bdf2Coeffs.
```

Note on the sketches here and in T2: in MathComp 2.6.0 plain `field` fails with "There are remaining goals" unless `done` closes the denominator side conditions, so `field; exact: den_neq0` would not run. The non-zero fact has to be in the context (or use `field by ...`), and `den` has to be unfolded so that `field` sees the denominator. Still unchecked.

A stronger form states it for `p : {poly R}` with `size p <= 3` and `deriv p`. The forward-looking version: when vstiff gains variable order, a Rocq function computing variable-step BDFk coefficients, proven exact on degree-k polynomials, is the reusable artifact.

### T2. Zero-stability of variable-step BDF2

**Math.** Since a₁ = 1 − a₀:

yₙ₊₂ − yₙ₊₁ = q(ωₙ)(yₙ₊₁ − yₙ) + δₙ,  where q(ω) = −a₀(ω), the parasitic root of the frozen-coefficient recurrence ([numerics/04-bdf.md](../numerics/04-bdf.md)).

Here δₙ is any perturbation: rounding, the stage residual, or the f-term for convergence proofs (for stiff f the f-term is not small, so this is the h → 0 notion of stability). For ω > 0, q(ω) < 1 ⟺ ω < 1 + √2 (exercise 4 of numerics/04-bdf.md finds this crossing).

If q(ωₙ) ≤ q* < 1 for all n, then:

|yₙ| ≤ |y₀| + (|y₁ − y₀| + Σₖ₍ₙ |δₖ|)/(1 − q*).

The bound applies componentwise, hence in the max norm. For vstiff (ω ≤ 2, by T3), q* = q(2), the value of −a₀ at ω = 2 (q increases with ω, so this is the worst case), and the factor is 1/(1 − q(2)); evaluate `Bdf2.coeffs 2.` to see it. The factor 1/(1 − q*) is attained by the homogeneous recurrence with a constant ratio. The threshold is sharp as a uniform condition: for a constant ratio ω ≥ 1 + √2, q ≥ 1 and the differences stop decaying (they grow geometrically above it). A grid with some ratios above it can still be stable, for example one alternating between 3 and 1/3, so the hypothesis is sufficient, not necessary, for an individual grid. At the threshold, q(1 + √2) = 1.

```coq
(* Sketch, not machine-checked; continues Section Bdf2Coeffs. *)
Definition q w := - a0 w.

Lemma bdf2_increment w yn ym : 0 < w ->
  a1 w * yn + a0 w * ym - yn = q w * (yn - ym).
Proof. by move=> w0; have := den_neq0 w0; rewrite /a1 /a0 /q /den; field. Qed.

Lemma q_lt1 w : 0 < w -> (q w < 1) = (w < 1 + Num.sqrt 2).
(* `Num.sqrt` exists only over `rcfType`, so state this one over `R : rcfType`, or keep a
   polynomial form of the threshold. *)

Theorem bdf2_zero_stable (w y d : nat -> R) (qs : R) :
  (forall n, 0 < w n) -> (forall n, q (w n) <= qs) -> 0 <= qs < 1 ->
  (forall n, y n.+2 = a1 (w n) * y n.+1 + a0 (w n) * y n + d n) ->
  forall n, `|y n| <= `|y 0| + (1 - qs)^-1 * (`|y 1 - y 0| + \sum_(k < n) `|d k|).
```

### T3. vstiff never attempts a BDF2 step with ω > 2 (cross-module)

**Claim.** Consider every call `Bdf2.step_with_error rhs h (After {h_prev; _}) at` made by `Adaptive.integrate (module Bdf2) (module Halving)`. In real arithmetic each satisfies h ≤ 2·h_prev, so ω ≤ 2, with one exception: a final step that the driver takes because the remainder is at most [`Clock.resolution t`](../../src/clock.ml). That step ends the run, and it exceeds 2·h_prev only if h_prev is below half of that resolution.

**Why.** The history's `h_prev` is the length of the last *accepted* step. A final step ends the integration, so `h_prev` equals the controller's proposal at that time. After an accept, the proposal is `dt` or `min(2 dt, dt_max)`, so at most 2·h_prev. Each rejection that does not end the run sets `dt := h/2` with h ≤ `dt`, so it only shrinks (a rejected remainder step may exceed `dt`, but half of it is below `Clock.resolution t`, so it ends the run). The attempted h is the proposal or a final step: the cut `t_end − t` (at most the proposal) or, since `dc3bb67`, a remainder at or below `Clock.resolution t` (16 to 32 ulp of `t`, since an ulp lies between eps|t|/2 and eps|t|) even when it exceeds the proposal. Snapping a non-final step to the floats is the identity in real arithmetic. The invariant spans [`Halving.accepted`](../../src/halving.ml)/`rejected`, [`Adaptive.go`](../../src/adaptive.ml) (cut and history threading) and [`Bdf2.step_with_error`](../../src/bdf2.ml) (history update).

**Binary64.** Until `dc3bb67` it also held in binary64:

- `2.*.dt` is exact or overflows, in which case `min` returns `dt_max ≤ max_float < 2·h_prev`;
- RN(h/2) ≤ h;
- the cut h = RN(t_end − t) ≤ `dt` by the code's own comparison;
- RN(h/h_prev) ≤ RN(2) = 2.

Since `dc3bb67` a non-final step is snapped to the floats, h = (t + dt) − t, which can differ from `dt` by up to half an ulp u of the new time `t + dt`. The first three points still hold for the proposal, but the fourth now gives only ω ≤ 2 + 3u/(2·h_prev): the snapped step is at most `dt` + u/2, and the previous proposal at most h_prev + u/2 (the ulp of `t` only grows along a run, so the same u covers the previous step). The excess is zero when `t + dt` is representable (every step of a dyadic grid, for example) and negligible for steps of many ulps. For h_prev ≥ 4u it keeps ω ≤ 2.375, still below 1 + √2. A step of a few ulps can have a larger ratio: a proposal of 1.4u snaps to u, and the doubled proposal 2.8u to 3u.

**Hypothesis the proof surfaces.** h_prev > 0, i.e. no accepted step of length 0. Since `dc3bb67` it holds: a step with h ≤ 0 is rejected as `Too_small` without calling the method, so only steps of positive length are accepted. At `f6d9b4e` a zero-length step could be accepted after `dt` underflowed to 0. From a normal `dt0` that takes about a thousand consecutive rejections at t = 0, where the 16 eps |t| floor vanishes (1055 for the default on a span of 1, 1072 for `dt0 = 1`, and the default `max_rejects` of 50 would stop the run first); a user-supplied `dt0 = 5e-324` gets there in one. With an `rhs` that is NaN for t > 0, the first attempt was rejected, the retry had h = 0 and was accepted, and the next BDF2 step had ω = 0/0 = NaN; the run ended `Error StepRejected 51` after 106 `rhs` calls. Since the fix the same run still ends `Error StepRejected 51`, through `Too_small` rejections, and the method is never called with h = 0.

**Empirical check.** `Bdf2` was wrapped to record ω. At `f6d9b4e` the maximum ω was exactly 2 on logistic, van der Pol and Robertson (358, 1362 and 842 attempts at ω = 2). Ratios were powers of two except when the `dt_max` cap engaged and on a cut final step (1 to 3 such attempts in each of the three runs checked). Since `dc3bb67` the ratios are no longer exact powers of two wherever `t + dt` rounds: ω can then exceed 2 by the rounding above.

```coq
(* Sketch, not machine-checked. Mirror of Halving and of the attempt loop of Adaptive.integrate
   over an ordered field (exact doubling/halving; snapping a step to the floats is the identity in R;
   dt_min mirrors Clock.resolution). *)
Record ctl := Ctl { dt : R; streak : nat; failures : nat }.

Definition accepted (dt_max : R) (c : ctl) : ctl :=
  let s := (streak c).+1 in
  Ctl (if s == 3 then Num.min (2 * dt c) dt_max else dt c) (s %% 3) 0.

Definition rejected (max_rejects : nat) (dt_min : R) (c : ctl) (h : R) : option ctl :=
  let f := (failures c).+1 in
  if (max_rejects < f)%N || (h / 2 < dt_min) then None else Some (Ctl (h / 2) 0 f).

Inductive hist := Start | After of R.            (* After h_prev *)

(* The last step takes the whole remainder, when it is within the proposal or within Clock.resolution. *)
Definition final (dt_min : R -> R) (t_end t : R) (c : ctl) :=
  (t_end - t <= dt c) || (t_end - t <= dt_min t).

Definition attempt dt_min t_end t c :=           (* proposal, or the last step *)
  if final dt_min t_end t c then t_end - t else dt c.

Inductive move (t_end dt_max : R) (dt_min : R -> R) (mr : nat) :
    R * ctl * hist -> R * ctl * hist -> Prop :=
| Accept t c hs : t < t_end -> 0 < attempt dt_min t_end t c ->
    move t_end dt_max dt_min mr (t, c, hs)
      (if final dt_min t_end t c then t_end else t + attempt dt_min t_end t c,
       accepted dt_max c, After (attempt dt_min t_end t c))
| Reject t c c' hs : t < t_end ->
    rejected mr (dt_min t) c (attempt dt_min t_end t c) = Some c' ->
    move t_end dt_max dt_min mr (t, c, hs) (t, c', hs).
(* A step with attempt <= 0 is never accepted: Adaptive rejects it as Too_small, which is the Reject case. *)

Definition ratio_ok dt_min t_end (s : R * ctl * hist) :=
  let: (t, c, hs) := s in
  0 <= dt c /\ forall hp, hs = After hp -> 0 < hp /\
    (t < t_end -> attempt dt_min t_end t c <= 2 * hp \/ t_end - t <= dt_min t).

Theorem ratio_le_2 t_end dt_max dt_min mr s s' :
  0 < dt_max -> ratio_ok dt_min t_end s -> move t_end dt_max dt_min mr s s' -> ratio_ok dt_min t_end s'.
(* Corollary with T2: apart from a last step within the resolution, every BDF2 step sequence vstiff
   produces has q <= q 2. *)
```

**Design value.** A future controller (PI, Gustafsson) gets a precise obligation: clamp the ratio below 1 + √2 for zero-stability, or below about 1.94 to stay inside the classical stiff energy estimates. This could become part of the [`Ode.Controller`](../../src/ode.mli) contract.

### T4. Termination of `Adaptive.integrate`: false at `f6d9b4e`, true since `dc3bb67`

**What failed at `f6d9b4e` (runs).** y′ = −y, y₀ = 1, `t0 = 1.`, `t_end = Float.succ 1.`, `~tol:1e-6`, all other options default:

- `dt_max = span/10 = 2⁻⁵²/10`, below half an ulp of 1, so `1. +. h = 1.`;
- every step was accepted, so `failures` kept resetting;
- `last` never fired, because `dt < t_end − t`.

After 3e6 `rhs` calls the largest `t` passed to `rhs` was exactly `0x1p+0`. With the defaults this happened whenever t_end − t₀ is below about five ulps of t₀. It also happened for any user `dt_max` below half an ulp of t (e.g. t₀ = 1e10, `dt_max = 5e-7`). With `dt_max = 2e-6` it progressed, and a span of 2⁻⁴⁰ finished in 58 steps. (A call budget of 3e6 is too small for the `dt_max = 2e-6` case, which needed about 5e5 steps; there `t` did advance.)

**A second failure: state and clock disagreed (`f6d9b4e`).** [`Adaptive`](../../src/adaptive.ml) advanced the state by `h` but the clock only to `RN(t + h)`. When `h` is a few ulps of `t` the two disagree by tens of percent per step and nothing failed, so a guard on `t + h = t` alone would have missed it. Take `y' = 1`, `y(t0) = 0`, `t0 = 1e15`, `t_end = t0 + 100` (ulp(t) = 0.125), `tol = 1e-6`. `dt_max = 0.19` returned `Ok` at `t = t_end` with `y = 76.36` (should be 100); `dt_max = 0.125` gave 100.24; `dt_max = 3`, below [`Halving`](../../src/halving.ml)'s own floor 16 eps |t| = 3.55, gave 100.08; the defaults (`dt_max = 10`) gave 100.32. With `t0 = 0` and `dt_max = 0.19` the same run returned 99.99999999999966.

**The fix (`dc3bb67`).** In `Adaptive.integrate`, every non-final step is snapped to the floats, h = (t + dt) − t, so the state advances by exactly what the clock does. A snapped step can be up to half an ulp of the new time `t + dt` longer or shorter than the controller's proposal, so it can exceed `dt_max` by that much. A step that cannot move `t` (h ≤ 0) is rejected with the new `Ode.rejection` constructor `Too_small`, without calling the method, and the controller decides what happens: `Halving` halves it, which is below its floor, so it gives up with `StepRejected`. A remainder t_end − t at or below [`Clock.resolution t`](../../src/clock.ml) is taken as the last step. The module `Clock` holds the floor, `Clock.resolution t = 16 eps |t|`, which `Halving`'s give-up test now calls (the rule is unchanged). Corpus lines 25 to 27 pin it: one ulp from `t = 1` ends `Ok` at `t_end`; `t0 = 1e10` with `dt_max = 5e-7`, below the resolution of `t`, ends `Error StepRejected 1`; `t0 = 1e15` over 100 with `dt0 = dt_max = 0.19` ends `Ok` at `t_end` with `y = 100`. Without the snapping and the `Too_small` rule, lines 25 and 26 hit the call budget (`no answer within 5e6 rhs calls`) and line 27 printed `y = 76`; with snapping but without the remainder rule, line 25 printed `Error StepRejected 1`. The `t0 = 1e15` runs of the previous paragraph that use the default `dt0`, whose first step 1e-4 is below half an ulp (0.0625) of `t`, now end in `StepRejected 1` instead of a wrong `Ok`. The rules are explained in [numerics/05-step-control.md](../numerics/05-step-control.md).

**Alternatives.** The smallest change, failing when `t + h = t` on a non-final step, catches exact stagnation only. The same check with the floor `Halving` applies to rejections (fail on a non-final step with h < 16 eps |t|) would also have removed the second failure, and so would a controller signature that sees the time on accepts. The fix took neither route: snapping makes the state follow the clock, no floor applies to accepted steps, and no signature changed apart from the new constructor.

**Theorem (current code).** For finite t₀ < t_end and a total, terminating `rhs`, the loop returns after at most (max_rejects + 1)·(N + 1) attempts, where N is the number of binary64 values in (t₀, t_end]. Two facts give it. Every accepted non-final step strictly increases `t` (h > 0 means t_next > t), so there are at most N accepted steps, the last one landing on `t_end`. And the controller ends any run of rejections (`Halving`: the (max_rejects + 1)-th rejection in a row gives up, or the floor does), so at most max_rejects rejections precede each accepted step, and a final run of at most max_rejects + 1 attempts ends in `Error`. Rejections include the `Too_small` ones, which pass through the controller without calling the method. N is about 2⁶² on (0, 1], so this is a termination argument, not a usable bound. The measure is lexicographic: floats left, then rejections left. Flocq's `Bsucc`/ulp lemmas or `PrimFloat.next_up_equiv` give the float order.

**Note.** Here floats make the proof *easier*. In the real-number model, termination needs a liveness hypothesis on the method: every attempt with h ≤ h* is accepted. Without it, the 16 eps |t| floor (`Clock.resolution`) allows Zeno sequences near t = 0.

### T5. Linear stability

- **Backward Euler.** For Re z ≤ 0, |1/(1 − z)| ≤ 1. Trivial.
- **Constant-step BDF2 is A-stable.** The roots ζ of its characteristic polynomial (the constant-step case, ω = 1, of the recurrence of [`Bdf2.coeffs`](../../src/bdf2.ml)) satisfy |ζ| ≤ 1 for Re z ≤ 0, with ζ = 1 a simple root only at z = 0.

  Proof idea, about a page: put w = 1/ζ and u = 1 − w = a + ib. The stability relation expresses z through u, so Re z is a quadratic expression in a and b. If |ζ| > 1 then (1 − a)² + b² < 1, so b² < 2a − a², and Re z is then positive. Needs complex numbers: MathComp's own `algC` (package `rocq-mathcomp-field`, with `'Re` and `'Im`), `R[i]` from `rocq-mathcomp-real-closed`, or pairs of reals.
- **Code path.** For linear f(t, y) = A y the forward-difference Jacobian is exactly A in real arithmetic. Damped Newton from any guess then takes one full step to the exact stage solution: the Armijo test passes with residual 0. So the real model of [`Bdf1.step`](../../src/bdf1.ml) returns exactly (I − hA)⁻¹yₙ, and `Bdf2` returns (I − βhA)⁻¹(a₁yₙ + a₀yₙ₋₁). This ties the textbook stability results to vstiff's actual algorithm. Needs T6a and a Newton mirror.
- **Variable steps.** Stability for λ ≤ 0 with ω ≤ 2 is T12; see 1.6. It is research-level to formalize.

### T6a. `Linalg.solve` is correct over an ordered field

Mirror `solve_augmented` (max-|·| pivot, row swap, Schur update on the trailing rows, back-substitution in the recursion) as a `Fixpoint` on `'M[R]_(n, 1 + n)`. MathComp's `cormen_lup` has the same recursion shape (`xrow 0 k`, `dlsubmx`, `ursubmx`, `drsubmx`) and already proves P·A = L·U.

```coq
(* Sketch, not machine-checked. *)
Theorem solve_sound n (A : 'M[R]_n) (b : 'cV[R]_n) x : solve A b = Some x -> A *m x = b.
Theorem solve_none n (A : 'M[R]_n) b : solve A b = None -> \det A = 0.
Theorem solve_some n (A : 'M[R]_n) b : \det A != 0 -> exists x, solve A b = Some x.
```

**In binary64:** `None` iff a computed pivot is exactly 0. That is neither necessary nor sufficient for singularity, which is a clean example of a real theorem that does not transfer. The pivot choice is irrelevant over a field. It matters only for T6b.

### T9. Linear invariants are conserved exactly, even with inexact Newton

**Math.** Suppose c·f(t, y) = 0 for all t and y; Robertson has c = (1, 1, 1). Then, in exact arithmetic:

- c·(f(y + hⱼeⱼ) − f(y))/hⱼ = 0, so c·J_fd = 0 and c·(I − γJ_fd) = c;
- each Newton correction satisfies c·dx = −c·G(x) = −(c·x − c·ψ).

If the guess satisfies c·x₀ = c·ψ, every damped iterate keeps c·x = c·ψ, whatever λ and whenever Newton stops. Both guesses do: yₙ for BDF1, and yₙ + ω(yₙ − yₙ₋₁) for BDF2. And c·ψ = (a₁ + a₀)c·yₙ = c·yₙ by T1.

Hence every accepted state satisfies c·yₙ = c·y₀. This explains why the corpus's `|y1 + y2 + y3 − 1| < 1e-8` holds: only rounding drifts the mass. Note that a1 + a0 ≠ 1 in binary64 for about half of all ω (see T1f).

```coq
(* Sketch, not machine-checked. *)
Variables (n : nat) (f : 'cV[R]_n -> 'cV[R]_n) (c : 'rV[R]_n).
Hypothesis cf0 : forall y, c *m f y = 0.
Definition jac_fd (hs : 'I_n -> R) y : 'M[R]_n :=
  \matrix_(i, j) ((f (y + hs j *: delta_mx j 0) - f y) i 0 / hs j).
Theorem newton_step_keeps_invariant gamma psi hs x dx lam :
  (forall j, hs j != 0) ->
  (1%:M - gamma *: jac_fd hs x) *m dx = - (x - psi - gamma *: f x) ->
  c *m x = c *m psi -> c *m (x + lam *: dx) = c *m psi.
```

**Value.** The theorem guards against regressions. Examples: a predictor that breaks c·x₀ = c·ψ (it then matters only when damping engages, since a full Newton step restores the invariant); Jacobian reuse across steps (still fine, c·J = 0); a quasi-Newton update that breaks c·J = 0.

### T10. IEEE safety and the two contract findings

- **Theorem.** If [`Halving.acceptable`](../../src/halving.ml) accepts (y, err) with err = y − be (After history) or err = ½(be − (yₙ + h f)) with y = be (Start), then every yᵢ is finite. Reason: an infinite yᵢ makes |errᵢ|/(1 + |yᵢ|) NaN (∞/∞), `Float.max` propagates NaN, and `NaN <= tol` is false. So `Adaptive.integrate … = Ok s ⇒ Vec.finite s.y`.
  - Proof sketch: case analysis over Flocq's special values, plus a model of `Float.max` and `Float.abs`.
- **Findings, both run at `f6d9b4e`.**
  - `Newton.solve (fun _ -> [|-1e298|]) (fun _ -> [|[|1.|]|]) [|max_float|] = Ok [|infinity|]`. The convergence branch returned `x + dx` unchecked. [`Stepper.fixed`](../../src/stepper.ml) could thus return `Ok` with a non-finite state if this happened on its last step. Extreme magnitudes only, but it contradicted the `.mli`'s implied contract. Fixed in `052fdc7`: the convergence branch returns `Error Nan` when `x + dx` is not finite, so `Newton.solve … = Ok x ⇒ Vec.finite x`; corpus line 24 pins it (it printed `Ok [inf]` before the fix).
  - `Stepper.fixed` accumulates `t` by repeated `+. h`. The last stage time is 0.99999999999999989 for `t_end = 1, dt = 0.1`, 1.0000000000000007 for `dt = 0.01`, and 2.9999999999999996 for `t_end = 3, dt = 0.3`. The code is unchanged; the `.mli` no longer says the last step "lands exactly" on `t_end` (finding 3).

### T11. The error estimate is first order

**Math.** For smooth y:

- the backward Euler local error (exact minus numerical) is −(h²/2)y″ + O(h³);
- the BDF2 local error (same convention) is −(β/6)h²(h + h_prev)y‴ + O(h⁴) (from T1's cubic residual);

so `est = y_BDF2 − y_BE` = e_BE − e_BDF2 = −(h²/2)y″ + O(h³). The controller therefore keeps *backward Euler's* local error below tol·(1 + |y|), which forces steps of order √tol rather than tol^(1/3). This quantifies the known "conservative estimate" gap and would validate a replacement estimator. Proving it needs 1-D Taylor with remainder (Coquelicot `Taylor_Lagrange`, or MathComp-Analysis `MVT` twice).

### T1f. Rounding of the computed coefficients

The code ([`Bdf2.coeffs`](../../src/bdf2.ml)) computes a₁ and a₀ with six roundings between them. Sampling 10⁵ values of ω ∈ (0, 2]:

- `a1 +. a0 <> 1.` for 50,753 of them;
- the worst gap was 5.55e-16 = 2.5 ulp(1).

Gappa can prove such a bound for an interval of ω, but not out of the box: plain interval arithmetic cannot see that a1 + a0 − 1 cancels (it gives an interval of order 1, not one of the size of the rounding errors), so the goal needs the exact counterparts of a1 and a0 and a rewriting hint. A hand analysis of the six roundings gives about 10u ≈ 5 ulp(1) as a worst case (u = 2⁻⁵³), above the sampled 2.5. Combined with T2's perturbed form, it bounds the drift of a constant solution under y′ = 0 by about 5·n·(few ulps). This is a cheap first floating-point theorem once the float toolchain exists.

```
# Gappa sketch, not checked.
@rnd = float<ieee_64, ne>;
w = rnd(w_);
# a1, a0: the expressions of Bdf2.coeffs with one rnd(...) per operation, in the same order (its shared
#   denominator is fused on arm64, exact either way)
# A1, A0: the same expressions without rounding, to be related to a1 and a0 by a rewriting hint (hint syntax not checked)
{ w in [1b-30, 2] -> a1 + a0 - 1 in ? }   # ask for the bound; a bound of 1b-50 (8u) is not guaranteed
```

### T7. Forward-difference Jacobian in binary64

**Statement.** For entry (i, j) let:

- h̃ⱼ = RN(ypⱼ − yⱼ), with ypⱼ = RN(yⱼ + RN(1e-8)·RN(1 + |yⱼ|)). On arm64 this is one fused rounding.
- εᵢ, εᵢ′ be the absolute evaluation errors of fᵢ at y and at y + h̃ⱼeⱼ.
- |∂ⱼ²fᵢ| ≤ M on the segment.

Then:

|J̃ᵢⱼ − ∂ⱼfᵢ(y)| ≤ (M/2)h̃ⱼ + (εᵢ + εᵢ′)/h̃ⱼ + (3u + O(u²))|J̃ᵢⱼ|.

And h̃ⱼ is the exact difference whenever |yⱼ| ≥ about 2.1e-8, by Sterbenz, whether fused or not. Otherwise its relative error is at most u.

**What it buys.** It makes the corpus comment's "rounding floor eps |λ₃y₃|/h ≈ 3e-5" a theorem, and justifies or improves the choice of the step. Libraries: Flocq (`Bplus_correct`, `Bfma_correct`, relative-error lemmas) plus 1-D Taylor. VCFloat can automate the per-expression rounding for a fixed f.

### T8. Damped Newton, local convergence

**Statement.** Use hypotheses that avoid multivariate calculus, which is MathComp-Analysis's thin spot:

- an affine-invariant remainder bound ‖G(x) − G(y) − G′(y)(x − y)‖ ≤ (L/2)‖x − y‖²;
- ‖J̃ − G′‖ ≤ η;
- ‖J̃⁻¹‖ ≤ β, with βη < 1/2.

Then near a root:

- the full step λ = 1 passes the Armijo test `‖G(x + dx)‖∞ ≤ (1 − 1e-4)‖G(x)‖∞`;
- iterates converge linearly at a rate of about βη + O(‖e‖) (βη/(1 − βη) if β bounds ‖(G′)⁻¹‖ of the exact Jacobian instead of ‖J̃⁻¹‖);
- the stopping test |dx|∞ ≤ 1e-10(1 + |x|∞) implies a computable residual bound.

Paşca's Kantorovich proof is prior art but predates MathComp 2. Global convergence of damped Newton is false in general (the corpus's `x0 = 0` case diverges by design), so only local statements are on offer.

### T6b, T12, and the impractical end-to-end theorem

- **T6b.** The Wilkinson/Higham backward-error theorem for GEPP is textbook, but large in Rocq. LAProof gives dot-product and FMA building blocks, but no elimination. The FMA-versus-not semantics differs per architecture, though [`Linalg`](../../src/linalg.ml) happens to be unfused on arm64. The growth factor ρₙ ≤ 2ⁿ⁻¹ makes the bound vacuous for large n anyway.
- **T12.** Formalizing the discrete-orthogonal-convolution-kernel stability proofs (Liao–Zhang, Li–Liao) is a research project. A cheaper alternative exists only if the controller's growth factor were capped at 1.9, inside Emmrich's 1.9104 energy bound. That changes bits and pinned counts, so it is a decision, not a recommendation.
- **End-to-end.** "The returned y is within C·tol of the exact solution" would need all of T4–T12 in binary64. It would also need an error estimator that *bounds* the error, which this one does not (it is a heuristic gap between two methods), plus stiffness assumptions. Not realistic.

---

## 4. Pilot plan

### 4.1 What and why

**Pilot = T1, then T2 + T3.** The deliverable is "every BDF2 step sequence vstiff produces is zero-stable with contraction q(2)", plus a mirror of [`Bdf2.coeffs`](../../src/bdf2.ml), [`Halving`](../../src/halving.ml) and [`Adaptive`](../../src/adaptive.ml)'s attempt logic that later theorems (T4, T9) build on. Why this pair:

- real arithmetic only, so it is immune to FMA and bit-level refactors;
- MathComp only;
- one theorem spans three modules, which tests can only sample;
- the result turns into a design rule for future controllers.

### 4.2 Prerequisites and cheap test pins

- **The T4 fix and its pins are in** (`dc3bb67`). [`Adaptive.integrate`](../../src/adaptive.ml) snaps steps to the floats and rejects steps that cannot move `t`. Corpus lines 25 to 27 pin one ulp at `t = 1` (using the existing [`Guard.budget`](../../test/guard.ml)), a `dt_max` below the resolution of `t`, and the `y' = 1` wrong-clock case, which a hang detector cannot see. The pins landed together with the fix, as the rule never to rewrite a passing expect line ([CONTRIBUTING.md](../../CONTRIBUTING.md), "Expect files are the contract") requires: before the fix the first two printed `no answer within 5e6 rhs calls`.
- **The Newton overflow is fixed** (`052fdc7`, corpus line 24), and the `.mli` of [`Stepper.fixed`](../../src/stepper.ml) no longer claims an exact landing on `t_end`, so the wording question raised by T10 is closed.
- **Optionally, add theorem-shadow tests in the corpus style:**
  - ω ≤ 2 on the adaptive cases, by wrapping [`Bdf2`](../../src/bdf2.ml) (a prototype did so with a recording wrapper and a global `ref`; in this repo's style the wrapper is a concrete `Bdf2` wrapper in `Guard`, the one test module allowed to raise, failing with `Invalid_argument` when h exceeds 2 h_prev by more than the snapping of T3 allows, a final step within the clock's resolution aside, which `Guard.run` already prints);
  - the Robertson invariant after every accepted step;
  - `|a1 +. a0 -. 1.| ≤ 3 ulp` over a grid.

  These test the real OCaml, not the mirror.

### 4.3 Toolchain

Not run. The dry-run solve of this set succeeded (1.2), and a second dry run also solved the 4.14.2 variant used below.

```sh
# Separate switch: MathComp needs dune < 3.24; the dev switch has dune 3.24.2.
# Compiler: see the first bullet below (5.5.0 and 4.14.2 both solve).
opam switch create vstiff-proofs ocaml-base-compiler.4.14.2 \
  --repos=default,rocq-released=https://rocq-prover.org/opam/released
opam install --switch=vstiff-proofs dune.3.23.1 rocq-core.9.2.0 rocq-stdlib.9.2.0 \
  rocq-mathcomp-ssreflect.2.6.0 rocq-mathcomp-algebra.2.6.0
# Phase 2 (analysis, floats):
opam install --switch=vstiff-proofs rocq-mathcomp-analysis.1.18.0 coq-flocq.4.2.2 \
  coq-coquelicot.3.4.5 coq-interval.4.11.5 coq-gappa.1.11.0 gappa.1.6.0   # gappa: brew gmp mpfr boost
# Only if T6b/T7 go that way (pulls CompCert 3.18 and VST 2.17):
opam install --switch=vstiff-proofs coq-vcfloat.2.4.2 rocq-laproof.2.0.1
```

- **Which OCaml for the proof switch.** Choosing 5.5.0 would let one switch parse 5.5 syntax and build vstiff, but that reason does not hold: the tripwire lives in the main workspace and runs in the dev switch or an OCaml 5.5 CI job (4.4), and Rocq-side differential checks run inside Rocq. Against 5.5.0: Rocq's `INSTALL.md` (V9.3.0) calls OCaml 5.x experimental and lists no test of 5.5; the `rocq/rocq-prover` Docker images are tagged `ocaml-4.14.2` (the MathComp images are built on them, unchecked); and on arm64 `native_compute` with OCaml 5 exists only from Rocq 9.3. So the command above uses 4.14.2 (dry runs: pilot set 51 packages, full set 86). 5.5.0 (56 and 91 packages) also solves if one switch matters more than the risk; 5.4.1 is the in-between. Do not add `coq-mathcomp-algebra-tactics`: its newest release requires MathComp < 2.6.
- **Why Rocq 9.2.0 and not 9.3.0.** 9.2 keeps VCFloat, LAProof and taylor_rocqs available. MathComp 2.6.0's changelog lists Rocq 9.0–9.2; 9.3 is allowed by opam and has a Docker image, but whether it builds was not checked.
- Editor: `vsrocq-language-server` 2.5.0 declares support for Rocq 9.0–9.3. Not solver-checked.

### 4.4 Layout

```
vstiff/
  dune                 ; new, one line: (data_only_dirs verif)   — keeps the main build Rocq-free
  dune-project src/ test/ ...
  sync/                ; in the MAIN workspace: the tripwire reads ../src, which a separate --root cannot see
    dune               ; executable canon (compiler-libs), and per mirrored module:
                       ;   (rule (targets halving.ml.out)
                       ;    (action (with-stdout-to halving.ml.out (run ./canon.exe %{dep:../src/halving.ml}))))
                       ;   (rule (alias runtest) (action (diff halving.ml.canon halving.ml.out)))
    canon.ml           ; section 2.1
    bdf2.ml.canon halving.ml.canon clock.ml.canon adaptive.ml.canon
  verif/               ; Rocq only, invisible to the main build
    dune-project       ; (lang dune 3.23) (using rocq 0.13)
    theories/
      dune             ; (rocq.theory (name Vstiff)
                       ;   (theories mathcomp.boot mathcomp.order mathcomp.ssreflect mathcomp.algebra))
      Bdf2Coeffs.v     ; T1 + mirror of Bdf2.coeffs
      ZeroStability.v  ; T2
      Controller.v     ; mirror of Halving + Adaptive attempt loop; T3
    README.md          ; which OCaml function each Rocq definition mirrors
```

- [`clock.ml`](../../src/clock.ml) is mirrored too: the floor of [`Halving`](../../src/halving.ml) is [`Clock.resolution`](../../src/clock.ml), so a change there changes the controller.
- Why `sync/` and not `verif/mirror` with `--root verif/mirror`: dune refuses dependencies outside the workspace root (`path outside the workspace: ../../src/halving.ml`), and a literal `../../src` path in an action does not exist under `_build`. Checked with dune 3.24.2 in a throwaway project: with the rules in the main workspace, comment-only edits pass, a code edit prints a diff, and `dune promote` accepts it; `(data_only_dirs verif)` at the root also hides a nested `verif/dune-project` and a `verif/dune` that the main project's dune language cannot interpret.
- The `theories` names follow MathComp 2.6's logical paths (`mathcomp.fingroup` is now `mathcomp.finite_group`). dune infers installed theory names from the `user-contrib` layout and never adds `Stdlib` implicitly (its docs), so list `Stdlib` if a file requires it. Whether it also wants the transitive closure (`mathcomp.finite_group`, `elpi`, `micromega_plugin`) listed is unchecked; confirm against `$(rocq c -where)/user-contrib`.
- `rocq makefile` with a `_RocqProject` is a fallback if dune's Rocq support gets in the way.
- After dune ≥ 3.24 becomes usable with MathComp, move to `(using rocq 0.14)`.

### 4.5 Keeping proofs and code in sync

1. **Mirror discipline.** Each Rocq definition names its OCaml counterpart, and keeps the OCaml's operation order where an IEEE instance will reuse it.
2. **Tripwire.** `dune build @runtest` at the repository root runs the `sync/` diffs of comment-free canonical sources (give them their own alias if they should not gate unrelated work). A code change to a mirrored module fails until someone re-checks the mirror and promotes. Comment edits never trigger it.
3. **Theorem-shadow tests.** These live in the normal corpus (4.2) and check the theorem's *conclusion* on the real binary.
4. **Optional, later: differential pins.** The `PrimFloat` instance of the mirror computes, inside Rocq, the same outputs the corpus prints (e.g. `coeffs` in hex for chosen ω, or [`Halving`](../../src/halving.ml) decision traces). dune Rocq expected-output tests pin those outputs, and CI diffs them against the OCaml lines.

### 4.6 CI

Two jobs. Sketch only; the repository has no CI configuration yet.

- **Proofs.** The container `mathcomp/mathcomp:2.6.0-rocq-prover-9.2` exists on Docker Hub (also `-rocq-prover-9.3`). Run `dune build --root verif`, or use [rocq-prover/docker-opam-action](https://github.com/rocq-prover/docker-opam-action) (formerly `coq-community/docker-coq-action`; the old names redirect). The image's OCaml (4.14.2) does not matter for pure Rocq. `(using rocq 0.13)` needs dune ≥ 3.23 in the image; its dune version was not checked.
- **Tripwire.** Use `ocaml/setup-ocaml` with 5.5.0 and run `dune build @runtest` in the main workspace (the same job as the corpus). Takes minutes.
- **Main build.** The Rocq side is invisible to it (`data_only_dirs`); it gains only `sync/`, which needs `compiler-libs` from the compiler itself.

### 4.7 Effort and milestones

| Step | Expert | Newcomer |
|---|---|---|
| Switch, layout, CI, tripwire | 1 d | 2 d |
| MathComp ramp-up: ssreflect tactics, `realFieldType`, bigops, `field` and `lra` | — | 5–10 d |
| T1 (+ cubic residual) | 0.5 d | 1–2 d |
| T2 (+ sharpness lemma) | 1–2 d | 3–5 d |
| T3 (mirror + invariant; real instance) | 2–4 d | 5–10 d |
| **Total pilot** | **~1–2 weeks** | **~3–6 weeks** |

Row sums: expert 4.5–7.5 days, newcomer 16–29 days. These are judgement, not data; the glue that turns a [`Halving`](../../src/halving.ml) trace into the ratio sequence T2 needs is the part most likely to run over.

**Exit criterion:** T1–T3 checked in CI, the tripwire green, and a one-page note mapping each theorem to the OCaml lines it covers.

### 4.8 After the pilot

- **Phase 2, real arithmetic:**
  - T4, now that the fix is in (it brings in Flocq for the float order);
  - T6a and the [`Jac`](../../src/jac.ml) affine-exactness lemma (it pins orientation; the corpus has a regression pin for that);
  - T9;
  - T5.
- **Phase 3, floats:**
  - T3 refinement to binary64;
  - T10;
  - T1f with Gappa;
  - T7.
- **Phase 4, only with a clear payoff:** T11 (before changing the estimator), T8, T6b, T12.
- **When variable order arrives:** generalize T1 to BDFk coefficients, and T2 to the corresponding ratio conditions. BDF3 zero-stability on variable grids is harder and the known ratio limits are small: about 1.476 (Calvo–Grande–Grigorieff, 1990) and 1.501 (Guglielmi–Zennaro, the best known), with (1 + √5)/2 ≈ 1.618 for constant ratios, per the introduction of the Li–Liao preprint. Li–Liao's 2.553 is a different, mesh-dependent stability notion for BDF3, not zero-stability.

---

## 5. Risks and limits

- **The mirror gap.** Theorems are about Rocq definitions, not the compiled OCaml. The tripwire catches drift; differential and shadow tests catch some slips. Nothing short of a verified translation from OCaml with real float semantics closes the gap, and no such tool exists today (2.3).
- **The real-versus-float gap.** Several statements true in ℝ fail in binary64:
  - "lands exactly on `t_end`";
  - a₁ + a₀ = 1;
  - `None` ⇔ singular;
  - "`Ok` means finite" for Newton.

  Termination goes the other way: true in binary64 since the clock fix, but needs liveness in ℝ. Every real theorem must be read with this in mind.
- **Platform dependence.** FMA contraction on arm64, power, riscv and s390x (not amd64) makes binary64 theorems architecture-specific, and flag-specific (`-unsafe`, inlining), unless the code pins the operation shape. Pinning it with `Sys.opaque_identity` or let-bound products changes the code. The corpus outputs did not change in the experiment, but that is not guaranteed for future problems.
- **User code.** The right-hand side is arbitrary OCaml. Every theorem takes f's behaviour as a hypothesis: c·f ≡ 0, Lipschitz bounds, evaluation-error bounds.
- **The heuristic estimator.** No theorem can say "error ≤ tol". The estimator controls backward Euler's local error (T11), not BDF2's global error.
- **Toolchain churn.** In 2026 alone: the Coq-to-Rocq rename, four dune Rocq-language versions, the removal of `(using coq)`, the `rocq-elpi` dune bound, VCFloat stuck below Rocq 9.3, and taylor_rocqs on 9.2 only. Expect a migration roughly every year. Pin versions in the switch and the Docker tag.
- **Maintenance.** Proofs at the binary64 level freeze operation order, which collides with the bit-identical-refactor workflow in both directions. Real-level proofs do not. Keep IEEE theorems few and about stable kernels.
- **Licensing.** vstiff is GPL-3.0-only. Proof-only dependencies are fine as long as they are not distributed with the library: CompCert via VCFloat, under its non-commercial terms; VST. Extracted code linking `rocq-primitive` (LGPL-2.1) would be a library dependency; check compatibility before Option B.

### Assurance bought at each level

| Level | What | Catches | Misses | Cost |
|---|---|---|---|---|
| L0 | Current expect corpus and soak | Regressions on 27 cases, bit-level reproducibility | Anything off the sampled paths (finding 1 was off them; lines 25 to 27 pin it now) | Done |
| L1 | Theorem-shadow property tests on the real binary | Violations of a stated property (ω ≤ 2, invariant drift, coefficient rounding) on more inputs | Unsampled inputs | Days |
| L2 | Real-arithmetic theorems on the mirror (pilot) | Design errors: wrong formula, unsafe ratio policy, broken invariant, non-termination in the model | Coding slips outside the mirror, all rounding effects | Weeks |
| L3 | Binary64 theorems on the mirror (Flocq, Gappa) | Special values, overflow, rounding bounds, float-only non-termination | Mirror gap, other architectures unless modelled | Weeks to months per theorem |
| L4 | Verified extraction or translation of the shipped kernels | Closes the mirror gap for extracted parts | FMA and backend semantics; collides with bit-identity and the OCaml API style | Months, plus generated code and dependencies |

**Recommendation:** L1 now; L2 as the pilot; selective L3 afterwards (T10, T1f, then T4 and T7). Skip L4 unless a kernel is rewritten anyway.

---

## Appendix A. Provenance of the measurements

The measurements in this document come from throwaway prototype scripts that are not kept in the repository, so the figures are not reproducible from the repository alone.

## Appendix B. Sources

- Rocq:
  - [releases](https://github.com/rocq-prover/rocq/releases);
  - [9.3 changelog](https://rocq-prover.org/doc/v9.3/refman/changes.html);
  - [INSTALL.md at V9.3.0](https://github.com/rocq-prover/rocq/blob/V9.3.0/INSTALL.md);
  - [using opam](https://rocq-prover.org/docs/using-opam);
  - [primitive float extraction, `ExtrOCamlFloats.v`](https://github.com/rocq-prover/stdlib/blob/V9.2.0/theories/extraction/ExtrOCamlFloats.v);
  - [`kernel/float64_63.ml`](https://github.com/rocq-prover/rocq/blob/V9.3.0/kernel/float64_63.ml).
- opam:
  - [opam-repository `rocq-runtime`](https://github.com/ocaml/opam-repository/tree/master/packages/rocq-runtime);
  - [Rocq opam archive](https://github.com/rocq-prover/opam).
- dune:
  - [CHANGES.md](https://github.com/ocaml/dune/blob/main/CHANGES.md);
  - [Rocq build language](https://dune.readthedocs.io/en/latest/rocq.html).
- OCaml:
  - [native-code compatibility notes, FMA (manual source)](https://github.com/ocaml/ocaml/blob/5.5.0/manual/src/cmds/native.etex);
  - [arm64 instruction selection](https://github.com/ocaml/ocaml/blob/5.5.0/asmcomp/arm64/selection.ml);
  - [discuss thread on FMA contraction](https://discuss.ocaml.org/t/x-y-z-is-systematically-compiled-into-an-fma-on-arm/9020).
- MathComp:
  - [MathComp 2.6.0](https://github.com/math-comp/math-comp/releases/tag/mathcomp-2.6.0) and [`algebra/matrix.v` (`cormen_lup`)](https://github.com/math-comp/math-comp/blob/mathcomp-2.6.0/algebra/matrix.v);
  - [MathComp-Analysis 1.18.0](https://github.com/math-comp/analysis/releases/tag/1.18.0).
- Floating-point libraries:
  - [Flocq](https://flocq.gitlabpages.inria.fr/) and [`PrimFloat.v`](https://gitlab.inria.fr/flocq/flocq/-/blob/flocq-4.2.2/src/IEEE754/PrimFloat.v);
  - [Coquelicot](https://coquelicot.gitlabpages.inria.fr/);
  - [CoqInterval](https://coqinterval.gitlabpages.inria.fr/);
  - [Gappa](https://gappa.gitlabpages.inria.fr/);
  - [VCFloat](https://verinum.org/vcfloat/);
  - [LAProof](https://github.com/VeriNum/LAProof);
  - [ValidSDP](https://github.com/validsdp/validsdp).
- Bridging tools:
  - [rocq-of-ocaml](https://github.com/formal-land/rocq-of-ocaml);
  - [CFML](https://github.com/charguer/cfml);
  - [Osiris/Horus, ICFP 2025](https://iris-project.org/pdfs/2025-icfp-osiris.pdf);
  - [rocq-verified-extraction](https://github.com/MetaRocq/rocq-verified-extraction) and [paper](https://doi.org/10.1145/3656379);
  - [Gospel](https://github.com/ocaml-gospel/gospel), [Ortac](https://github.com/ocaml-gospel/ortac), [Cameleer](https://github.com/ocaml-gospel/cameleer).
- Prior work:
  - [wave equation, method error](https://arxiv.org/abs/1005.0824) and [Trusting computations](https://arxiv.org/abs/1212.6641);
  - [Runge–Kutta round-off in Coq (JAR 2023)](https://link.springer.com/article/10.1007/s10817-023-09686-y);
  - [verified leapfrog (NSV 2022)](https://link.springer.com/chapter/10.1007/978-3-031-21222-2_9);
  - [VCFloat2 (CPP 2024)](https://dl.acm.org/doi/10.1145/3636501.3636953);
  - [Jacobi (CICM 2023)](https://link.springer.com/chapter/10.1007/978-3-031-42753-4_14);
  - [Lax equivalence (NFM 2021)](https://arxiv.org/abs/2103.13534);
  - [Paşca, Newton](http://www-sop.inria.fr/marelle/Ioana.Pasca/code/);
  - [Picard in Coq (ITP 2013)](https://doi.org/10.1007/978-3-642-39634-2_34);
  - [Taylor models in Coq (ITP 2024)](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ITP.2024.30);
  - [taylor_rocqs (CPP 2026)](https://doi.org/10.1145/3779031.3779097);
  - [Immler, verified ODE solver (JAR 2018)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6044317/);
  - [Isabelle numerical methods (2025)](https://arxiv.org/abs/2511.20550).
- BDF2 theory:
  - [Akrivis et al., BIT 2024](https://www.cs.uoi.gr/~akrivis/ACHYZ.pdf), which summarizes Grigorieff 1983, Crouzeix–Lisbona 1984, Becker, Emmrich and Hairer–Nørsett–Wanner;
  - [Liao–Zhang, arXiv 1912.11182](https://arxiv.org/abs/1912.11182);
  - [Li–Liao, SINUM 2022](https://doi.org/10.1137/21M1462398).
