# Mutation baseline

A run of [mutate](README.md) on the kernel and the solver: `src/numerics/{vec,linalg,newton,jac}.ml` and `src/{bdf1,bdf2,stage,halving,adaptive,stepper,clock,check}.ml`, against `corpus` and `soak`. No survivor was fixed. From the repository root, with dune on the `PATH`:

```sh
./_build/default/tools/mutate/mutate.exe run --results results.tsv \
  src/numerics/{vec,linalg,newton,jac}.ml src/{bdf1,bdf2,stage,halving,adaptive,stepper,clock,check}.ml
```

**Total wall time: 22 min 05 s**, for 450 mutants on 4 workers (an arm64 Mac with 10 cores, OCaml 5.5.0, dune 3.24.2). The unmodified tests took 8.5 s, so a mutant had 25.5 s. The machine was shared with other work during the run.

The run killed 52 mutants by the timeout. A second run of just those, on an idle machine (6 min 13 s, a limit of 27.1 s), timed out 50 of them again; `src/stage.ml:18:59:float_div2` and `src/halving.ml:28:71:fdiv_to_fmul` ended with a diff a little after the first limit. The table uses the second outcome for these 52 and the first for the other 398.

All 9 stillborn mutants swap the body of a match arm for one that names a variable the target arm does not bind (`Unbound value`).

Every survivor has a verdict line:

- **Equivalent**: no input tells the mutant from the original; the line gives the argument.
- **Looks equivalent**: it differs only by rounding, or where two computed floats are exactly equal, which no natural test reaches; the line says what it takes.
- **Not equivalent**: an input tells it from the original, and the line gives it. The numbers were measured by applying the mutant (`mutate apply ID`) to a scratch copy of the library and running that input.

| Module | Mutants | Killed | Survived | Stillborn |
| --- | ---: | ---: | ---: | ---: |
| src/adaptive.ml | 61 | 43 | 14 | 4 |
| src/bdf1.ml | 1 | 1 | 0 | 0 |
| src/bdf2.ml | 47 | 45 | 0 | 2 |
| src/check.ml | 37 | 34 | 3 | 0 |
| src/clock.ml | 7 | 7 | 0 | 0 |
| src/halving.ml | 66 | 63 | 3 | 0 |
| src/numerics/jac.ml | 25 | 25 | 0 | 0 |
| src/numerics/linalg.ml | 63 | 59 | 3 | 1 |
| src/numerics/newton.ml | 96 | 63 | 31 | 2 |
| src/numerics/vec.ml | 7 | 7 | 0 | 0 |
| src/stage.ml | 10 | 10 | 0 | 0 |
| src/stepper.ml | 30 | 26 | 4 | 0 |
| **Total** | 450 | 383 | 58 | 9 |

Of the 383 killed mutants, 50 were killed by the timeout.

## Survivors

### src/adaptive.ml

- `src/adaptive.ml:11:98:int_minus1`: `50` -> `49`
  Not equivalent: the default `max_rejects` is unpinned: a run whose steps all fail (rhs nan after t0, `dt0 = dt_max = 1e3`) ends `Error StepRejected 51`, and `Error StepRejected 50` with this default.
- `src/adaptive.ml:11:98:int_plus1`: `50` -> `51`
  Not equivalent: the default `max_rejects` is unpinned: the run of the previous entry ends `Error StepRejected 52` with this default instead of 51.
- `src/adaptive.ml:13:22:fsub_to_fadd`: `p.t_end -. p.t0` -> `p.t_end +. p.t0`
  Not equivalent: the default steps come from `t_end + t0` instead of the span, and no corpus case with t0 other than 0 depends on them (each gives both steps or ends at once): `y' = 0` on [5, 6] takes 45 accepted steps instead of 58.
- `src/adaptive.ml:17:61:fdiv_to_fmul`: `span /. 10.` -> `span *. 10.`
  Not equivalent: the default `dt_max = span / 10` is a known gap (docs/numerics/06-the-corpus.md, section 10): no corpus step grows that far; with a cap of ten spans `y' = 0` on [0, 1] takes 55 accepted steps instead of 58.
- `src/adaptive.ml:17:64:float_div10`: `10.` -> `1.`
  Not equivalent: the default `dt_max` gap: a cap of one span instead of a tenth moves no corpus line; `y' = 0` on [0, 1] takes 55 accepted steps instead of 58.
- `src/adaptive.ml:17:64:float_div2`: `10.` -> `5.`
  Not equivalent: the default `dt_max` gap: a cap of a fifth of the span moves no corpus line; `y' = 0` on [0, 1] takes 56 accepted steps instead of 58.
- `src/adaptive.ml:17:64:float_zero`: `10.` -> `0.`
  Not equivalent: the default `dt_max` gap: `span / 0.` is infinity, no cap at all, and moves no corpus line; `y' = 0` on [0, 1] takes 55 accepted steps instead of 58.
- `src/adaptive.ml:33:21:ge_to_gt`: `dt >= remaining` -> `dt > remaining`
  Not equivalent: when dt equals the remaining span the step is no longer the last, so t advances by `t + dt`, which can round away from t_end: with `dt0 = dt_max = t_end - t0` from 9.58 to 641.005 the run ends one ulp from t_end instead of at it, against the contract of `Adaptive.solution`.
- `src/adaptive.ml:33:47:le_to_lt`: `remaining <= (Clock.resolution at.t)` -> `remaining < (Clock.resolution at.t)`
  Not equivalent: the remainder rule needs `<=`: a span of exactly `Clock.resolution t0` (from 1 to 1 plus 16 epsilon) ends `Error StepRejected 1` with this, and `Ok` with the original.
- `src/adaptive.ml:37:12:le_to_lt`: `h <= 0.` -> `h < 0.`
  Equivalent: h is 0 only for a non-final step at a t other than 0 (a final step has h = remaining > 0, and at t = 0 a step is dt > 0), and there the second operand of the `||`, `(not last) && h < Clock.resolution t`, is true for 0 as well.
- `src/adaptive.ml:37:18:or_to_and`: `(h <= 0.) || ((not last) && (h < (Cl...` -> `(h <= 0.) && ((not last) && (h < (Cl...`
  Not equivalent: this rejects a step only when $`h \le 0`$, so a step below the resolution of $`t`$ that still moves $`t`$ is taken: a known gap (docs/numerics/06-the-corpus.md, section 10), $`y' = 1`$ from $`t = 10^{15}`$ over 100 with `dt0 = dt_max = 0.19` ends `Ok` instead of `Error StepRejected 1`.
- `src/adaptive.ml:37:38:lt_to_le`: `h < (Clock.resolution at.t)` -> `h <= (Clock.resolution at.t)`
  Not equivalent: a step of exactly the resolution of t is taken by the original and rejected with `<=`: `y' = 1` from t = 1 over two resolutions with `dt0 = dt_max = Clock.resolution 1` is `Ok` (2 steps), and `Error StepRejected 1` with this.
- `src/adaptive.ml:42:11:arm_body_from_2`: `retry (Ode.Solver e)` -> `retry Ode.Too_large`
  Not equivalent: the driver tells the controller `Solver e` after a failed step and `Too_large` after a rejected estimate, and only `Halving`, which ignores the reason, is tested; a controller that records the reasons sees `Too_large` four times instead of `Solver` for a run of failed solves.
- `src/adaptive.ml:45:27:and_to_or`: `(Vec.finite p.y0) && (Vec.finite f0)` -> `(Vec.finite p.y0) || (Vec.finite f0)`
  Not equivalent: a known gap (docs/numerics/06-the-corpus.md, section 10): the start check needs both halves only when y0 is not finite and rhs t0 y0 is; with `y0 = nan` and a constant right-hand side it ends `Error StepRejected 51` instead of `Error Nan`.

### src/check.ml

- `src/check.ml:2:30:fsub_to_fadd`: `p.t_end -. p.t0` -> `p.t_end +. p.t0`
  Not equivalent: no case gives `Stepper.fixed` a t_end before t0 with dt > 0 (the `dt = 0` case never reaches the span test, and only `Adaptive` has a case for t_end < t0); with t0 = 1 and t_end = 0 the original raises `Invalid_argument` and this runs on, ending `Error Diverged`.
- `src/check.ml:3:19:gt_to_ge`: `p.t_end > p.t0` -> `p.t_end >= p.t0`
  Not equivalent: an empty span with a tiny dt at a large t, as in `Stepper.fixed ~dt:1e-7` with t0 = t_end = 1e10, is `Ok` (nothing to do) in the original and raises `dt is below the resolution of t` with this.
- `src/check.ml:3:32:lt_to_le`: `dt < (Clock.resolution (Floa...` -> `dt <= (Clock.resolution (Floa...`
  Not equivalent: a `dt` exactly equal to the resolution of the larger time is accepted by the original (backward Euler from 1e10 to 1e10 + 1 gives 0.368) and raises `dt is below the resolution of t` with `<=`; no case sits on that boundary.

### src/halving.ml

- `src/halving.ml:28:103:le_to_lt`: `...loat.abs yi))) err y)) <= c.tol` -> `...loat.abs yi))) err y)) < c.tol`
  Not equivalent: a boundary of a public function: `Halving.acceptable` with the scaled error equal to tol is `true`, and `false` with this; no case builds one.
- `src/halving.ml:33:27:int_plus1`: `1` -> `2`
  Equivalent: the streak runs 0, 2, 1 instead of 0, 1, 2 and still reaches 3 on the third accept, then restarts at 0; a rejection resets it to 0 in both and only `accepted` reads it (the proposals are the same on a mixed run of accepts and rejections).
- `src/halving.ml:50:15:gt_to_ge`: `failures > c.max_rejects` -> `failures >= c.max_rejects`
  Not equivalent: the `max_rejects` give-up is a known gap (docs/numerics/06-the-corpus.md, section 10): with `max_rejects = 3` a run of rejections ends `Error StepRejected 4`, and `Error StepRejected 3` with this; every corpus case ends by the step floor first.

### src/numerics/linalg.ml

- `src/numerics/linalg.ml:10:68:gt_to_ge`: `...oat.abs ((a.(j)).(0))) > (Float.abs ((a.(i)).(0)...` -> `...oat.abs ((a.(j)).(0))) >= (Float.abs ((a.(i)).(0)...`
  Looks equivalent: on a tie in $`\lvert a_{i0} \rvert`$ either row is a valid pivot. On 3000 random integer systems it changes the answer of 1069: the 1064 nonsingular ones by at most 1.1e-13 relative, below every printed digit, and 5 singular ones, where `linalg.mli` leaves the answer open (`None` or a huge vector).
- `src/numerics/linalg.ml:12:10:int_minus1`: `1` -> `0`
  Equivalent: `best 0 0` first compares row 0 with itself, which never replaces it, so it is `best 0 1` with one more step (bit-identical on 3000 random integer systems).
- `src/numerics/linalg.ml:15:26:int_minus1`: `0` -> `(-1)`
  Equivalent: `swapped` is only called with `i + 1`, never with 0, so its `i = 0` branch is dead and `i = -1` is as dead (bit-identical on 3000 random integer systems).

### src/numerics/newton.ml

- `src/numerics/newton.ml:6:11:float_div2`: `1e-10` -> `5e-11`
  Not equivalent: the tolerance is pinned within a decade (`1e-9` and `1e-11` are killed) but not within a factor of two. `x^3 = 0` from 0.05 prints `Ok 1.76e-10` in 47 steps, and `Ok 7.84e-11` in 49 with this tolerance.
- `src/numerics/newton.ml:6:11:float_x2`: `1e-10` -> `2e-10`
  Not equivalent: same gap as `5e-11`: `x^3 = 0` from 0.05 prints `Ok 3.97e-10` in 45 steps instead of `Ok 1.76e-10` in 47.
- `src/numerics/newton.ml:7:16:int_minus1`: `50` -> `49`
  Not equivalent: the only case that meets the limit, `x^3 = 0` from 1, needs 55 steps, so every limit from 1 to 54 ends the same. `x^3 = 0` from 0.1 needs 49 steps: `Ok` with limit 50, `Error Diverged` with 49.
- `src/numerics/newton.ml:7:16:int_plus1`: `50` -> `51`
  Not equivalent: `x^3 = 0` from 0.15 needs 50 steps: `Error Diverged` with limit 50, `Ok` with 51. The limit is pinned from above (100 is killed) but not at 51.
- `src/numerics/newton.ml:10:19:float_div10`: `1.` -> `0.1`
  Not equivalent: the damping floor is pinned from above only (the corpus's deep line search needs 1/512, none of its cases needs 1/1024). A deeper floor shows only when a search is exhausted: with a Jacobian of the wrong sign, 15 evaluations of f instead of 12.
- `src/numerics/newton.ml:10:19:float_div2`: `1.` -> `0.5`
  Not equivalent: the damping floor is pinned from above only (no corpus case needs 1/1024); with a Jacobian of the wrong sign the search takes 13 evaluations of f instead of 12.
- `src/numerics/newton.ml:10:19:float_x2`: `1.` -> `2.`
  Not equivalent: a floor of 1/512 still admits the factor 1/512 that the deep line search needs, because the test is `lambda < min_damping`. `x^2 - 2` from 5e-4 needs 1/1024: `Ok` with the original, `Error Diverged` with this floor.
- `src/numerics/newton.ml:10:19:float_zero`: `1.` -> `0.`
  Not equivalent: with no floor the search ends only when lambda is so small that x does not move, and then wastes the iterations up to the limit: with a Jacobian of the wrong sign, 5298 evaluations of f instead of 12, the same `Error Diverged`.
- `src/numerics/newton.ml:10:25:float_div2`: `1024.` -> `512.`
  Not equivalent: a floor of 1/512 still admits the factor 1/512 that the deep line search needs, because the test is `lambda < min_damping`. `x^2 - 2` from 5e-4 needs 1/1024: `Ok` with the original, `Error Diverged` with this floor.
- `src/numerics/newton.ml:10:25:float_x10`: `1024.` -> `10240.`
  Not equivalent: the damping floor is pinned from above only (no corpus case needs 1/1024); with a Jacobian of the wrong sign the search takes 15 evaluations of f instead of 12.
- `src/numerics/newton.ml:10:25:float_x2`: `1024.` -> `2048.`
  Not equivalent: the damping floor is pinned from above only (no corpus case needs 1/1024); with a Jacobian of the wrong sign the search takes 13 evaluations of f instead of 12.
- `src/numerics/newton.ml:14:14:float_div10`: `1e-4` -> `1e-5`
  Looks equivalent: a constant of `1e-5` accepts a different step only when it cuts $`|f|`$ by a fraction between `1e-5 lambda` and `1e-4 lambda`, a band that none of the corpus, `x^2 - 5`, `x^2 - 2`, `x^3` or a wrong-sign Jacobian lands in; only a crafted case with a badly scaled Jacobian would.
- `src/numerics/newton.ml:14:14:float_div2`: `1e-4` -> `5e-5`
  Looks equivalent: a constant of `5e-5` differs from `1e-4` only for a step that cuts $`|f|`$ by a fraction between the two (about 0.005% of lambda), a band no natural case lands in.
- `src/numerics/newton.ml:14:14:float_x10`: `1e-4` -> `1e-3`
  Looks equivalent: a constant of `1e-3` differs from `1e-4` only for a step that cuts $`|f|`$ by a fraction between the two, a band no natural case lands in.
- `src/numerics/newton.ml:14:14:float_x2`: `1e-4` -> `2e-4`
  Looks equivalent: a constant of `2e-4` differs from `1e-4` only for a step that cuts $`|f|`$ by a fraction between the two, a band no natural case lands in; `x` with a Jacobian of 1e4 from 1, where every step cuts $`|f|`$ by exactly `1e-4 lambda`, ends `Error Diverged` after 12 evaluations of f instead of 51.
- `src/numerics/newton.ml:14:14:float_zero`: `1e-4` -> `0.`
  Not equivalent: the Armijo constant is a known gap (docs/numerics/06-the-corpus.md, section 10). With 0 a step that leaves $`|f|`$ unchanged is accepted: `x^2 - 5` from 1 takes 6 evaluations of f instead of 7.
- `src/numerics/newton.ml:22:15:ge_to_gt`: `k >= max_iter` -> `k > max_iter`
  Not equivalent: the same boundary as the limit itself: `x^3 = 0` from 0.15 needs 50 steps and ends `Error Diverged` with the original, `Ok` with this.
- `src/numerics/newton.ml:26:9:arm_body_from_4`: `Error Fail.Diverged` -> `let rec damp lambda = if lambda < min_damping then Error Fail.Diverged e...`
  Not equivalent: for the test's f every trial point is infinite and the search ends in the same `Error Diverged`, after eleven evaluations at infinity. For an f that is finite there, such as `1 / (1 + x^2)` with a Jacobian of 1e-320, it is `Error Nan`.
- `src/numerics/newton.ml:29:38:le_to_lt`: `(Vec.norm_inf dx) <= (tol *. (1. +. (Vec.nor...` -> `(Vec.norm_inf dx) < (tol *. (1. +. (Vec.nor...`
  Looks equivalent: it differs only when $`|dx|`$ equals $`tol (1 + |x|)`$ exactly, which only a crafted case reaches: `x - 1e-10` from 0 takes 2 evaluations of f instead of 1, with the same answer.
- `src/numerics/newton.ml:29:49:float_div10`: `1.` -> `0.1`
  Not equivalent: the absolute part of the stopping test is unpinned: `x^3 = 0` from 0.05 ends `Error Diverged` after 50 steps with `0.1 + |x|`, and `Ok 1.76e-10` with `1 + |x|`.
- `src/numerics/newton.ml:29:49:float_div2`: `1.` -> `0.5`
  Not equivalent: the absolute part of the stopping test is unpinned: `x^3 = 0` from 0.05 prints `Ok 7.84e-11` with `0.5 + |x|` and `Ok 1.76e-10` with `1 + |x|`.
- `src/numerics/newton.ml:29:49:float_x2`: `1.` -> `2.`
  Not equivalent: the absolute part of the stopping test is unpinned: `x^2 - 5` from 1 takes 6 evaluations of f with `2 + |x|` instead of 7.
- `src/numerics/newton.ml:29:49:float_zero`: `1.` -> `0.`
  Not equivalent: with `0. + |x|` the stopping test is purely relative, which a root at 0 never meets: `x^3 = 0` from 0.05 ends `Error Diverged`, where the original prints `Ok 1.76e-10`.
- `src/numerics/newton.ml:36:23:lt_to_le`: `lambda < min_damping` -> `lambda <= min_damping`
  Not equivalent: `lambda <= min_damping` rejects the factor 1/1024 itself, so the depth is 1/512; `x^2 - 2` from 5e-4 needs 1/1024: `Ok` with the original, `Error Diverged` with this.
- `src/numerics/newton.ml:40:53:le_to_lt`: `(Vec.norm_inf fx') <= ((1. -. (armijo *. lamb...` -> `(Vec.norm_inf fx') < ((1. -. (armijo *. lamb...`
  Looks equivalent: it differs only when $`\lvert f(x') \rvert = (1 - \mathrm{armijo}\,\lambda)\,\lvert f(x) \rvert`$ exactly, which only a crafted case reaches: `x` with a Jacobian of 1e4 from 1 ends `Error Diverged` after 12 evaluations of f instead of 51.
- `src/numerics/newton.ml:40:60:fsub_to_fadd`: `1. -. (armijo *. lambda)` -> `1. +. (armijo *. lambda)`
  Not equivalent: the Armijo gap again: `(1 + armijo lambda) r` accepts a step that leaves $`|f|`$ unchanged, and `x^2 - 5` from 1 takes 6 evaluations of f instead of 7.
- `src/numerics/newton.ml:40:71:fmul_to_fdiv`: `armijo *. lambda` -> `armijo /. lambda`
  Not equivalent: `armijo / lambda` is the same constant at $`\lambda = 1`$ and 0.1 at $`\lambda = 1/1024`$; $`x^2 - 2`$ from $`5 \times 10^{-4}`$ needs $`\lambda = 1/1024`$, and ends `Error Diverged` instead of `Ok`.
- `src/numerics/newton.ml:40:106:int_plus1`: `1` -> `2`
  Not equivalent: counting each step twice halves the limit to 25: `x^3 = 0` from 0.05 or 0.1 ends `Error Diverged` after 26 evaluations of f, where the original prints `Ok`.
- `src/numerics/newton.ml:41:36:float_x2`: `2.` -> `4.`
  Not equivalent: the damping factor is pinned only loosely: dividing by 4 instead of 2 passes every case, and `x^2 - 2` from 5e-4 then takes 11 evaluations of f instead of 16.
- `src/numerics/newton.ml:45:11:int_minus1`: `0` -> `(-1)`
  Not equivalent: starting the step counter at -1 allows 51 steps, the same gap as the limit: `x^3 = 0` from 0.15 ends `Ok` instead of `Error Diverged`.
- `src/numerics/newton.ml:45:11:int_plus1`: `0` -> `1`
  Not equivalent: starting the step counter at 1 allows 49 steps, the same gap as the limit: `x^3 = 0` from 0.1 ends `Error Diverged` instead of `Ok`.

### src/stepper.ml

- `src/stepper.ml:11:22:fsub_to_fadd`: `p.t_end -. p.t0` -> `p.t_end +. p.t0`
  Not equivalent: every fixed-step case starts at t0 = 0 or fails an argument check first, so the span `t_end + t0` is never told from `t_end - t0`: backward Euler on [1, 2] with dt = 0.1 gives -0.070 instead of 0.386 (30 steps, the last one backwards), and an empty span at t = 1e10 with dt = 1e-7 never ends.
- `src/stepper.ml:12:27:and_to_or`: `(Vec.finite p.y0) && (Vec.finite f0)` -> `(Vec.finite p.y0) || (Vec.finite f0)`
  Not equivalent: a known gap (docs/numerics/06-the-corpus.md, section 10): with `y0 = nan` and a finite right-hand side an empty span is `Error Nan` in the original and `Ok [nan]` with this; the corpus lines have a right-hand side that is nan too.
- `src/stepper.ml:15:17:int_minus1`: `1` -> `0`
  Not equivalent: the rule of at least one step is a known gap (docs/numerics/06-the-corpus.md, section 10): with `dt` longer than the span the original takes one step (backward Euler on [0, 1], dt = 4, gives 0.5), and this takes none and returns y0 (1).
- `src/stepper.ml:15:17:int_plus1`: `1` -> `2`
  Not equivalent: `max 2` takes two steps where the original takes one if dt is longer than 2/3 of the span: backward Euler on [0, 1] with dt = 4 gives 0.444 instead of 0.5; no case has dt that long.

## What the survivors say

The 47 survivors that are not equivalent make a short list of gaps. Some are already in [docs/numerics/06-the-corpus.md](../../docs/numerics/06-the-corpus.md), section 10: the Armijo constant, the defaults and the give-up of the drivers, the start checks, the rule of at least one step and the rule for short steps. The rest are new. Each is a limit or a boundary that the corpus meets from one side only.

- **Newton's limits (25).**
  - The iteration limit is pinned from above only: 49 and 51, `>` for `>=`, a counter that starts at -1 or 1, or counts twice. `x^3 = 0` from 0.1 (49 steps) and from 0.15 (50 steps) tell them apart.
  - The damping floor is pinned from above only (8 mutants), and so is the damping factor. `x^2 - 2` from 5e-4 needs 1/1024.
  - The Armijo constant (0, `1 +`, `/ lambda`), the tolerance and the absolute part of the stopping test (6 mutants; `x^3 = 0` from 0.05), and the step that is not finite.
- **Defaults and boundaries of the adaptive driver (13).**
  - The defaults: `max_rejects` (2), `dt_max` (4), and the default steps when t0 is not 0 (1).
  - The clock rules at their boundaries: a step of exactly the resolution, a span of exactly the resolution, a step equal to the remaining span, and the rule that rejects a step below the resolution (4).
  - The reason the driver gives the controller, and one half of the start check (2).
- **The fixed-step driver and its checks (7).** Nothing runs it from a t0 other than 0, with a t_end before t0, with a step longer than the span, or with an empty span at a large t; a start check with one half; a step exactly at the resolution.
- **Halving (2).** `acceptable` at an error equal to tol, and the `max_rejects` give-up.

The seven that look equivalent are the tie-break of the pivot, the Armijo constant within a factor of ten (4), and two equalities of computed floats in Newton. The four equivalent ones are in `Linalg.pivot`, in `Linalg.swapped`, in the streak of `Halving` and in the `h <= 0.` test of `Adaptive`. `Jac`, `Stage`, `Bdf1`, `Bdf2`, `Clock` and `Vec` have no survivor.
