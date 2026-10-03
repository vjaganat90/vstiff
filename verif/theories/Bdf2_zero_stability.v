(* T2 of docs/plans/formal-verification.md: variable-step BDF2 is zero-stable
   on every grid whose step ratios stay at most some ws < 1 + sqrt 2, with an
   explicit bound in terms of the start and the perturbations, and no constant
   ratio at or above 1 + sqrt 2 is. *)
From mathcomp Require Import boot order algebra.
From Vstiff Require Import Bdf2_weights.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import Order.Theory GRing.Theory Num.Theory.
Local Open Scope ring_scope.

Section Increment.

Variable R : fieldType.

(* The parasitic root of the recurrence with the ratio frozen at w. *)
Definition q (w : R) := - a0 w.

Lemma qE (w : R) : q w = w ^+ 2 / (1 + 2 * w).
Proof. by rewrite /q /a0 /den mulNr opprK expr2. Qed.

(* Since a1 = 1 - a0, a step moves y by q times the previous move. *)
Lemma bdf2_increment (w y0 y1 : R) : den w != 0 ->
  a1 w * y1 + a0 w * y0 - y1 = q w * (y1 - y0).
Proof. by rewrite /q /a1 /a0 /den => nz; field. Qed.

End Increment.

Section Contraction.

Variable R : realFieldType.

Lemma den_gt0 (w : R) : 0 <= w -> 0 < den w.
Proof. by rewrite /den => w0; lra. Qed.

Lemma q_ge0 (w : R) : 0 <= w -> 0 <= q w.
Proof. by move=> w0; rewrite qE divr_ge0 ?sqr_ge0 //; lra. Qed.

(* q grows with the ratio, so the largest ratio is the worst case. *)
Lemma q_le (w ws : R) : 0 <= w -> w <= ws -> q w <= q ws.
Proof.
move=> w0 wws; have d0 : 0 < 1 + 2 * w by lra.
have d1 : 0 < 1 + 2 * ws by lra.
rewrite !qE ler_pdivrMr // mulrAC ler_pdivlMr //.
have : 0 <= w * ws * (ws - w) by rewrite !mulr_ge0 ?subr_ge0 //; lra.
by nra.
Qed.

(* With every ratio in [0, ws], a step difference is at most q ws times the
   one before, plus the perturbation. *)
Lemma bdf2_step_difference (ws : R) (w y d : nat -> R) :
  (forall n, 0 <= w n <= ws) ->
  (forall n, y n.+2 = a1 (w n) * y n.+1 + a0 (w n) * y n + d n) ->
  forall n, `|y n.+2 - y n.+1| <= q ws * `|y n.+1 - y n| + `|d n|.
Proof.
move=> wb run n; have /andP [w0 wws] := wb n.
rewrite run addrAC bdf2_increment ?lt0r_neq0 ?den_gt0 //.
apply: (le_trans (ler_normD _ _)); rewrite lerD2r normrM ger0_norm ?q_ge0 //.
by rewrite ler_wpM2r ?normr_ge0 ?q_le.
Qed.

(* The bound of T2: the start and the perturbations, amplified by at most
   1 / (1 - q ws), whatever the grid. *)
Theorem bdf2_stable (ws : R) (w y d : nat -> R) :
  q ws < 1 -> (forall n, 0 <= w n <= ws) ->
  (forall n, y n.+2 = a1 (w n) * y n.+1 + a0 (w n) * y n + d n) ->
  forall n,
    `|y n| <= `|y 0| + (`|y 1 - y 0| + \sum_(k < n) `|d k|) / (1 - q ws).
Proof.
move=> q1 wb run.
have q0 : 0 <= q ws by have /andP [w0 wws] := wb 0%N; apply: q_ge0; lra.
pose S n := \sum_(k < n) `|y k.+1 - y k|.
pose D n := \sum_(k < n) `|d k|.
have S_le n : (1 - q ws) * S n <= `|y 1 - y 0| + D n.
  case: n => [|n]; first by rewrite /S /D !big_ord0 mulr0 addr0 normr_ge0.
  have SS : S n <= S n.+1 by rewrite /S big_ord_recr /= lerDl normr_ge0.
  have DD : D n <= D n.+1 by rewrite /D big_ord_recr /= lerDl normr_ge0.
  have rec : S n.+1 <= `|y 1 - y 0| + q ws * S n + D n.
    rewrite /S /D big_ord_recl /= -addrA lerD2l mulr_sumr -big_split /=.
    by apply: ler_sum => k _; apply: bdf2_step_difference.
  by have := ler_wpM2l q0 SS; nra.
have y_le n : `|y n| <= `|y 0| + S n.
  elim: n => [|n ih]; first by rewrite /S big_ord0 addr0.
  have := lerB_dist (y n.+1) (y n); rewrite lerBlDl => h.
  apply: (le_trans h); rewrite /S big_ord_recr /= addrA lerD2r.
  exact: ih.
move=> n; apply: (le_trans (y_le n)).
rewrite lerD2l ler_pdivlMr ?subr_gt0 // mulrC.
exact: S_le.
Qed.

(* Without perturbations the step differences contract geometrically... *)
Theorem bdf2_contraction (ws : R) (w y : nat -> R) :
  (forall n, 0 <= w n <= ws) ->
  (forall n, y n.+2 = a1 (w n) * y n.+1 + a0 (w n) * y n) ->
  forall n, `|y n.+1 - y n| <= q ws ^+ n * `|y 1 - y 0|.
Proof.
move=> wb run; have q0 : 0 <= q ws.
  by have /andP [w0 wws] := wb 0%N; apply: q_ge0; lra.
have run0 n : y n.+2 = a1 (w n) * y n.+1 + a0 (w n) * y n + 0.
  by rewrite addr0.
elim=> [|n ih]; first by rewrite expr0 mul1r.
have := bdf2_step_difference wb run0 n; rewrite normr0 addr0 => h.
by apply: (le_trans h); rewrite exprS -mulrA ler_wpM2l.
Qed.

(* ...and the solution stays within its start, amplified by 1 / (1 - q ws). *)
Theorem bdf2_bounded (ws : R) (w y : nat -> R) :
  q ws < 1 -> (forall n, 0 <= w n <= ws) ->
  (forall n, y n.+2 = a1 (w n) * y n.+1 + a0 (w n) * y n) ->
  forall n, `|y n| <= `|y 0| + `|y 1 - y 0| / (1 - q ws).
Proof.
move=> q1 wb run n.
have run0 k : y k.+2 = a1 (w k) * y k.+1 + a0 (w k) * y k + 0.
  by rewrite addr0.
have zero : \sum_(k < n) `|0 : R| = 0 by apply: big1 => k _; exact: normr0.
by have := bdf2_stable q1 wb run0 n; rewrite zero addr0.
Qed.

End Contraction.

Section Threshold.

Variable R : rcfType.

Lemma sqr_sqrt2 : Num.sqrt 2 ^+ 2 = 2 :> R.
Proof. by rewrite sqr_sqrtr. Qed.

Lemma sqrt2_gt1 : 1 < Num.sqrt 2 :> R.
Proof. by have := sqrtr_ge0 (2 : R); have := sqr_sqrt2; nra. Qed.

(* Below 1 + sqrt 2 the frozen recurrence contracts, at or above it not. *)
Lemma q_lt1 (w : R) : 0 <= w -> (q w < 1) = (w < 1 + Num.sqrt 2).
Proof.
move=> w0; have s2 := sqr_sqrt2; have s1 := sqrt2_gt1.
have d0 : 0 < 1 + 2 * w by lra.
rewrite qE ltr_pdivrMr // mul1r.
by apply/idP/idP => h; nra.
Qed.

Lemma q_threshold : q (1 + Num.sqrt 2) = 1 :> R.
Proof.
have s2 := sqr_sqrt2; have s1 := sqrt2_gt1.
have nz : 1 + 2 * (1 + Num.sqrt 2) != 0 :> R by apply: lt0r_neq0; lra.
by rewrite qE; field: s2.
Qed.

(* T2: every grid whose ratios lie in [0, ws], ws < 1 + sqrt 2. *)
Theorem bdf2_zero_stable (ws : R) (w y d : nat -> R) :
  ws < 1 + Num.sqrt 2 -> (forall n, 0 <= w n <= ws) ->
  (forall n, y n.+2 = a1 (w n) * y n.+1 + a0 (w n) * y n + d n) ->
  q ws < 1 /\ forall n,
    `|y n| <= `|y 0| + (`|y 1 - y 0| + \sum_(k < n) `|d k|) / (1 - q ws).
Proof.
move=> wsl wb run.
have q1 : q ws < 1 by have /andP [w0 wws] := wb 0%N; rewrite q_lt1 //; lra.
by split=> //; exact: bdf2_stable.
Qed.

(* Sharpness: at a constant ratio of 1 + sqrt 2 or more, q w >= 1, and every
   solution that moves at its first step drifts at least linearly. *)
Theorem bdf2_sharp (w : R) (y : nat -> R) : 1 + Num.sqrt 2 <= w ->
  (forall n, y n.+2 = a1 w * y n.+1 + a0 w * y n) ->
  1 <= q w /\ forall n, n%:R * `|y 1 - y 0| <= `|y n - y 0|.
Proof.
move=> wl run; have s1 := sqrt2_gt1.
have w0 : 0 <= w by lra.
have q1 : 1 <= q w by rewrite leNgt q_lt1 // -leNgt.
have e n : y n.+1 - y n = q w ^+ n * (y 1 - y 0).
  elim: n => [|n ih]; first by rewrite expr0 mul1r.
  rewrite run bdf2_increment ?lt0r_neq0 ?den_gt0 // ih.
  by rewrite exprS mulrA.
have drift n : y n - y 0 = (\sum_(k < n) q w ^+ k) * (y 1 - y 0).
  elim: n => [|n ih]; first by rewrite big_ord0 mul0r subrr.
  by rewrite big_ord_recr /= mulrDl -ih -e; ring.
have sum_ge n : n%:R <= \sum_(k < n) q w ^+ k.
  elim: n => [|n ih]; first by rewrite big_ord0 mulr0n.
  by rewrite big_ord_recr /= mulrSr lerD ?exprn_ege1.
split=> // n; have s0 : 0 <= \sum_(k < n) q w ^+ k.
  by apply: sumr_ge0 => k _; rewrite exprn_ge0 //; lra.
by rewrite (drift n) normrM (ger0_norm s0) ler_wpM2r ?normr_ge0 ?sum_ge.
Qed.

End Threshold.
