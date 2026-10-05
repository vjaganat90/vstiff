(* The weights of variable-step BDF2, as Bdf2.coeffs computes them for the
   step ratio w = h / h_prev, over any field, and T1: the formula is exact on
   quadratics, and its residual on a cubic is a fixed multiple of the
   coefficient of t^3. *)
From mathcomp Require Import boot order algebra.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory.
Local Open Scope ring_scope.

Section Bdf2.

Variable R : fieldType.

(* Each definition keeps the operation order of the OCaml it mirrors, for a
   later IEEE instance: a0 is (-.(omega *. omega)) /. d there too. *)

(* Mirrors d in Bdf2.coeffs (src/bdf2.ml). *)
Definition den (w : R) := 1 + 2 * w.
(* Mirrors the field a1 of Bdf2.coeffs. *)
Definition a1 (w : R) := (1 + w) * (1 + w) / den w.
(* Mirrors the field a0 of Bdf2.coeffs. *)
Definition a0 (w : R) := - (w * w) / den w.
(* Mirrors the field beta of Bdf2.coeffs. *)
Definition beta (w : R) := (1 + w) / den w.

(* The history weights sum to one, so the formula is exact on constants. *)
Lemma weights_sum_to_one (w : R) : den w != 0 -> a1 w + a0 w = 1.
Proof. by rewrite /a1 /a0 /den => nz; field. Qed.

(* The residual of the step from tn to tn + h on p, after a step hp, for the
   ratio w = h / hp: the value of p minus what the formula gives. It is the
   local truncation error when p is the solution. *)
Definition residual (w hp tn : R) (p : {poly R}) :=
  let h := w * hp in
  p.[tn + h] -
    (a1 w * p.[tn] + a0 w * p.[tn - hp] + beta w * h * p^`().[tn + h]).

(* On a cubic the residual is its coefficient of t^3 times one constant. *)
Theorem bdf2_residual_on_cubics (w hp tn : R) (p : {poly R}) :
  den w != 0 -> (size p <= 4)%N ->
  let h := w * hp in
  residual w hp tn p = - (p`_3 * (beta w * h ^+ 2 * (h + hp))).
Proof.
rewrite /residual /a1 /a0 /beta /den => nz sp /=.
have sd : (size p^`() <= 3)%N.
  by rewrite (leq_trans (size_poly _ _)) // -subn1 leq_subLR.
rewrite !(horner_coef_wide _ sp) (horner_coef_wide _ sd).
by rewrite !big_ord_recr !big_ord0 /= !coef_deriv; field.
Qed.

(* T1: exact on every polynomial of degree at most 2. *)
Theorem bdf2_exact_on_quadratics (w hp tn : R) (p : {poly R}) :
  den w != 0 -> (size p <= 3)%N ->
  let h := w * hp in
  p.[tn + h] =
    a1 w * p.[tn] + a0 w * p.[tn - hp] + beta w * h * p^`().[tn + h].
Proof.
move=> nz sp /=; apply/eqP; rewrite -subr_eq0.
have := bdf2_residual_on_cubics hp tn nz (leq_trans sp (leqnSn 3)).
by rewrite /residual /= => ->; rewrite (leq_sizeP _ _ sp) // mul0r oppr0.
Qed.

(* T1: the residual on t^3 is -beta h^2 (h + hp). *)
Theorem bdf2_cubic_residual (w hp tn : R) : den w != 0 ->
  let h := w * hp in residual w hp tn 'X^3 = - (beta w * h ^+ 2 * (h + hp)).
Proof.
move=> nz /=; rewrite (bdf2_residual_on_cubics hp tn nz) ?size_polyXn //.
by rewrite coefXn mul1r.
Qed.

(* The same value on the cubic (t - (tn + h))^3, with the origin at the
   new time point. *)
Lemma bdf2_cubic_residual_shifted (w hp : R) : den w != 0 ->
  let h := w * hp in
  a1 w * h ^+ 3 + a0 w * (h + hp) ^+ 3 = - (beta w * h ^+ 2 * (h + hp)).
Proof. by rewrite /a1 /a0 /beta /den => nz /=; field. Qed.

End Bdf2.
