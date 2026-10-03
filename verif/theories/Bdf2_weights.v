(* The weights of variable-step BDF2, as Bdf2.coeffs computes them for the
   step ratio w = h / h_prev, over any field. See
   docs/plans/formal-verification.md for the plan this starts. *)
From mathcomp Require Import boot order algebra.

Set Implicit Arguments.
Unset Strict Implicit.
Unset Printing Implicit Defensive.

Import GRing.Theory.
Local Open Scope ring_scope.

Section Bdf2.

Variable R : fieldType.

Definition den (w : R) := 1 + w *+ 2.
Definition a1 (w : R) := (1 + w) ^+ 2 / den w.
Definition a0 (w : R) := - (w ^+ 2) / den w.

(* The history weights sum to one, so the formula is exact on constants. *)
Lemma weights_sum_to_one (w : R) : den w != 0 -> a1 w + a0 w = 1.
Proof.
move=> nz; rewrite /a1 /a0 -mulrDl.
have -> : (1 + w) ^+ 2 + - w ^+ 2 = den w.
  by rewrite /den sqrrD expr1n mul1r addrK addrC.
by rewrite divff.
Qed.

End Bdf2.
