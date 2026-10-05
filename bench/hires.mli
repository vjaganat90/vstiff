(** HIRES: eight species of a plant's response to light, with one nonlinear term, from [(1, 0, 0, 0, 0, 0, 0, 0.0057)]
    to [t = 321.8122]. Moderately stiff and the standard first test of a stiff solver. *)

(** The right-hand side, spelled as in bench/compare/problems.py. *)
val rhs : Vstiff.Ode.rhs

val y0 : float array
val problem : Vstiff.Ode.problem
