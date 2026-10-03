(** The problems of the bench, each with its reference and the tolerances of the roadmap's Appendix B. *)

type case = {
  name : string;  (** as bench/compare/problems.py names it *)
  problem : Vstiff.Ode.problem;
  reference : Reference.t;
  tols : float list;  (** each is run as [rtol = atol = tol] *)
  box : (float * float) array;  (** per component, the range of the states that the transcription check draws *)
}

(** Every problem of the bench, in the order of Appendix B. *)
val cases : case list
