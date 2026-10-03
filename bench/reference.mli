(** Reference states computed outside vstiff, with scipy, by bench/compare/references.py. Each keeps only the digits on
    which independent solvers agree; reference.ml gives, next to every value, the solvers, versions, tolerances and the
    agreement. *)

(** The state [y] at [t_end] of one problem. *)
type t = { t_end : float; y : float array }

(** Robertson from [(1, 0, 0)] at [t = 1e4]. *)
val robertson : t

(** HIRES at [t = 321.8122]. *)
val hires : t

(** van der Pol with [mu = 1000] from [(2, 0)] at [t = 2000]. *)
val van_der_pol : t

(** The Brusselator on 40 cells (n = 80) at [t = 10]. *)
val brusselator_80 : t
