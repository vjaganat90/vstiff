(** The vocabulary every integrator shares. *)

(** The right-hand side [f] of [y' = f(t, y)]. *)
type rhs = float -> Vec.t -> Vec.t

(** Integrate [y' = rhs t y] from [y(t0) = y0] up to [t_end]. *)
type problem = { rhs : rhs; t0 : float; t_end : float; y0 : Vec.t }

(** A point [(t, y)] on a solution. *)
type point = { t : float; y : Vec.t }
