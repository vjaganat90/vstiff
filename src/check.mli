(** Argument checks for the drivers. This is the only module in the library
    that raises: a bad argument is a programming error, not a numerical
    failure. *)

(** @raise Invalid_argument unless [dt > 0] and [t_end >= t0]. *)
val fixed : dt:float -> Ode.problem -> unit

(** @raise Invalid_argument if [t_end < t0], or if the span is not empty and
    [dt0] is not positive. *)
val adaptive : dt0:float -> Ode.problem -> unit
