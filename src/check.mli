(** Argument checks for the drivers, the only library module that raises on purpose; a [nan] argument raises too. *)

(** @raise Invalid_argument unless [dt > 0] and [t_end >= t0], or when the span
    is not empty and [dt] is below {!Clock.resolution} of the larger of
    [|t0|] and [|t_end|]. *)
val fixed : dt:float -> Ode.problem -> unit

(** @raise Invalid_argument unless [t_end >= t0], and [dt0 > 0] when the span is not empty. *)
val adaptive : dt0:float -> Ode.problem -> unit
