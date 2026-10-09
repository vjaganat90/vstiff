(** Argument checks for the drivers, the only library module that raises on purpose; a [nan] argument raises too. *)

(** @raise Invalid_argument unless [dt > 0] and [t_end >= t0], or when the span
    is not empty and [dt] is below {!Clock.resolution} of the larger of
    [|t0|] and [|t_end|]. *)
val fixed : dt:float -> Ode.problem -> unit

(** @raise Invalid_argument unless [t_end >= t0], [dt0 > 0] when the span is not empty, and [tol > 0]. *)
val adaptive : dt0:float -> tol:float -> Ode.problem -> unit

(** [output ~caller problem f0] checks [f0 = rhs t0 y0], which each driver evaluates once.
    @raise Invalid_argument, naming [caller], unless [f0] is as long as [y0]. *)
val output : caller:string -> Ode.problem -> float array -> unit
