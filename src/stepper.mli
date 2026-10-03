(** Fixed-step driving for one-step methods. *)

(** [fixed ~dt problem advance start] applies [advance h] to [start] once per
    step, [round ((t_end - t0) / dt)] times (at least once when
    [t_end > t0]), with [h] chosen so the last step lands on [t_end]. Stops at
    the first failure.
    @raise Invalid_argument unless [dt > 0] and [t_end >= t0]. *)
val fixed : dt:float -> Ode.problem -> (float -> 's -> ('s, Fail.t) result) -> 's -> ('s, Fail.t) result
