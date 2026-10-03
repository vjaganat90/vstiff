(** Backward Euler: [y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})]. *)

(** [step rhs h at] is the state at [at.t + h]. *)
val step : Ode.rhs -> float -> Ode.point -> (Vec.t, Fail.t) result

(** The final state after equal steps of about [dt].
    @raise Invalid_argument as {!Stepper.fixed}. *)
val integrate : dt:float -> Ode.problem -> (Vec.t, Fail.t) result
