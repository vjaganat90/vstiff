(** Variable-step BDF2. *)

(** The step solves [y_{n+1} = a1 y_n + a0 y_{n-1} + beta h f(t_{n+1}, y_{n+1})]. *)
type coeffs = { a1 : float; a0 : float; beta : float }

(** [coeffs omega] for the step ratio [omega = h / h_prev]: the weights that
    make the formula exact on quadratics through the last three points. *)
val coeffs : float -> coeffs

(** The previous step: its length and the state it started from. *)
type history = { dt_prev : float; y_prev : Vec.t }

(** [step rhs h history at] is the state at [at.t + h]. *)
val step : Ode.rhs -> float -> history -> Ode.point -> (Vec.t, Fail.t) result

(** The final state after equal steps of about [dt], the first one by
    backward Euler since BDF2 needs a previous step.
    @raise Invalid_argument as {!Stepper.fixed}. *)
val integrate : dt:float -> Ode.problem -> (Vec.t, Fail.t) result
