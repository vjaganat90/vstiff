(** Fixed-step integration. *)

(** [fixed (module M) ~dt problem] is the state at [t_end] after [n = round ((t_end - t0) / dt)] steps of [M]
    (at least one unless the span is empty), or the first failure. Step [k] ends at [t0 + k h], [h = span / n], and
    the last one at [t_end] itself; each step is the difference of its end times, so the state advances by exactly
    what the clock does. [Error Nan] if [y0] or [rhs t0 y0] is not finite.
    @raise Invalid_argument unless [dt > 0] and [t_end >= t0], if [dt] is below {!Clock.resolution} of [t], or if
    [rhs t0 y0] is not as long as [y0]. *)
val fixed : (module M : Ode.Method) -> dt:float -> Ode.problem -> (Vec.t, Fail.t) result
