(** Fixed-step integration. *)

(** [fixed (module M) ~dt problem] is the state at [t_end] after [round ((t_end - t0) / dt)] equal steps of [M]
    that cover the span (at least one unless it is empty), or the first failure. [dt] is only a target, and must not
    make that count exceed [max_int]: the conversion is unspecified there.
    @raise Invalid_argument unless [dt > 0] and [t_end >= t0]. *)
val fixed : (module M : Ode.Method) -> dt:float -> Ode.problem -> (Vec.t, Fail.t) result
