(** Fixed-step integration. *)

(** [fixed (module M) ~dt problem] is the state at [t_end] after
    [round ((t_end - t0) / dt)] equal steps of method [M] (at least one when
    [t_end > t0]), sized so the last one lands exactly on [t_end]; or the
    first failure.
    @raise Invalid_argument unless [dt > 0] and [t_end >= t0]. *)
val fixed : (module M : Ode.Method) -> dt:float -> Ode.problem -> (Vec.t, Fail.t) result
