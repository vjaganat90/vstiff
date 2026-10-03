(** Adaptive integration with step rejection. *)

(** The end of a run: the state [y] at time [t], exactly [t_end], and what the
    controller reports in [stats]. The type depends on the controller, which is
    why [integrate] takes its modules as modular explicits (docs/ocaml.md). *)
type 'stats solution = { t : float; y : Vec.t; stats : 'stats }

(** [integrate (module M) (module C) ~tol problem] advances [problem.y0] from [t0] to [t_end], the last step cut to
    land on [t_end]. Every other step is snapped to the floats, [h = (t + dt) - t], so the state advances by exactly
    what the clock does; a step that cannot move [t] is rejected as [Too_small] without calling [M].
    [M] steps with an error estimate; [C] decides whether to keep a step and reports its [stats].
    [tol] limits each step's estimate, not the final error. Defaults: [dt0 = 1e-6 (t_end - t0)], capped at
    [dt_max = (t_end - t0) / 10]; [max_rejects = 50]. [Error Nan] if [y0] or [rhs t0 y0] is not finite, and [C]'s
    [Error] when it gives up.
    @raise Invalid_argument if [t_end < t0], the span is not empty and [dt0] or [dt_max] is not positive, [tol] is
    not positive, or [rhs t0 y0] is not as long as [y0]. *)
val integrate :
  (module M : Ode.Embedded) ->
  (module C : Ode.Controller) ->
  ?dt0:float ->
  ?dt_max:float ->
  ?max_rejects:int ->
  tol:float ->
  Ode.problem ->
  (C.stats solution, Fail.t) result
