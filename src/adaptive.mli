(** Adaptive integration with step rejection. *)

type 'stats solution = { t : float; y : Vec.t; stats : 'stats }

(** [integrate (module M) (module C) ~tol problem] advances [problem.y0] from
    [t0] to [t_end]. Each attempt takes a step of [M] with its error estimate;
    [C] decides whether to keep it and what to try next, and reports its
    [stats] at the end. A failed step counts as a rejection.

    [dt0] defaults to [1e-6 (t_end - t0)], [dt_max] to [(t_end - t0) / 10],
    [max_rejects] to [50]. Returns [Error Nan] if [y0] or [rhs t0 y0] is not
    finite, and [C]'s [Error] when it gives up.
    @raise Invalid_argument if [t_end < t0] or a step bound is not positive. *)
val integrate :
  (module M : Ode.Embedded) ->
  (module C : Ode.Controller) ->
  ?dt0:float ->
  ?dt_max:float ->
  ?max_rejects:int ->
  tol:float ->
  Ode.problem ->
  (C.stats solution, Fail.t) result
