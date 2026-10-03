(** Adaptive BDF2 with step rejection. *)

type solution = { t : float; y : Vec.t; accepted_steps : int; rejected_steps : int }

(** [integrate ~tol problem] advances [problem.y0] from [t0] to [t_end]. A step
    is accepted when its local error estimate, the gap between backward Euler
    and BDF2 from the same history, is at most [tol] in the norm
    [max_i |e_i| / (1 + |y_i|)]; otherwise the step is halved. After three
    accepts in a row the step doubles, up to [dt_max].

    [dt0] defaults to [1e-6 (t_end - t0)], [dt_max] to [(t_end - t0) / 10],
    [max_rejects] to [50]. A failed Newton solve counts as a rejection.
    Returns [Error (StepRejected n)] after [n > max_rejects] rejections in a
    row or once a halved step falls below [16 eps |t|], and [Error Nan] if
    [y0] or [rhs t0 y0] is not finite.
    @raise Invalid_argument if [t_end < t0] or a step bound is not positive. *)
val integrate :
  ?dt0:float -> ?dt_max:float -> ?max_rejects:int -> tol:float -> Ode.problem -> (solution, Fail.t) result
