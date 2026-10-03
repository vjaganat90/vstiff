(** The step-size policy: a step is kept when its error estimate is at most
    [tol] in the norm [max_i |err_i| / (1 + |y_i|)]. A rejection, for either
    reason, halves the step that failed; three accepts in a row double the
    step, up to [dt_max]. The policy gives up with [StepRejected n] after
    [n > max_rejects] rejections in a row, or once a halved step falls below
    [16 eps |t|]. *)

type stats = { accepted_steps : int; rejected_steps : int }

include Ode.Controller with type stats := stats
