(** A step is kept when [max_i |err_i| / (1 + |y_i|) <= tol]; a rejection halves
    the step that failed, and three accepts in a row double it, up to [dt_max].
    Gives up with [StepRejected n] after [n > max_rejects] rejections in a row,
    or once a halved step is below [16 eps |t|] ([eps] is [Float.epsilon]) or is 0. *)

(** Counts over the whole run; every rejection counts, whatever its reason. *)
type stats = { accepted_steps : int; rejected_steps : int }

(* [with type stats := stats] replaces the abstract [stats] of the signature by
   the record above: docs/ocaml.md. *)
include Ode.Controller with type stats := stats
