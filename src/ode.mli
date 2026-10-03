(** The vocabulary every integrator shares: what a problem is, and what a
    method and a step-size controller must provide. *)

(** The right-hand side [f] of [y' = f(t, y)]. *)
type rhs = float -> Vec.t -> Vec.t

(** Integrate [y' = rhs t y] from [y(t0) = y0] up to [t_end]. *)
type problem = { rhs : rhs; t0 : float; t_end : float; y0 : Vec.t }

(** A point [(t, y)] on a solution. *)
type point = { t : float; y : Vec.t }

(** Why an adaptive step was rejected. *)
type rejection =
  | Too_large  (** The local error estimate exceeded the tolerance. *)
  | Solver of Fail.t  (** The method could not take the step. *)

(** A method that advances a solution one step at a time. *)
module type Method = sig
  (** What the method carries from one step to the next. *)
  type history

  (** The history before the first step. *)
  val start : history

  (** [step rhs h history at] is the state at [at.t + h] and the history for
      the step after it. *)
  val step : rhs -> float -> history -> point -> (Vec.t * history, Fail.t) result
end

(** A method that also estimates the local error of each step. *)
module type Embedded = sig
  include Method

  (** As [step], with an estimate of the step's local error in the middle. *)
  val step_with_error : rhs -> float -> history -> point -> (Vec.t * Vec.t * history, Fail.t) result
end

(** A step-size policy for adaptive integration. *)
module type Controller = sig
  (** The policy's state. *)
  type t

  (** What the policy reports once the integration is done. *)
  type stats

  val init : tol:float -> dt0:float -> dt_max:float -> max_rejects:int -> t

  (** The step size to try next. *)
  val proposal : t -> float

  (** Whether a step that reached [y] with local error estimate [err] is good
      enough to keep. *)
  val acceptable : t -> y:Vec.t -> err:Vec.t -> bool

  (** The state after a step was accepted. *)
  val accepted : t -> t

  (** The state after a step of length [h] from time [at] was rejected, or
      [Error] when the policy gives up. *)
  val rejected : t -> rejection -> at:float -> h:float -> (t, Fail.t) result

  val stats : t -> stats
end
