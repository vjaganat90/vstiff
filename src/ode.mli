(** The problem, and the contracts that methods and step-size controllers meet,
    as module types (OCaml interfaces; docs/ocaml.md). Declarations only, so
    there is no ode.ml (src/dune says modules_without_implementation). Who calls
    what: docs/architecture.md; terms: docs/glossary.md. *)

(** The right-hand side [f] of [y' = f(t, y)]. [rhs t y] returns a fresh vector as
    long as [y] on every call and does not modify [y]: the library keeps earlier
    results while it calls [rhs] again. *)
type rhs = float -> Vec.t -> Vec.t

(** Integrate [y' = rhs t y] from [y(t0) = y0] up to [t_end]. *)
type problem = { rhs : rhs; t0 : float; t_end : float; y0 : Vec.t }

(** A point [(t, y)] on a solution. *)
type point = { t : float; y : Vec.t }

(** Why an adaptive step was rejected. The controller is told so that it can
    treat the causes differently; [Halving] does not. *)
type rejection =
  | Too_large  (** The local error estimate exceeded the tolerance. *)
  | Solver of Fail.t  (** The method could not take the step. *)
  | Too_small  (** The step is below the resolution of [t], and might not move [t] at all: the method was not called. *)

(** A method that advances a solution one step at a time: [Stepper.fixed] runs
    any [Method], [Adaptive.integrate] an [Embedded] one. [step] is pure, so a
    driver can drop a step and retry from the same point and history. *)
module type Method = sig
  (** What the method remembers between steps: nothing for backward Euler, the
      previous step and state for BDF2. Drivers never look inside. *)
  type history

  (** The history before the first step. *)
  val start : history

  (** [step rhs h history at] is the state at [at.t + h] and the history for the
      next step, or [Error] if the step cannot be taken (for BDF, Newton failed).
      [history] is [start] or the one returned by the step that reached [at]. *)
  val step : rhs -> float -> history -> point -> (Vec.t * history, Fail.t) result
end

(** A method that also estimates the local error of each step: the error the
    step makes when it starts from the exact solution. [include Method] copies
    the items of [Method] here (docs/ocaml.md). *)
module type Embedded = sig
  include Method

  (** As [step], with the error estimate in the middle: [(y, err, history)]. *)
  val step_with_error : rhs -> float -> history -> point -> (Vec.t * Vec.t * history, Fail.t) result
end

(** A step-size policy for adaptive integration. For each attempt the driver
    takes the [proposal] (cutting the last step to land on [t_end]) and tries
    it. If the method returns a state, [acceptable] decides, and the driver
    calls [accepted] or [rejected]; if the method fails it calls [rejected]. *)
module type Controller = sig
  (** The policy's state: a value that each transition replaces. *)
  type t

  (** What the policy reports once the integration is done. The result type of
      [Adaptive.integrate] depends on it. *)
  type stats

  (** [init ~tol ~dt0 ~dt_max ~max_rejects] starts a run: [tol] bounds the error
      estimate, [dt0] is the first step to try, [dt_max] the longest step to
      propose, [max_rejects] how many rejections in a row to put up with. *)
  val init : tol:float -> dt0:float -> dt_max:float -> max_rejects:int -> t

  (** The step size to try next. *)
  val proposal : t -> float

  (** Whether a step that reached [y] with local error estimate [err] is good
      enough to keep. [y] and [err] have the same length. *)
  val acceptable : t -> y:Vec.t -> err:Vec.t -> bool

  (** The state after a step was accepted. *)
  val accepted : t -> t

  (** The state after the step of length [h] from time [at] was rejected, or
      [Error] when the policy gives up. [h] is the step attempted: shorter than
      the [proposal] when the driver cut it. *)
  val rejected : t -> rejection -> at:float -> h:float -> (t, Fail.t) result

  (** The [stats] so far. *)
  val stats : t -> stats
end
