(** Halve the step on rejection; double it, up to [dt_max], after three
    accepts in a row. *)

type stats = { accepted_steps : int; rejected_steps : int }

type t = {
  tol : float;
  dt_max : float;
  max_rejects : int;
  dt : float;  (** The step to try next. *)
  streak : int;  (** Accepts since [dt] last changed. *)
  failures : int;  (** Rejections in a row. *)
  stats : stats;
}

let init ~tol ~dt0 ~dt_max ~max_rejects =
  { tol; dt_max; max_rejects; dt = dt0; streak = 0; failures = 0; stats = { accepted_steps = 0; rejected_steps = 0 } }

let proposal c = c.dt

(* The largest component of [err], each measured against 1 + |y_i|. *)
let acceptable c ~y ~err =
  Array.fold_left Float.max 0. (Array.map2 (fun ei yi -> Float.abs ei /. (1. +. Float.abs yi)) err y) <= c.tol

let accepted c =
  let streak = c.streak + 1 in
  {
    c with
    dt = (if streak = 3 then Float.min (2. *. c.dt) c.dt_max else c.dt);
    streak = streak mod 3;
    failures = 0;
    stats = { c.stats with accepted_steps = c.stats.accepted_steps + 1 };
  }

(* Below this a step can move t by only a few units in the last place.
   Needing one means the solution is singular there, or the right-hand side
   fails just ahead, and halving further would never get past it. *)
let dt_min t = 16. *. Float.epsilon *. Float.abs t

(* Too large or unsolvable, a rejection halves the step that failed. *)
let rejected c (_ : Ode.rejection) ~at ~h =
  let failures = c.failures + 1 and dt = h /. 2. in
  if failures > c.max_rejects || dt < dt_min at then Error (Fail.StepRejected failures)
  else Ok { c with dt; streak = 0; failures; stats = { c.stats with rejected_steps = c.stats.rejected_steps + 1 } }

let stats c = c.stats
