(* Every proposal is at most twice the last accepted step, so BDF2's step ratio
   omega = h / h_prev stays at most 2, inside its limit 1 + sqrt 2, up to the
   driver snapping steps of a few ulps to the floats. The estimate
   only decides accept or reject; production controllers also use its size.
   See docs/numerics/05-step-control.md. *)

type stats = { accepted_steps : int; rejected_steps : int }

type t = {
  tol : float;
  dt_max : float;
  max_rejects : int;
  dt : float;
  streak : int;  (* Accepts towards the next doubling: 0, 1 or 2. A rejection resets it. *)
  failures : int;  (* Rejections in a row. *)
  stats : stats;
}

let init ~tol ~dt0 ~dt_max ~max_rejects =
  { tol; dt_max; max_rejects; dt = dt0; streak = 0; failures = 0; stats = { accepted_steps = 0; rejected_steps = 0 } }

let proposal c = c.dt

(* 1 + |y_i| makes tol relative for large components and absolute for small
   ones, and is never zero. A nan makes the norm nan, and nan <= tol is false,
   so the step is rejected. Unequal lengths raise Invalid_argument. *)
let acceptable c ~y ~err =
  Array.fold_left Float.max 0. (Array.map2 (fun ei yi -> Float.abs ei /. (1. +. Float.abs yi)) err y) <= c.tol

(* [{ c with ... }] copies c with those fields replaced (docs/ocaml.md). The
   streak restarts at every third accept, even when dt_max stops dt growing. *)
let accepted c =
  let streak = c.streak + 1 in
  {
    c with
    dt = (if streak = 3 then Float.min (2. *. c.dt) c.dt_max else c.dt);
    streak = streak mod 3;
    failures = 0;
    stats = { c.stats with accepted_steps = c.stats.accepted_steps + 1 };
  }

(* Every reason halves the step that failed, not the proposal: the driver may
   have cut the last step, and halving the proposal could retry the same one.
   A shorter step also eases a failed solve, as I - gamma J nears I. *)
let rejected c (_ : Ode.rejection) ~at ~h =
  let failures = c.failures + 1 and dt = h /. 2. in
  (* An accept resets failures, so a run closing in on a blow-up or a failing rhs
     may never reach max_rejects: the floor ends it. *)
  (* At t = 0 the floor is 0, so a halved step that underflows to 0 ends the run too. *)
  if failures > c.max_rejects || dt < Clock.resolution at || not (dt > 0.) then Error (Fail.StepRejected failures)
  else Ok { c with dt; streak = 0; failures; stats = { c.stats with rejected_steps = c.stats.rejected_steps + 1 } }

let stats c = c.stats
