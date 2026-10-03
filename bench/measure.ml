open Vstiff

type figures = { steps : int; rejected : int; rhs_calls : int; error : float; scd : float }
type outcome = Done of figures | Failed of Fail.t
type t = { problem : string; tol : float; outcome : outcome }

let largest f ~reference y = Array.fold_left Float.max 0. (Array.map2 f reference y)
let error ~reference y = largest (fun r v -> Float.abs (v -. r) /. (1. +. Float.abs r)) ~reference y

(* A zero reference component cannot be divided by; its absolute difference stands in. *)
let scd ~reference y =
  -.Float.log10 (largest (fun r v -> if r = 0. then Float.abs v else Float.abs (v -. r) /. Float.abs r) ~reference y)

(* Counting through Instrument keeps the mutable cell out of the solver; one fresh cell per run. *)
let run (case : Suite.case) ~tol =
  let rhs, calls = Instrument.count case.problem.rhs in
  let outcome =
    match Adaptive.integrate (module Bdf2) (module Halving) ~tol { case.problem with rhs } with
    | Ok s ->
        let reference = case.reference.y in
        Done
          {
            steps = s.stats.accepted_steps;
            rejected = s.stats.rejected_steps;
            rhs_calls = calls ();
            error = error ~reference s.y;
            scd = scd ~reference s.y;
          }
    | Error e -> Failed e
  in
  { problem = case.name; tol; outcome }
