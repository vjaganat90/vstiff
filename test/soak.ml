(* Every corpus problem, ten rounds each: each round must pass its criterion
   and every round must reproduce the first bit for bit. *)
open Vstiff

let rounds = 10

module type Case = sig
  type t

  val name : string
  val run : unit -> (t, Fail.t) result
  val pass : t -> bool
end

(* Every round's result. The result type depends on the module passed in, which
   takes a modular explicit: (module C : Case) -> (C.t, Fail.t) result list. *)
let repeat (module C : Case) : (C.t, Fail.t) result list = List.init rounds (fun _ -> C.run ())

let soak (module C : Case) =
  let results = repeat (module C) in
  let passed = List.length (List.filter (function Ok r -> C.pass r | Error _ -> false) results) in
  let identical = List.for_all (( = ) (List.hd results)) results in
  Printf.printf "soak %s x%d: passed %d/%d, identical: %b\n" C.name rounds passed rounds identical

let max_error exact y = Vec.norm_inf (Vec.sub y exact)
let bdf2_halving = Adaptive.integrate (module Bdf2) (module Halving)

module Canary = struct
  type t = float

  let name = "canary"

  let run () =
    Result.map (max_error (Problems.Canary.exact 1.)) (Stepper.fixed (module Bdf1) ~dt:2e-6 Problems.Canary.problem)

  let pass e = e < 1e-6
end

module LogisticOrder = struct
  type t = float * float

  let name = "logistic order"

  let error dt =
    Result.map (max_error (Problems.Logistic.exact 5.)) (Stepper.fixed (module Bdf2) ~dt Problems.Logistic.problem)

  let run () = Result.bind (error 0.01) (fun coarse -> Result.map (fun fine -> (coarse, fine)) (error 0.005))
  let pass (coarse, fine) = 3.5 <= coarse /. fine && coarse /. fine <= 4.5
end

module VanDerPol = struct
  type t = Halving.stats Adaptive.solution

  let name = "van der Pol"
  let run () = bdf2_halving ~tol:1e-4 Problems.VanDerPol.problem
  let pass (s : t) = s.t = 2000. && s.stats.rejected_steps >= 1 && Vec.finite s.y
end

module Robertson = struct
  type t = Halving.stats Adaptive.solution

  let name = "robertson"
  let run () = bdf2_halving ~tol:1e-6 Problems.Robertson.problem

  let pass (s : t) =
    Float.abs (s.y.(0) -. Refs.robertson_y1_at_1e4) < 1e-3 && Float.abs (Array.fold_left ( +. ) 0. s.y -. 1.) < 1e-8
end

let () =
  soak (module Canary);
  soak (module LogisticOrder);
  soak (module VanDerPol);
  soak (module Robertson)
