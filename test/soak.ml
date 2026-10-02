(* Every corpus problem, ten rounds each: each round must pass its criterion
   and every round must reproduce the first bit for bit. *)
open Vstiff

let rounds = 10

module type CASE = sig
  type t

  val name : string
  val run : unit -> (t, Fail.t) result
  val pass : t -> bool
end

let soak (module C : CASE) =
  let results = List.init rounds (fun _ -> C.run ()) in
  let passed = List.length (List.filter (function Ok r -> C.pass r | Error _ -> false) results) in
  let identical = List.for_all (( = ) (List.hd results)) results in
  Printf.printf "soak %s x%d: passed %d/%d, identical: %b\n" C.name rounds passed rounds identical

let max_error exact y = Vec.norm_inf (Vec.sub y exact)

module Canary = struct
  open Problems.Canary

  type t = float

  let name = "canary"
  let run () = Result.map (max_error (exact 1.)) (Bdf1.integrate ~rhs ~t0:0. ~t_end:1. ~dt:2e-6 y0)
  let pass e = e < 1e-6
end

module Logistic_order = struct
  open Problems.Logistic

  type t = float * float

  let name = "logistic order"
  let error dt = Result.map (max_error (exact 5.)) (Bdf2.integrate ~rhs ~t0:0. ~t_end:5. ~dt y0)
  let run () = Result.bind (error 0.01) (fun coarse -> Result.map (fun fine -> (coarse, fine)) (error 0.005))
  let pass (coarse, fine) = 3.5 <= coarse /. fine && coarse /. fine <= 4.5
end

module Van_der_pol = struct
  open Problems.Van_der_pol

  type t = Adaptive.solution

  let name = "van der Pol"
  let run () = Adaptive.integrate ~tol:1e-4 ~rhs ~t0:0. ~t_end:2000. y0
  let pass (s : t) = s.t = 2000. && s.rejected >= 1 && Vec.finite s.y
end

module Robertson = struct
  open Problems.Robertson

  type t = Adaptive.solution

  let name = "robertson"
  let run () = Adaptive.integrate ~tol:1e-6 ~rhs ~t0:0. ~t_end:1e4 y0

  let pass (s : t) =
    Float.abs (s.y.(0) -. Refs.robertson_y1_at_1e4) < 1e-3
    && Float.abs (Array.fold_left ( +. ) 0. s.y -. 1.) < 1e-8
end

let () =
  List.iter soak
    [ (module Canary : CASE); (module Logistic_order); (module Van_der_pol); (module Robertson) ]
