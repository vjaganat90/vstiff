open Vstiff

(** The corpus problems. Each exposes only its right-hand side, initial data
    and, where one exists, the exact solution. *)

(** y' = -Λ y with Λ = diag(1, 100, 1e4): a swapped Jacobian row puts a fast
    rate on a slow component and wrecks the answer. *)
module Canary = struct
  let lambda = [| 1.; 100.; 1e4 |]
  let rhs _t y = Array.mapi (fun i yi -> -.lambda.(i) *. yi) y
  let y0 = [| 1.; 1.; 1. |]
  let exact t = Array.map (fun l -> exp (-.l *. t)) lambda
  let problem = { Ode.rhs; t0 = 0.; t_end = 1.; y0 }
end

(** y' = y (1 - y), y(0) = 0.1: smooth, with a closed form. *)
module Logistic = struct
  let rhs _t y = [| y.(0) *. (1. -. y.(0)) |]
  let y0 = [| 0.1 |]
  let exact t = [| 1. /. (1. +. (((1. /. y0.(0)) -. 1.) *. exp (-.t))) |]
  let problem = { Ode.rhs; t0 = 0.; t_end = 5.; y0 }
end

(** van der Pol, μ = 1000: slow drifts broken by jumps on a 1/μ time scale. *)
module Van_der_pol = struct
  let mu = 1000.
  let rhs _t y = [| y.(1); (mu *. (1. -. (y.(0) *. y.(0))) *. y.(1)) -. y.(0) |]
  let y0 = [| 2.; 0. |]
  let problem = { Ode.rhs; t0 = 0.; t_end = 2000.; y0 }
end

(** Robertson's chemical kinetics: rates spanning nine orders of magnitude,
    and y1 + y2 + y3 is invariant. *)
module Robertson = struct
  let rhs _t y =
    let a = 0.04 *. y.(0) and b = 1e4 *. y.(1) *. y.(2) and c = 3e7 *. y.(1) *. y.(1) in
    [| b -. a; a -. b -. c; c |]

  let y0 = [| 1.; 0.; 0. |]
  let problem = { Ode.rhs; t0 = 0.; t_end = 1e4; y0 }
end
