open Vstiff

(** The corpus problems. Each module has its right-hand side [rhs], initial state [y0], the [problem] record the
    integrators take and, where a closed form exists, [exact]. What each proves: docs/numerics/06-the-corpus.md. *)

(** Three independent decays: y' = -Lambda y, Lambda = diag(1, 100, 1e4), y(0) = (1, 1, 1), exact
    y_i = exp(-lambda_i t). Stiff because of the fast rate. J = -Lambda is diagonal and known exactly, so a
    mix-up between components shows and a transposed J does not. *)
module Canary = struct
  let lambda = [| 1.; 100.; 1e4 |]
  let rhs _t y = Array.mapi (fun i yi -> -.lambda.(i) *. yi) y
  let y0 = [| 1.; 1.; 1. |]
  let exact t = Array.map (fun l -> exp (-.l *. t)) lambda
  let problem = { Ode.rhs; t0 = 0.; t_end = 1.; y0 }
end

(** Growth that levels off: y' = y (1 - y), y(0) = 0.1, exact y = 1 / (1 + 9 exp(-t)). Smooth and not stiff, so it
    measures accuracy: the order test and the step-control counts run on it. *)
module Logistic = struct
  let rhs _t y = [| y.(0) *. (1. -. y.(0)) |]
  let y0 = [| 0.1 |]
  let exact t = [| 1. /. (1. +. (((1. /. y0.(0)) -. 1.) *. exp (-.t))) |]
  let problem = { Ode.rhs; t0 = 0.; t_end = 5.; y0 }
end

(** van der Pol with mu = 1000: y1' = y2, y2' = mu (1 - y1^2) y2 - y1, from (2, 0). A relaxation oscillation: slow
    drifts broken by jumps near t = 807 and t = 1614 (leading-order theory). No closed form and no reference
    value, so it tests step control. *)
module VanDerPol = struct
  let mu = 1000.
  let rhs _t y = [| y.(1); (mu *. (1. -. (y.(0) *. y.(0))) *. y.(1)) -. y.(0) |]
  let y0 = [| 2.; 0. |]
  let problem = { Ode.rhs; t0 = 0.; t_end = 2000.; y0 }
end

(** Robertson's chemical kinetics from y(0) = (1, 0, 0): y1' = -0.04 y1 + 1e4 y2 y3,
    y2' = 0.04 y1 - 1e4 y2 y3 - 3e7 y2^2, y3' = 3e7 y2^2. A classic stiff problem: rate constants from 0.04 to
    3e7, and y2 stays tiny (peak about 3.6e-5). The right-hand sides sum to zero and BDF weights sum to one, so
    y1 + y2 + y3 stays 1 to round-off. *)
module Robertson = struct
  let rhs _t y =
    let a = 0.04 *. y.(0) and b = 1e4 *. y.(1) *. y.(2) and c = 3e7 *. y.(1) *. y.(1) in
    [| b -. a; a -. b -. c; c |]

  let y0 = [| 1.; 0.; 0. |]
  let problem = { Ode.rhs; t0 = 0.; t_end = 1e4; y0 }
end
