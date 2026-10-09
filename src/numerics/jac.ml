(* Truncation error, about (step / 2) |f''|, falls with step; round-off, about eps |f| / step, grows as step falls
   (eps = Float.epsilon). For |f| and |f''| of order 1 they meet near sqrt eps = 1.5e-8: 1e-8 is a round number close
   by; 1 + |y_j| makes it relative for large y_j, nonzero near 0. docs/numerics/03-jacobians-and-floating-point.md *)
let step yj = 1e-8 *. (1. +. Float.abs yj)

let forward (f : Vec.t -> Vec.t) (y : Vec.t) : Linalg.matrix =
  let fy = f y in
  let column j =
    let yp = Array.mapi (fun k yk -> if k = j then yk +. step yk else yk) y in
    (* Divide by the perturbation stored in [yp], not by step: y_j + step is rounded, so the two can differ. *)
    (1. /. (yp.(j) -. y.(j)), f yp)
  in
  let cols = Array.init (Array.length y) column in
  Array.mapi (fun i fyi -> Array.map (fun (inv, fp) -> inv *. (fp.(i) -. fyi)) cols) fy
