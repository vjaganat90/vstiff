(** Forward-difference Jacobian. *)

(** Perturbation for component [yj]: relative for large values, absolute near zero. *)
let step yj = 1e-8 *. (1. +. Float.abs yj)

(** [forward f y] is the matrix [J] with [J.(i).(j)] ≈ ∂f_i/∂y_j at [y]. *)
let forward (f : Vec.t -> Vec.t) (y : Vec.t) : float array array =
  let fy = f y in
  let column j =
    let yp = Array.mapi (fun k yk -> if k = j then yk +. step yk else yk) y in
    (* Divide by the perturbation actually representable in [yp], not the nominal one. *)
    Vec.scale (1. /. (yp.(j) -. y.(j))) (Vec.sub (f yp) fy)
  in
  let cols = Array.init (Array.length y) column in
  Array.init (Array.length fy) (fun i -> Array.map (fun col -> col.(i)) cols)
