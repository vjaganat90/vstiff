(** Damped Newton for [f x = 0] with a caller-supplied Jacobian. *)

(* Converged once |dx|_inf <= tol (1 + |x|_inf). *)
let tol = 1e-10
let max_iter = 50

(* Smallest line-search factor tried before giving up. *)
let min_damping = 1. /. 1024.

(* Armijo constant: a damped step must shrink the residual by this fraction of
   the decrease the linear model promises. *)
let armijo = 1e-4

let solve (f : Vec.t -> Vec.t) (jac : Vec.t -> Linalg.matrix) (x0 : Vec.t) : (Vec.t, Fail.t) result =
  let rec iterate k x fx =
    let r = Vec.norm_inf fx in
    if not (Vec.finite x && Vec.finite fx) then Error Fail.Nan
    else if r = 0. then Ok x
    else if k >= max_iter then Error Fail.Diverged
    else
      match Linalg.solve (jac x) (Vec.scale (-1.) fx) with
      | None -> Error Fail.Diverged
      | Some dx when not (Vec.finite dx) -> Error Fail.Diverged
      | Some dx when Vec.norm_inf dx <= tol *. (1. +. Vec.norm_inf x) -> Ok (Vec.add x dx)
      | Some dx ->
          (* Backtrack: halve the step until the residual drops enough. *)
          let rec damp lambda =
            if lambda < min_damping then Error Fail.Diverged
            else
              let x' = Vec.axpy lambda dx x in
              let fx' = f x' in
              if Vec.finite fx' && Vec.norm_inf fx' <= (1. -. (armijo *. lambda)) *. r then iterate (k + 1) x' fx'
              else damp (lambda /. 2.)
          in
          damp 1.
  in
  iterate 0 x0 (f x0)
