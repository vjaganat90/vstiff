(** Damped Newton for [f x = 0] with a caller-supplied Jacobian. *)

type config = {
  tol : float;  (** Stop once [|dx|_inf <= tol (1 + |x|_inf)]. *)
  max_iter : int;
  min_damping : float;  (** Smallest line-search factor tried before giving up. *)
}

let default = { tol = 1e-10; max_iter = 50; min_damping = 1. /. 1024. }

(* Armijo constant: a damped step must shrink the residual by this fraction of
   the decrease the linear model promises. *)
let armijo = 1e-4

let solve ?(config = default) (f : Vec.t -> Vec.t) (jac : Vec.t -> Linalg.matrix) (x0 : Vec.t) :
    (Vec.t, Fail.t) result =
  let { tol; max_iter; min_damping } = config in
  let rec damp x r dx lambda =
    if lambda < min_damping then Error Fail.Diverged
    else
      let x' = Vec.axpy lambda dx x in
      let fx' = f x' in
      if Vec.finite fx' && Vec.norm_inf fx' <= (1. -. (armijo *. lambda)) *. r then Ok (x', fx')
      else damp x r dx (lambda /. 2.)
  in
  let rec iterate k x fx =
    if not (Vec.finite x && Vec.finite fx) then Error Fail.Nan
    else if Vec.norm_inf fx = 0. then Ok x
    else if k >= max_iter then Error Fail.Diverged
    else
      match Linalg.solve (jac x) (Vec.scale (-1.) fx) with
      | None -> Error Fail.Diverged
      | Some dx when not (Vec.finite dx) -> Error Fail.Diverged
      | Some dx when Vec.norm_inf dx <= tol *. (1. +. Vec.norm_inf x) -> Ok (Vec.add x dx)
      | Some dx -> Result.bind (damp x (Vec.norm_inf fx) dx 1.) (fun (x', fx') -> iterate (k + 1) x' fx')
  in
  iterate 0 x0 (f x0)
