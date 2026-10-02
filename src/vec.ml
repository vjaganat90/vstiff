(** Dense vectors as immutable-by-convention [float array]s: every operation
    allocates its result and never writes to an argument. *)

type t = float array

let map2 f a b = Array.init (Array.length a) (fun i -> f a.(i) b.(i))
let add = map2 ( +. )
let sub = map2 ( -. )
let scale s = Array.map (( *. ) s)

(** [axpy a x y] is [a x + y]. *)
let axpy a x y = map2 (fun xi yi -> (a *. xi) +. yi) x y

let dot a b = Array.fold_left ( +. ) 0. (map2 ( *. ) a b)
let norm_inf = Array.fold_left (fun m v -> Float.max m (Float.abs v)) 0.
let finite = Array.for_all Float.is_finite
