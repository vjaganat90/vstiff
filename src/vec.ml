(** Dense vectors as immutable-by-convention [float array]s: every operation
    allocates its result and never writes to an argument. *)

type t = float array

let add a b = Array.init (Array.length a) (fun i -> a.(i) +. b.(i))
let sub a b = Array.init (Array.length a) (fun i -> a.(i) -. b.(i))
let scale s = Array.map (( *. ) s)

(** [axpy a x y] is [a x + y]. *)
let axpy a x y = Array.init (Array.length x) (fun i -> (a *. x.(i)) +. y.(i))

let dot a b = Array.fold_left ( +. ) 0. (Array.init (Array.length a) (fun i -> a.(i) *. b.(i)))
let norm_inf = Array.fold_left (fun m v -> Float.max m (Float.abs v)) 0.
let finite = Array.for_all Float.is_finite
