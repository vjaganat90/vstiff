(* a.(i) indexes an array; float operators end in a dot; [( *. ) s] is a partial application: docs/ocaml.md. *)

type t = float array

let add a b = Array.init (Array.length a) (fun i -> a.(i) +. b.(i))
let sub a b = Array.init (Array.length a) (fun i -> a.(i) -. b.(i))
let scale s = Array.map (( *. ) s)

let axpy a x y = Array.init (Array.length x) (fun i -> (a *. x.(i)) +. y.(i))

let dot a b = Array.fold_left ( +. ) 0. (Array.init (Array.length a) (fun i -> a.(i) *. b.(i)))

(* Float.max returns nan if either argument is nan; the generic max can drop it. *)
let norm_inf = Array.fold_left (fun m v -> Float.max m (Float.abs v)) 0.
let finite = Array.for_all Float.is_finite
