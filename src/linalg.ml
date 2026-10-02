(** Dense linear solve by Gaussian elimination with partial pivoting, written
    as a recursion on the trailing submatrix so no array is ever mutated. *)

(* Index of the row with the largest |a.(i).(0)|. *)
let pivot a =
  let best i j = if Float.abs a.(j).(0) > Float.abs a.(i).(0) then j else i in
  Seq.fold_left best 0 (Seq.init (Array.length a) Fun.id)

let swap p i = if i = 0 then p else if i = p then 0 else i

(** [solve a b] is [Some x] with [a x = b], or [None] when elimination meets
    an exactly zero pivot. *)
let rec solve (a : float array array) (b : Vec.t) : Vec.t option =
  match Array.length b with
  | 0 -> Some [||]
  | n ->
      let p = pivot a in
      let row = a.(p) and rhs = b.(p) in
      if row.(0) = 0. then None
      else
        let rest = Array.init (n - 1) (fun i -> swap p (i + 1)) in
        let l i = a.(i).(0) /. row.(0) in
        let a' =
          Array.map (fun i -> Array.init (n - 1) (fun j -> a.(i).(j + 1) -. (l i *. row.(j + 1)))) rest
        and b' = Array.map (fun i -> b.(i) -. (l i *. rhs)) rest in
        Option.map
          (fun x ->
            let tail = Array.sub row 1 (n - 1) in
            Array.append [| (rhs -. Vec.dot tail x) /. row.(0) |] x)
          (solve a' b')
