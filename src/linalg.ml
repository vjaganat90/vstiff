(** Dense linear solve by Gaussian elimination with partial pivoting, written
    as a recursion on the trailing submatrix so no array is ever mutated. *)

type matrix = float array array

(* Index of the row with the largest |a.(i).(0)|; the first of equals. *)
let pivot (a : matrix) =
  let rec best i j =
    if j = Array.length a then i else best (if Float.abs a.(j).(0) > Float.abs a.(i).(0) then j else i) (j + 1)
  in
  best 0 1

(* Which original row sits at position [i] once rows 0 and [p] trade places. *)
let swapped p i = if i = 0 then p else if i = p then 0 else i

(* [rows] are augmented: each holds its coefficients, then its right-hand side.
   Eliminating column 0 leaves the same problem with one unknown fewer. *)
let rec solve_augmented (rows : matrix) : Vec.t option =
  match Array.length rows with
  | 0 -> Some [||]
  | n ->
      let p = pivot rows in
      let top = rows.(p) in
      if top.(0) = 0. then None
      else
        let reduce row =
          let m = row.(0) /. top.(0) in
          Array.init n (fun j -> row.(j + 1) -. (m *. top.(j + 1)))
        in
        let rest = Array.init (n - 1) (fun i -> reduce rows.(swapped p (i + 1))) in
        Option.map
          (fun x -> Array.append [| (top.(n) -. Vec.dot (Array.sub top 1 (n - 1)) x) /. top.(0) |] x)
          (solve_augmented rest)

let solve (a : matrix) (b : Vec.t) : Vec.t option =
  solve_augmented (Array.mapi (fun i row -> Array.append row [| b.(i) |]) a)
