(* Gaussian elimination as a recursion (let rec: docs/ocaml.md) on fresh arrays, so nothing is mutated.
   Worked example and the pivoting argument: docs/numerics/02-newton.md. *)

type matrix = float array array

(* Partial pivoting: the pivot row has the largest |a.(i).(0)|, the first one on a tie. A zero top entry
   cannot be divided by, and a tiny one gives huge multipliers that round-off turns into a wrong answer. *)
let pivot (a : matrix) =
  let rec best i j =
    if j = Array.length a then i else best (if Float.abs a.(j).(0) > Float.abs a.(i).(0) then j else i) (j + 1)
  in
  best 0 1

(* Which original row sits at position [i] once rows 0 and [p] trade places. *)
let swapped p i = if i = 0 then p else if i = p then 0 else i

(* [rows] are augmented, [A | b]: coefficients, then the right-hand side. Exchanging rows reorders equations,
   not unknowns, so there is nothing to undo at the end. *)
let rec solve_augmented (rows : matrix) : Vec.t option =
  match Array.length rows with
  | 0 -> Some [||]
  | n ->
      let p = pivot rows in
      let top = rows.(p) in
      (* The pivot has the largest absolute value in its column, so zero means the whole column is zero. *)
      if top.(0) = 0. then None
      else
        let reduce row =
          let m = row.(0) /. top.(0) in
          Array.init n (fun j -> row.(j + 1) -. (m *. top.(j + 1)))
        in
        let rest = Array.init (n - 1) (fun i -> reduce rows.(swapped p (i + 1))) in
        (* Option.map: f inside Some, None passes through (docs/ocaml.md); a deeper zero pivot surfaces as None. *)
        Option.map
          (fun x -> Array.append [| (top.(n) -. Vec.dot (Array.sub top 1 (n - 1)) x) /. top.(0) |] x)
          (solve_augmented rest)

let solve (a : matrix) (b : Vec.t) : Vec.t option =
  solve_augmented (Array.mapi (fun i row -> Array.append row [| b.(i) |]) a)
