(** Dense linear solves. *)

(** Row-major: [a.(i).(j)] is row [i], column [j]. *)
type matrix = float array array

(** [solve a b] is [Some x] with [a x = b], by Gaussian elimination with
    partial pivoting, or [None] when a pivot is exactly zero. *)
val solve : matrix -> Vec.t -> Vec.t option
