(** Dense linear solves. *)

(** An array of rows: [a.(i).(j)] is row [i], column [j]. *)
type matrix = float array array

(** [solve a b] is [Some x] with [a x = b] (Gaussian elimination, partial pivoting), or [None] when a pivot is
    exactly zero. A singular [a] gives [None] unless rounding hides it; a nearly singular one is not detected, and
    [x] may be huge or non-finite. [a] must be square, [b] as long as [a] (not checked). *)
val solve : matrix -> Vec.t -> Vec.t option
