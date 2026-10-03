(** Vectors as [float array]s, never mutated, so sharing is safe. [add], [sub], [axpy] and [dot] use the first
    vector's length: a shorter second one raises [Invalid_argument], extra entries of a longer one are ignored. *)

type t = float array

val add : t -> t -> t
val sub : t -> t -> t
val scale : float -> t -> t

(** [axpy a x y] is [a x + y]. *)
val axpy : float -> t -> t -> t

val dot : t -> t -> float

(** Largest absolute entry, [|v|_inf]; [nan] if any entry is [nan], so a [nan] fails any [<= tol] test. *)
val norm_inf : t -> float

(** Every entry is neither infinite nor [nan]. *)
val finite : t -> bool
