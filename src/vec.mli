(** Vectors as [float array]s that are never mutated: every operation
    allocates its result. Binary operations expect equal lengths. *)

type t = float array

val add : t -> t -> t
val sub : t -> t -> t
val scale : float -> t -> t

(** [axpy a x y] is [a x + y]. *)
val axpy : float -> t -> t -> t

val dot : t -> t -> float

(** Largest absolute entry; [nan] if any entry is [nan]. *)
val norm_inf : t -> float

(** Every entry is neither infinite nor [nan]. *)
val finite : t -> bool
