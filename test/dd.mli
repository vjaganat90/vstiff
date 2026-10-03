(** Double-double arithmetic, the oracle of properties that must see past float rounding: a number is the unevaluated
    sum [hi + lo] of two floats, about 106 bits. Sums and products of a few finite floats come out within a few
    [eps^2] of their magnitude, as long as nothing overflows or underflows. *)

type t = { hi : float; lo : float }

val of_float : float -> t

(** The nearest float, [hi]. *)
val to_float : t -> float

val add : t -> t -> t
val sub : t -> t -> t
val mul : t -> t -> t

(** [dot a b] is the sum of [a.(i) *. b.(i)] over the indices of [a], every product exact. *)
val dot : float array -> float array -> t
