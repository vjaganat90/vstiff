(** Counts calls by wrapping the [rhs] a method is given, so methods stay pure; the library's only mutable state. *)

(** [count rhs] is [rhs] wrapped to count its calls, and a function that reads the count. *)
val count : Ode.rhs -> Ode.rhs * (unit -> int)
