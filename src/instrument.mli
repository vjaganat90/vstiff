(** Instrumentation by wrapping what a method is given, so that methods stay
    pure. This is the only module in the library with mutable state. *)

(** [count rhs] is [rhs] wrapped to count its calls, and a function that reads
    the count. *)
val count : Ode.rhs -> Ode.rhs * (unit -> int)
