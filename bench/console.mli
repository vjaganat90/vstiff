(** The only module of the bench with effects: it reads the command line, prints, reads the CPU clock and exits.
    Everything else is pure, as the library is. *)

(** The arguments after the program name. *)
val args : unit -> string list

(** [print s] writes [s] and a newline on standard output and flushes it; [warn] does the same on standard error. *)
val print : string -> unit

val warn : string -> unit

(** [cpu f] is [f ()] and the CPU seconds it took ([Sys.time]), which are not reproducible. *)
val cpu : (unit -> 'a) -> 'a * float

(** [exit status] ends the program. *)
val exit : int -> 'a
