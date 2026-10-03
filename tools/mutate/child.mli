(** Child processes and the clock, the tool's process effects. Each child leads a process group of its own, so killing
    it takes along what it started (dune starts the tests). Nothing blocks: one loop can watch several children. *)

type t

(** [spawn ~dir ~log argv] starts the command [argv] (looked up in [PATH]) in the directory [dir], with its output and
    errors written to the file [log], which is emptied first. *)
val spawn : dir:string -> log:string -> string list -> t

(** The exit status once the child has ended, [None] while it runs. *)
val poll : t -> Unix.process_status option

(** Seconds since [spawn], or the child's whole run time once it has ended. *)
val elapsed : t -> float

(** Kills the child and its group, and waits for it. *)
val kill : t -> unit

(** Waits a few milliseconds between two polls. *)
val pause : unit -> unit

(** The wall clock, in seconds. *)
val now : unit -> float
