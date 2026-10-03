(** The tool's file system effects, kept in one module so that the rest stays pure. *)

(** The whole contents of a file. *)
val read : string -> string
