(** The tool's file system effects, kept in one module so that the rest stays pure. *)

(** The whole contents of a file. *)
val read : string -> string

val exists : string -> bool

(** [write path ~contents] replaces the file, [append] adds to its end; both create it if it does not exist. *)
val write : string -> contents:string -> unit

val append : string -> contents:string -> unit

(** [copy_tree ~src ~dst] copies the directory [src] to a new directory [dst], leaving out what dune ignores: the
    entries whose names start with a dot or an underscore ([.git], [_build]). *)
val copy_tree : src:string -> dst:string -> unit

(** [temp_dir prefix] creates an empty directory with a name that starts with [prefix] in the system's temporary
    directory ([TMPDIR]), and returns its path. *)
val temp_dir : string -> string

(** Deletes a file or a directory with everything in it. A symbolic link is removed, never followed. *)
val remove_tree : string -> unit
