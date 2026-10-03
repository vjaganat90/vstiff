(** What running the tests against one mutant showed, and the file the runner keeps its outcomes in. Pure. *)

type t =
  | Killed of string
      (** A test failed or crashed, or the run hit the timeout. The text says where: a file, an exception, [timeout]. *)
  | Survived  (** Every test passed. *)
  | Stillborn of string  (** The mutant does not compile. The text is the compiler's first error line. *)
  | Equivalent of string  (** Listed as equivalent to the original, with the reason; not run. *)

(** A mutant that ran past the time limit. *)
val timeout : t

(** [killed], [survived], [stillborn] or [equivalent]. *)
val name : t -> string

(** The text of [Killed], [Stillborn] and [Equivalent]; empty for [Survived]. *)
val detail : t -> string

(** [killed ~log] is [Killed] with where the failing build log [log] points: an uncaught exception, else the file of
    the first [File] line (a diff names the [.expected] file), else [failed]. *)
val killed : log:string -> t

(** [stillborn ~log] is [Stillborn] with the first [Error] line of the compiler's [log]. *)
val stillborn : log:string -> t

(** One line of a results file: tab-separated file, line, column, operator, outcome, detail, original and mutated,
    with backslash, tab and newline escaped. No newline at the end. *)
val to_line : Mutant.t * t -> string

(** Reads a line written by [to_line]; [None] for any other line. *)
val of_line : string -> (Mutant.t * t) option

(** [equivalents text] reads a list of known equivalent mutants: one per line, the id and then the reason; blank lines
    and lines starting with [#] are skipped. *)
val equivalents : string -> (string * string) list
