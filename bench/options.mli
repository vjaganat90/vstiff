(** The command line of the bench. *)

type mode =
  | Table  (** run every row and print it *)
  | Csv  (** run every row and print it as CSV, the columns of bench/compare/scipy_rows.py *)
  | Help  (** print {!usage} *)

type t = {
  mode : mode;
  cpu : bool;  (** add the CPU time of each row, which is not reproducible *)
  only : string option;  (** run only the problem of this name *)
}

(** [parse args] reads the arguments after the program name; [Error] carries a message for the user. *)
val parse : string list -> (t, string) result

val usage : string
