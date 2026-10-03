(** The command line of the bench. *)

type mode =
  | Table  (** run every row and print it, with its verdict against the golden table *)
  | Csv  (** run every row and print it as CSV, the columns of bench/compare/scipy_rows.py *)
  | Pin  (** run every row and print the source of golden.ml that pins them *)
  | Transcription  (** run nothing; print the right-hand sides at seeded points, for bench/compare/transcription.py *)
  | Help  (** print {!usage} *)

type t = {
  mode : mode;
  cpu : bool;  (** add the CPU time of each row, which is not reproducible *)
  only : string option;  (** run only the problem of this name *)
}

(** [parse args] reads the arguments after the program name; [Error] carries a message for the user. [Pin] and
    [Transcription] cover every problem and take neither [cpu] nor [only]: a pin of some rows would drop the others. *)
val parse : string list -> (t, string) result

val usage : string
