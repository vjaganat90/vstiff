(** The text the bench prints, as strings: nothing here prints. *)

(** One run, with its CPU time when it was asked for, and its verdict against the golden table. *)
type entry = { row : Measure.t; cpu : float option; verdict : Gate.verdict }

(** What the rows of the table are, as one line above the header. *)
val title : string

(** The header of the table, with a column for the CPU time if the options ask for it, and one line per entry, which
    ends with the verdict. *)
val table_header : Options.t -> string

val table_line : entry -> string

(** The header of the CSV, and one line per entry. The columns are those of bench/compare/scipy_rows.py: [problem,
    solver, rtol, steps, rejected, rhs, nfev, njev, error, scd], with [nfev] and [njev] empty, then [cpu_s] when the
    options ask for the CPU time. A failed run has empty figures. *)
val csv_header : Options.t -> string

val csv_line : entry -> string
