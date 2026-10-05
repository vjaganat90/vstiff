(** The gate of the bench: a run is judged against the golden table, which pins the figures of an earlier one.
    Its scd must not be lower and its cost, the right-hand-side calls, not more than 10 % higher. Both are bands, not
    equalities, so they hold from one platform to the next. *)

(** How many percent more right-hand-side calls than the pinned run a row may use. *)
val cost_limit_percent : int

(** How much scd a row may lose against the pinned figure, in digits: the pin is rounded, and other platforms round
    differently. *)
val scd_slack : float

type verdict =
  | Within_band  (** no more costly and no less accurate than the pin, within the band *)
  | Regressed of string  (** the figures that moved, for the reader *)
  | Unpinned  (** the pinned table has no row for this problem and tolerance *)

(** [verdict pinned row] judges [row] against its row of [pinned]. A failed run regresses against a pinned one. *)
val verdict : Measure.t list -> Measure.t -> verdict

(** The text that says what a verdict is. *)
val describe : verdict -> string

(** [source rows] is the OCaml source of golden.ml that pins [rows], which [bench --pin] prints, or [None] when a row
    failed or has a figure that is not finite: that is not a figure to pin. *)
val source : Measure.t list -> string option
