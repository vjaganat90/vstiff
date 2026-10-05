(** One run of vstiff's adaptive driver, measured. *)

(** What a successful run costs and how accurate it is. *)
type figures = {
  steps : int;  (** accepted steps *)
  rejected : int;  (** rejected steps, whatever the reason *)
  rhs_calls : int;  (** every call of the right-hand side, counted through [Vstiff.Instrument] *)
  error : float;  (** [max_i |y_i - ref_i| / (1 + |ref_i|)] at [t_end]: the mixed error *)
  scd : float;  (** [-log10 max_i |y_i - ref_i| / |ref_i|]: the correct digits, as in the IVP Test Set *)
}

(** A run finishes with its figures, or fails with the error the driver returns. *)
type outcome = Done of figures | Failed of Vstiff.Fail.t

(** The run of [problem] at [tol]. *)
type t = { problem : string; tol : float; outcome : outcome }

(** [error ~reference y] is the mixed error of [y]; [reference] and [y] have the same length. *)
val error : reference:float array -> float array -> float

(** [scd ~reference y] is the number of correct digits of [y], [infinity] for an exact answer. A zero component of
    [reference] counts as an absolute difference. *)
val scd : reference:float array -> float array -> float

(** [run case ~tol] integrates [case.problem] with [Bdf2] and [Halving] at [rtol = atol = tol], every other argument
    left at its default, and measures the run against [case.reference]. *)
val run : Suite.case -> tol:float -> t
