(** Why a solve or an integration stopped. Numerical failures are values,
    never exceptions. *)

type t =
  | Diverged  (** Newton found no root: singular Jacobian, no descent, or out of iterations. *)
  | StepRejected of int  (** Step control gave up after this many rejections in a row. *)
  | Nan  (** The state or the right-hand side stopped being finite. *)

val to_string : t -> string
