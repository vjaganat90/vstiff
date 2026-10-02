(** Why an integration or a nonlinear solve stopped. *)
type t =
  | Diverged  (** Newton found no root: singular Jacobian, no descent, or out of iterations. *)
  | StepRejected of int  (** Step control gave up after this many consecutive rejections. *)
  | Nan  (** The state or the right-hand side stopped being finite. *)

let to_string = function
  | Diverged -> "Diverged"
  | StepRejected n -> Printf.sprintf "StepRejected %d" n
  | Nan -> "Nan"
