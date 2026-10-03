(** Why a solve or an integration stopped. Numerical failures are returned as [Error], never raised. *)

type t =
  | Diverged  (** [Newton.solve] gave up: no usable step, no residual decrease, or out of iterations. *)
  | StepRejected of int  (** [Halving] gave up after this many rejections in a row. *)
  | Nan  (** [Newton.solve] or [Adaptive.integrate] found a value that had to be finite and was not. *)

(** The constructor name, followed by the count for [StepRejected]: [Diverged], [StepRejected 46], [Nan]. *)
val to_string : t -> string

(** Chain results, stopping at the first [Error]: [let* x = r in e] is [Result.bind r (fun x -> e)],
    and [let+ x = r in e] is [Result.map (fun x -> e) r]. *)
module Syntax : sig
  val ( let* ) : ('a, 'e) result -> ('a -> ('b, 'e) result) -> ('b, 'e) result
  val ( let+ ) : ('a, 'e) result -> ('a -> 'b) -> ('b, 'e) result
end
