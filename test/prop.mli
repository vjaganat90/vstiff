(** Property-based tests. A property is a module with a generator and a predicate; a check runs it on generated cases
    and shrinks the first case that fails. Checks return values and [Report] prints them. docs/testing.md,
    Properties. *)

(** [holds] must be true on every case that [gen] draws; [show] prints a case for a failure report. *)
module type Property = sig
  type t

  val name : string
  val gen : t Gen.t
  val show : t -> string
  val holds : t -> bool
end

(** A property that must hold on at least the share [threshold] (between 0 and 1) of the cases. *)
module type Statistical = sig
  include Property

  val threshold : float
end

(** What a check found. *)
type outcome =
  | Passed of int  (** Every case held; the number of cases. *)
  | Met of { threshold : float; count : int }  (** At least the share [threshold] of [count] cases held. *)
  | Failed of { seed : int; case : string; shrinks : int; why : string }
      (** The run seed, the failing case as [show] prints it after [shrinks] shrink steps, and [why]: which case
          failed first, or how many held, and the exception if the shrunk case raises. *)

(** [check (module P) ~count ~seed] runs [P] on [count] cases, drawn from a state made from [seed] and [P.name],
    and stops at the first that fails. A case that raises fails. *)
val check : (module P : Property) -> count:int -> seed:int -> outcome

(** As [check], for a statistical property: it runs every case, and fails when too few held. *)
val check_statistical : (module P : Statistical) -> count:int -> seed:int -> outcome

(** [ok 1000], [ok, at least 99% of 1000], or [FAILED] with why, the seed and the shrunk case. *)
val to_string : outcome -> string
