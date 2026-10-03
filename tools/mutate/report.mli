(** The markdown report of a run. Pure. *)

(** How the run went: [baseline] is the time of the unmodified tests, [timeout] the limit each mutant got and [wall]
    the time of the whole run, all in seconds. *)
type run = { workers : int; baseline : float; timeout : float; wall : float }

(** [render ?run entries] counts the outcomes per module, lists each survivor with its id and the expression before and
    after, and lists the mutants known to be equivalent. *)
val render : ?run:run -> (Mutant.t * Outcome.t) list -> string

(** [window a b] shortens two different strings to the part where they differ, with a little text around it, so that
    one changed operator in a long expression stays visible. *)
val window : string -> string -> string * string
