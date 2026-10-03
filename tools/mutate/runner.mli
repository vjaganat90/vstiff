(** Runs a repository's tests against mutants. Each worker has a scratch copy of the repository outside it; a mutant
    replaces one file of the copy, is type-checked (does not compile: stillborn), then built and tested with
    [dune build @runtest] (a failure or the timeout: killed; none: survived), and the file is put back. *)

(** [root] is the repository, [workers] the number of copies and mutants run at once, [timeout_factor] the multiple of
    the baseline run's time that a mutant may take, and [keep] leaves the copies in place. *)
type config = { root : string; workers : int; timeout_factor : float; keep : bool }

(** Seconds: the unmodified tests' run, the timeout derived from it, and the whole run. *)
type timing = { baseline : float; timeout : float; wall : float }

(** [run config mutants say on_outcome] copies the repository, checks that the unmodified copies pass their tests and
    that printing each mutated file unchanged does too, then runs the mutants, handing every outcome to [on_outcome]
    as soon as it is known and progress messages to [say]. [Error] when the unmodified tests fail, printing changes
    what the code does, or the run is interrupted (SIGINT or SIGTERM; nothing is left behind). *)
val run : config -> Mutant.t list -> (string -> unit) -> (Mutant.t -> Outcome.t -> unit) -> (timing, string) result
