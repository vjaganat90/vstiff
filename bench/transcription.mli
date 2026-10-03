(** The right-hand sides of the bench problems at seeded points, for the transcription check (plan 7): every problem
    exists in OCaml and in Python, and before any comparison the two right-hand sides must agree. The Python side,
    bench/compare/transcription.py, reads these lines and evaluates its own. *)

(** The number of points per problem. *)
val points : int

(** One line per point and problem: [problem, t, y_1, ..., y_n, f_1, ..., f_n], with [f = rhs t y] and every number
    at 17 digits, which a double survives. The points are the same on every platform: the generator is pure integer
    arithmetic. *)
val lines : unit -> string list
