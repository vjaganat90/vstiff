(** Newton's method for [f x = 0], with a line search that damps steps which would overshoot.

    Each iteration solves [jac(x) dx = -f(x)] for [dx] and moves by [lambda dx], where [lambda] is the first of
    1, 1/2, ..., 1/1024 that cuts the residual norm to [1 - 1e-4 lambda] times its old value or less. *)

(** [solve f jac x0] solves [f x = 0] from [x0]; [jac x] is the Jacobian of [f] at [x]. Returns [Ok (x + dx)]
    once [|dx|_inf <= 1e-10 (1 + |x|_inf)] (or [Ok x] if [f x] is exactly 0), else [Error Diverged] (singular
    Jacobian, non-finite step, no [lambda] passes, or 50 iterations) or [Error Nan] (non-finite [x], [f x] or
    [x + dx]). *)
val solve : (Vec.t -> Vec.t) -> (Vec.t -> Linalg.matrix) -> Vec.t -> (Vec.t, Fail.t) result
