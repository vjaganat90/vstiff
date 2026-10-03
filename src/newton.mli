(** Damped Newton for [f x = 0]. *)

(** [solve f jac x0] starts at [x0] and returns the root, [Error Diverged]
    (singular Jacobian, no residual decrease, or out of iterations) or
    [Error Nan] (non-finite iterate or residual). Converged once
    [|dx|_inf <= 1e-10 (1 + |x|_inf)]; at most 50 iterations; the line search
    halves the step down to [1/1024]. *)
val solve : (Vec.t -> Vec.t) -> (Vec.t -> Linalg.matrix) -> Vec.t -> (Vec.t, Fail.t) result
