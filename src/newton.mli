(** Damped Newton for [f x = 0]. *)

type config = {
  tol : float;  (** Converged once [|dx|_inf <= tol (1 + |x|_inf)]. *)
  max_iter : int;
  min_damping : float;  (** Smallest line-search factor tried before giving up. *)
}

(** [tol = 1e-10], [max_iter = 50], [min_damping = 1/1024]. *)
val default : config

(** [solve f jac x0] starts at [x0] and returns the root, [Error Diverged]
    (singular Jacobian, no residual decrease, or out of iterations) or
    [Error Nan] (non-finite iterate or residual). *)
val solve : ?config:config -> (Vec.t -> Vec.t) -> (Vec.t -> Linalg.matrix) -> Vec.t -> (Vec.t, Fail.t) result
