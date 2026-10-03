(** The implicit equation every BDF step reduces to. *)

(** [x = psi + gamma f(t, x)]. *)
type equation = { t : float; gamma : float; psi : Vec.t }

(** [solve rhs eq guess] solves [eq] for [x] by damped Newton from [guess],
    with a forward-difference Jacobian. *)
val solve : Ode.rhs -> equation -> Vec.t -> (Vec.t, Fail.t) result
