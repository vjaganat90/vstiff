(** The implicit equation every BDF step reduces to (docs/numerics/04-bdf.md). *)

(** [x = psi + gamma f(t, x)] for the new state [x]: [t] is the time of the new
    point, [psi] comes from the past, [gamma] weighs the new slope. *)
type equation = { t : float; gamma : float; psi : Vec.t }

(** [solve rhs eq guess] solves [eq] for [x] by damped Newton from [guess], with
    a forward-difference Jacobian; a [guess] near the root saves iterations.
    Fails as [Newton.solve] does: [Error Diverged] or [Error Nan]. *)
val solve : Ode.rhs -> equation -> Vec.t -> (Vec.t, Fail.t) result
