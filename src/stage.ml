(** The implicit stage equation every BDF step reduces to:
    [x = psi + gamma f(t, x)], solved by damped Newton on
    [G x = x - psi - gamma f(t, x)] with [G' = I - gamma J],
    [J] by forward differences. *)

type equation = { t : float; gamma : float; psi : Vec.t }

let solve (rhs : Ode.rhs) { t; gamma; psi } guess =
  let f = rhs t in
  let residual x = Vec.sub (Vec.sub x psi) (Vec.scale gamma (f x)) in
  let jacobian x =
    Array.mapi
      (fun i row -> Array.mapi (fun j v -> (if i = j then 1. else 0.) -. (gamma *. v)) row)
      (Jac.forward f x)
  in
  Newton.solve residual jacobian guess
