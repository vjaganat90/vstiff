(** The implicit stage equation every BDF step reduces to:
    [x = psi + gamma f(t, x)], solved by damped Newton on
    [G x = x - psi - gamma f(t, x)] with [G' = I - gamma J],
    [J] by forward differences. *)

let solve ~(rhs : Stepper.rhs) ~t ~gamma ~psi guess =
  let f = rhs t in
  let g x = Vec.sub (Vec.sub x psi) (Vec.scale gamma (f x)) in
  let jac x =
    Array.mapi
      (fun i row -> Array.mapi (fun j v -> (if i = j then 1. else 0.) -. (gamma *. v)) row)
      (Jac.forward f x)
  in
  Newton.solve ~f:g ~jac guess
