(* Newton solves G(x) = x - psi - gamma f(t, x) = 0. The Jacobian of G is
   I - gamma J, J being the Jacobian of f. Without gamma J the update would be
   x <- psi + gamma f(t, x), fixed-point iteration, which converges only for
   short steps; stiff problems need long ones. See docs/numerics/02-newton.md, section 7. *)

type equation = { t : float; gamma : float; psi : Vec.t }

let solve (rhs : Ode.rhs) { t; gamma; psi } guess =
  (* [rhs t] is [f] with the time fixed: partial application, docs/ocaml.md. *)
  let f = rhs t in
  let residual x = Vec.sub (Vec.sub x psi) (Vec.scale gamma (f x)) in
  let jacobian x =
    Array.mapi
      (fun i row -> Array.mapi (fun j v -> (if i = j then 1. else 0.) -. (gamma *. v)) row)
      (Jac.forward f x)
  in
  Newton.solve residual jacobian guess
