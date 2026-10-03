(** Backward Euler: [y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})]. *)

open Ode

type history = unit

let start = ()

let step rhs h () at =
  Result.map (fun y -> (y, ())) (Stage.solve rhs { Stage.t = at.t +. h; gamma = h; psi = at.y } at.y)
