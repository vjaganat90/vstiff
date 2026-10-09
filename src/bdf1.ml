(* open Ode makes at.t mean Ode.point's t; Stage.t names a record type; history is unit: docs/ocaml.md, Records. *)
open Ode

type history = unit

let start = ()

let step rhs h () at =
  Result.map (fun y -> (y, ())) (Stage.solve rhs { Stage.t = at.t +. h; gamma = h; psi = at.y } at.y)
