(* Where [coeffs] comes from: the quadratic through [y_{n-1}], [y_n] and the
   unknown [y_{n+1}] must satisfy the ODE at [t_{n+1}], which fixes the weights
   from [omega = h / h_prev] alone (docs/numerics/04-bdf.md). *)

open Ode

(* let+ and let* take an Ok value and pass an Error through: docs/ocaml.md. *)
open Fail.Syntax

type coeffs = { a1 : float; a0 : float; beta : float }

let coeffs omega =
  let d = 1. +. (2. *. omega) in
  { a1 = (1. +. omega) *. (1. +. omega) /. d; a0 = -.(omega *. omega) /. d; beta = (1. +. omega) /. d }

(* h_prev is the step just taken, y_prev the state it started from. After carries an inline record: named
   fields without a type of their own (docs/ocaml.md). *)
type history = Start | After of { h_prev : float; y_prev : Vec.t }

let start = Start

(* Only y_0 exists at the start, so the first step is backward Euler. Its local
   error is O(h^2), not O(h^3), but it is made once and the method does not
   amplify it, so the run stays second order. *)
let backward_euler rhs h at =
  let+ y, _ = Bdf1.step rhs h Bdf1.start at in
  y

let bdf2 rhs h ~h_prev ~y_prev at =
  (* Nothing checks omega against the zero-stability limit 1 + sqrt 2: the
     controller must. Halving proposes at most 2; snapping a step of a few ulps
     to the floats can push the taken ratio past it. *)
  let omega = h /. h_prev in
  let { a1; a0; beta } = coeffs omega in
  let psi = Vec.add (Vec.scale a1 at.y) (Vec.scale a0 y_prev) in
  (* Newton starts from the line through the last two points, extended over
     the new step: closer to y_{n+1} than y_n when the solution is smooth. *)
  let guess = Vec.axpy omega (Vec.sub at.y y_prev) at.y in
  Stage.solve rhs { Stage.t = at.t +. h; gamma = beta *. h; psi } guess

let step rhs h history at =
  let+ y =
    match history with
    | Start -> backward_euler rhs h at
    | After { h_prev; y_prev } -> bdf2 rhs h ~h_prev ~y_prev at
  in
  (y, After { h_prev = h; y_prev = at.y })

(* The gap to backward Euler is, to leading order, backward Euler's local error:
   an overestimate for BDF2. The first step is backward Euler itself, so half
   its gap to explicit Euler serves instead: their errors are +-h^2 y'' / 2. *)
let step_with_error rhs h history at =
  let* be = backward_euler rhs h at in
  let next = After { h_prev = h; y_prev = at.y } in
  match history with
  | Start -> Ok (be, Vec.scale 0.5 (Vec.sub be (Vec.axpy h (rhs at.t at.y) at.y)), next)
  | After { h_prev; y_prev } ->
      let+ y = bdf2 rhs h ~h_prev ~y_prev at in
      (y, Vec.sub y be, next)
