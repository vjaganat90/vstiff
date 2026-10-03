(** Variable-step BDF2. With [omega = h / h_prev] the step solves
    [y_{n+1} = a1 y_n + a0 y_{n-1} + beta h f(t_{n+1}, y_{n+1})], the
    coefficients making the formula exact on quadratics through
    [t_{n-1}], [t_n], [t_{n+1}]. *)

open Ode
open Fail.Syntax

type coeffs = { a1 : float; a0 : float; beta : float }

let coeffs omega =
  let d = 1. +. (2. *. omega) in
  { a1 = (1. +. omega) *. (1. +. omega) /. d; a0 = -.(omega *. omega) /. d; beta = (1. +. omega) /. d }

(* Before the first step there is no previous point, so that step is
   backward Euler. *)
type history = Start | After of { h_prev : float; y_prev : Vec.t }

let start = Start

let backward_euler rhs h at =
  let+ y, _ = Bdf1.step rhs h Bdf1.start at in
  y

let bdf2 rhs h ~h_prev ~y_prev at =
  let omega = h /. h_prev in
  let { a1; a0; beta } = coeffs omega in
  let psi = Vec.add (Vec.scale a1 at.y) (Vec.scale a0 y_prev) in
  (* Newton starts from the line through the last two points. *)
  let guess = Vec.axpy omega (Vec.sub at.y y_prev) at.y in
  Stage.solve rhs { Stage.t = at.t +. h; gamma = beta *. h; psi } guess

let step rhs h history at =
  let+ y =
    match history with
    | Start -> backward_euler rhs h at
    | After { h_prev; y_prev } -> bdf2 rhs h ~h_prev ~y_prev at
  in
  (y, After { h_prev = h; y_prev = at.y })

(* The estimate is the gap to backward Euler from the same point. On the
   first step it is half the gap between backward Euler and its explicit
   Euler predictor, the leading term of backward Euler's local error. *)
let step_with_error rhs h history at =
  let* be = backward_euler rhs h at in
  let next = After { h_prev = h; y_prev = at.y } in
  match history with
  | Start -> Ok (be, Vec.scale 0.5 (Vec.sub be (Vec.axpy h (rhs at.t at.y) at.y)), next)
  | After { h_prev; y_prev } ->
      let+ y = bdf2 rhs h ~h_prev ~y_prev at in
      (y, Vec.sub y be, next)
