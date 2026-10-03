(** Variable-step BDF2. With [omega = dt / dt_prev] the step solves
    [y_{n+1} = a1 y_n + a0 y_{n-1} + beta dt f(t_{n+1}, y_{n+1})], the
    coefficients making the formula exact on quadratics through
    [t_{n-1}], [t_n], [t_{n+1}]. *)

type coeffs = { a1 : float; a0 : float; beta : float }

let coeffs omega =
  let d = 1. +. (2. *. omega) in
  { a1 = (1. +. omega) *. (1. +. omega) /. d; a0 = -.(omega *. omega) /. d; beta = (1. +. omega) /. d }

type history = { dt_prev : float; y_prev : Vec.t }

let step rhs dt { dt_prev; y_prev } (at : Ode.point) =
  let omega = dt /. dt_prev in
  let { a1; a0; beta } = coeffs omega in
  let psi = Vec.add (Vec.scale a1 at.y) (Vec.scale a0 y_prev) in
  (* Newton starts from the line through the last two points. *)
  let guess = Vec.axpy omega (Vec.sub at.y y_prev) at.y in
  Stage.solve rhs { Stage.t = at.t +. dt; gamma = beta *. dt; psi } guess

(* Fixed steps; the first has no history and is backward Euler. *)
let integrate ~dt (p : Ode.problem) =
  let advance h ((at : Ode.point), prev) =
    let next = match prev with None -> Bdf1.step p.rhs h at | Some history -> step p.rhs h history at in
    Result.map (fun y -> ({ Ode.t = at.t +. h; y }, Some { dt_prev = h; y_prev = at.y })) next
  in
  Result.map (fun ((at : Ode.point), _) -> at.y) (Stepper.fixed ~dt p advance ({ Ode.t = p.t0; y = p.y0 }, None))
