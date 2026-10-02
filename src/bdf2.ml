(** Variable-step BDF2. With [omega = dt / dt_prev] the step solves
    [y_{n+1} = a1 y_n + a0 y_{n-1} + beta dt f(t_{n+1}, y_{n+1})], the
    coefficients making the formula exact on quadratics through
    [t_{n-1}], [t_n], [t_{n+1}]. *)

type coeffs = { a1 : float; a0 : float; beta : float }

let coeffs omega =
  let d = 1. +. (2. *. omega) in
  { a1 = (1. +. omega) *. (1. +. omega) /. d; a0 = -.(omega *. omega) /. d; beta = (1. +. omega) /. d }

let step ~rhs ~t ~dt ~dt_prev ~y_prev y =
  let omega = dt /. dt_prev in
  let { a1; a0; beta } = coeffs omega in
  let psi = Vec.add (Vec.scale a1 y) (Vec.scale a0 y_prev) in
  (* Newton starts from the line through the last two points. *)
  let guess = Vec.axpy omega (Vec.sub y y_prev) y in
  Stage.solve ~rhs ~t:(t +. dt) ~gamma:(beta *. dt) ~psi guess

(** Fixed steps; the first has no history and is backward Euler. *)
module Fixed = struct
  type state = { t : float; y : Vec.t; prev : (float * Vec.t) option }

  let bdf2 = step
  let init ~t0 y = { t = t0; y; prev = None }

  let step ~rhs ~dt { t; y; prev } =
    let next =
      match prev with
      | None -> Bdf1.step ~rhs ~t ~dt y
      | Some (dt_prev, y_prev) -> bdf2 ~rhs ~t ~dt ~dt_prev ~y_prev y
    in
    Result.map (fun y' -> { t = t +. dt; y = y'; prev = Some (dt, y) }) next

  let y s = s.y
end

let integrate ~rhs ~t0 ~t_end ~dt y0 =
  Result.map Fixed.y (Stepper.fixed (module Fixed) ~rhs ~t0 ~t_end ~dt y0)
