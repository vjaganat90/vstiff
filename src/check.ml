let fixed ~dt (p : Ode.problem) =
  if not (dt > 0. && p.t_end -. p.t0 >= 0.) then invalid_arg "Stepper.fixed: need dt > 0 and t_end >= t0"
  else if p.t_end > p.t0 && dt < Clock.resolution (Float.max (Float.abs p.t0) (Float.abs p.t_end)) then
    invalid_arg "Stepper.fixed: dt is below the resolution of t"

let adaptive ~dt0 ~tol (p : Ode.problem) =
  let span = p.t_end -. p.t0 in
  if not (span >= 0.) then invalid_arg "Adaptive.integrate: t_end is before t0"
  else if span > 0. && not (dt0 > 0.) then invalid_arg "Adaptive.integrate: dt0 and dt_max must be positive"
  else if not (tol > 0.) then invalid_arg "Adaptive.integrate: tol must be positive"

let output ~caller (p : Ode.problem) f0 =
  if Array.length f0 <> Array.length p.y0 then invalid_arg (caller ^ ": rhs returns a vector of the wrong length")
