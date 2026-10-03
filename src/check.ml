let fixed ~dt (p : Ode.problem) =
  if not (dt > 0. && p.t_end -. p.t0 >= 0.) then invalid_arg "Stepper.fixed: need dt > 0 and t_end >= t0"

let adaptive ~dt0 (p : Ode.problem) =
  let span = p.t_end -. p.t0 in
  if not (span >= 0.) then invalid_arg "Adaptive.integrate: t_end is before t0"
  else if span > 0. && not (dt0 > 0.) then invalid_arg "Adaptive.integrate: dt0 and dt_max must be positive"
