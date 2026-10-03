(** Backward Euler: [y_{n+1} = y_n + dt f(t_{n+1}, y_{n+1})]. *)

let step rhs dt (at : Ode.point) = Stage.solve rhs { Stage.t = at.t +. dt; gamma = dt; psi = at.y } at.y

let integrate ~dt (p : Ode.problem) =
  let advance h (at : Ode.point) = Result.map (fun y -> { Ode.t = at.t +. h; y }) (step p.rhs h at) in
  Result.map (fun (at : Ode.point) -> at.y) (Stepper.fixed ~dt p advance { Ode.t = p.t0; y = p.y0 })
