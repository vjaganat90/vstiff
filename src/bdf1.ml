(** Backward Euler: [y_{n+1} = y_n + dt f(t_{n+1}, y_{n+1})]. *)

let step ~rhs ~t ~dt y = Stage.solve ~rhs ~t:(t +. dt) ~gamma:dt ~psi:y y

module Fixed = struct
  type state = { t : float; y : Vec.t }

  let backward_euler = step
  let init ~t0 y = { t = t0; y }
  let step ~rhs ~dt { t; y } = Result.map (fun y -> { t = t +. dt; y }) (backward_euler ~rhs ~t ~dt y)
  let y s = s.y
end

let integrate ~rhs ~t0 ~t_end ~dt y0 =
  Result.map Fixed.y (Stepper.fixed (module Fixed) ~rhs ~t0 ~t_end ~dt y0)
