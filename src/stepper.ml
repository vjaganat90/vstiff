(** A fixed-step driver for any one-step method. *)

(** [fixed ~dt problem advance start] applies [advance h] to [start]
    [round ((t_end - t0) / dt)] times (at least once when [t_end > t0]) with
    [h] chosen so the steps land exactly on [t_end], stopping at the first
    failure.
    @raise Invalid_argument unless [dt > 0] and [t_end >= t0]. *)
let fixed ~dt (p : Ode.problem) advance start =
  let span = p.t_end -. p.t0 in
  if not (dt > 0. && span >= 0.) then invalid_arg "Stepper.fixed: need dt > 0 and t_end >= t0"
  else
    let n = if span = 0. then 0 else max 1 (Float.to_int (Float.round (span /. dt))) in
    let h = span /. float_of_int (max n 1) in
    let rec go k s = if k = n then Ok s else Result.bind (advance h s) (go (k + 1)) in
    go 0 start
