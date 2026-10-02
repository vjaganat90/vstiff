(** One-step-at-a-time integrators and a fixed-step driver over any of them. *)

type rhs = float -> Vec.t -> Vec.t

module type S = sig
  type state

  val init : t0:float -> Vec.t -> state
  val step : rhs:rhs -> dt:float -> state -> (state, Fail.t) result
  val y : state -> Vec.t
end

(** [fixed (module S) ~rhs ~t0 ~t_end ~dt y0] takes [round ((t_end - t0) / dt)]
    equal steps (at least one when [t_end > t0]), so the last one lands on
    [t_end], and returns the stepper's final state or the first failure.
    @raise Invalid_argument unless [dt > 0] and [t_end >= t0]. *)
let fixed (module S : S) ~rhs ~t0 ~t_end ~dt y0 : (S.state, Fail.t) result =
  let span = t_end -. t0 in
  if not (dt > 0. && span >= 0.) then invalid_arg "Stepper.fixed: need dt > 0 and t_end >= t0"
  else
    let n = if span = 0. then 0 else max 1 (Float.to_int (Float.round (span /. dt))) in
    let h = span /. float_of_int (max n 1) in
    let rec go k s = if k = n then Ok s else Result.bind (S.step ~rhs ~dt:h s) (go (k + 1)) in
    go 0 (S.init ~t0 y0)
