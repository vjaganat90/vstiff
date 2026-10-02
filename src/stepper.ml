(** One-step-at-a-time integrators and a fixed-step driver over any of them. *)

type rhs = float -> Vec.t -> Vec.t

module type S = sig
  type state

  val init : t0:float -> Vec.t -> state
  val step : rhs:rhs -> dt:float -> state -> (state, Fail.t) result
  val y : state -> Vec.t
end

(** [fixed (module S) ~rhs ~t0 ~t_end ~dt y0] takes [round ((t_end - t0) / dt)]
    equal steps, so the last one lands on [t_end], and returns the stepper's
    final state or the first failure. *)
let fixed (module S : S) ~rhs ~t0 ~t_end ~dt y0 : (S.state, Fail.t) result =
  let n = max 1 (Float.to_int (Float.round ((t_end -. t0) /. dt))) in
  let h = (t_end -. t0) /. float_of_int n in
  let rec go k s = if k = n then Ok s else Result.bind (S.step ~rhs ~dt:h s) (go (k + 1)) in
  go 0 (S.init ~t0 y0)
