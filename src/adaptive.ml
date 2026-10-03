(** Adaptive BDF2. Every attempt takes the step with backward Euler and with
    BDF2 from the same history; their gap estimates the local error of the
    lower-order member and bounds that of BDF2, which is the one kept.

    Step control: reject and halve [dt] when the estimate exceeds [tol];
    double it, capped at [dt_max], after three accepts in a row. *)

type solution = { t : float; y : Vec.t; accepted : int; rejected : int }

type state = {
  t : float;
  y : Vec.t;
  prev : Bdf2.history option;  (** The last accepted step, once there is one. *)
  dt : float;
  streak : int;  (** Accepts since [dt] last changed. *)
  accepted : int;
  rejected : int;
  failures : int;  (** Consecutive rejections. *)
}

let ( let* ) = Result.bind

(* Largest component of [e], each measured against 1 + |y_i|. *)
let scaled e y = Array.fold_left Float.max 0. (Vec.map2 (fun ei yi -> Float.abs ei /. (1. +. Float.abs yi)) e y)

(* The step's result and its local error estimate. Without history the
   estimate is half the gap between backward Euler and its explicit Euler
   predictor, which is the leading term of backward Euler's local error. *)
let attempt rhs { t; y; prev; _ } h =
  let at = { Ode.t; y } in
  let* be = Bdf1.step rhs h at in
  match prev with
  | None -> Ok (be, Vec.scale 0.5 (Vec.sub be (Vec.axpy h (rhs t y) y)))
  | Some history ->
      let* bdf2 = Bdf2.step rhs h history at in
      Ok (bdf2, Vec.sub bdf2 be)

(* Below this a step can move [t] by only a few units in the last place.
   Needing one means the solution is singular there, or the right-hand side
   fails just ahead, and halving further would never get past it. *)
let dt_min t = 16. *. Float.epsilon *. Float.abs t

(** [integrate ~tol problem] advances [problem.y0] from [t0] to [t_end].
    [dt0] defaults to [1e-6 (t_end - t0)] and [dt_max] to [(t_end - t0) / 10].
    A failed Newton solve counts as a rejection. Returns [Error (StepRejected n)]
    after [n > max_rejects] rejections in a row, or as soon as a halved step
    falls below [16 eps |t|]; [Error Nan] if [y0] or [rhs t0 y0] is not finite.
    @raise Invalid_argument if [t_end < t0] or a step bound is not positive. *)
let integrate ?dt0 ?dt_max ?(max_rejects = 50) ~tol (p : Ode.problem) : (solution, Fail.t) result =
  let { Ode.rhs; t0; t_end; y0 } = p in
  let span = t_end -. t0 in
  let dt_max = Option.value dt_max ~default:(span /. 10.) in
  let dt0 = Float.min dt_max (Option.value dt0 ~default:(1e-6 *. span)) in
  let rec go (s : state) =
    if s.t >= t_end then Ok { t = s.t; y = s.y; accepted = s.accepted; rejected = s.rejected }
    else
      let last = s.dt >= t_end -. s.t in
      let h = if last then t_end -. s.t else s.dt in
      match attempt rhs s h with
      | Ok (y, est) when scaled est y <= tol ->
          let streak = s.streak + 1 in
          go
            {
              t = (if last then t_end else s.t +. h);
              y;
              prev = Some { Bdf2.dt_prev = h; y_prev = s.y };
              dt = (if streak = 3 then Float.min (2. *. s.dt) dt_max else s.dt);
              streak = streak mod 3;
              accepted = s.accepted + 1;
              rejected = s.rejected;
              failures = 0;
            }
      | Ok _ | Error _ ->
          (* Halve the step that failed, which is shorter than [s.dt] when it
             was cut to land on [t_end]. *)
          let failures = s.failures + 1 and dt = h /. 2. in
          if failures > max_rejects || dt < dt_min s.t then Error (Fail.StepRejected failures)
          else go { s with dt; streak = 0; rejected = s.rejected + 1; failures }
  in
  if not (span >= 0.) then invalid_arg "Adaptive.integrate: t_end is before t0"
  else if span > 0. && not (dt0 > 0.) then invalid_arg "Adaptive.integrate: dt0 and dt_max must be positive"
  else if not (Vec.finite y0 && Vec.finite (rhs t0 y0)) then Error Fail.Nan
  else go { t = t0; y = y0; prev = None; dt = dt0; streak = 0; accepted = 0; rejected = 0; failures = 0 }
