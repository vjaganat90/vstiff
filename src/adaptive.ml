(** Adaptive BDF2. Every attempt takes the step with backward Euler and with
    BDF2 from the same history; their gap estimates the local error of the
    lower-order member and bounds that of BDF2, which is the one kept.

    Step control: reject and halve [dt] when the estimate exceeds [tol];
    double it, capped at [dt_max], after three accepts in a row. *)

type solution = { t : float; y : Vec.t; accepted : int; rejected : int }

type state = {
  t : float;
  y : Vec.t;
  prev : (float * Vec.t) option;  (** Previous accepted step and the state it started from. *)
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
let attempt ~rhs { t; y; prev; _ } h =
  let* be = Bdf1.step ~rhs ~t ~dt:h y in
  match prev with
  | None -> Ok (be, Vec.scale 0.5 (Vec.sub be (Vec.axpy h (rhs t y) y)))
  | Some (dt_prev, y_prev) ->
      let* bdf2 = Bdf2.step ~rhs ~t ~dt:h ~dt_prev ~y_prev y in
      Ok (bdf2, Vec.sub bdf2 be)

let integrate ?dt0 ?dt_max ?(max_rejects = 50) ~tol ~rhs ~t0 ~t_end y0 : (solution, Fail.t) result =
  let span = t_end -. t0 in
  let dt_max = Option.value dt_max ~default:(span /. 10.) in
  let rec go (s : state) =
    if s.t >= t_end then Ok { t = s.t; y = s.y; accepted = s.accepted; rejected = s.rejected }
    else if s.failures > max_rejects then Error (Fail.StepRejected s.failures)
    else
      let last = s.dt >= t_end -. s.t in
      let h = if last then t_end -. s.t else s.dt in
      match attempt ~rhs s h with
      | Ok (y, est) when scaled est y <= tol ->
          let streak = s.streak + 1 in
          go
            {
              t = (if last then t_end else s.t +. h);
              y;
              prev = Some (h, s.y);
              dt = (if streak = 3 then Float.min (2. *. s.dt) dt_max else s.dt);
              streak = streak mod 3;
              accepted = s.accepted + 1;
              rejected = s.rejected;
              failures = 0;
            }
      | Ok _ | Error _ ->
          go { s with dt = s.dt /. 2.; streak = 0; rejected = s.rejected + 1; failures = s.failures + 1 }
  in
  if not (Vec.finite y0 && Vec.finite (rhs t0 y0)) then Error Fail.Nan
  else
    go
      {
        t = t0;
        y = y0;
        prev = None;
        dt = Float.min dt_max (Option.value dt0 ~default:(1e-6 *. span));
        streak = 0;
        accepted = 0;
        rejected = 0;
        failures = 0;
      }
