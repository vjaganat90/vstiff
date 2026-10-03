(** Adaptive BDF2. Every attempt takes the step with backward Euler and with
    BDF2 from the same history; their gap estimates the local error of the
    lower-order member and bounds that of BDF2, which is the one kept.

    Step control: reject and halve the step when the estimate exceeds [tol];
    double it, capped at [dt_max], after three accepts in a row. *)

type solution = { t : float; y : Vec.t; accepted_steps : int; rejected_steps : int }

type state = {
  at : Ode.point;
  prev : Bdf2.history option;  (** The last accepted step, once there is one. *)
  dt : float;  (** The next step to try. *)
  streak : int;  (** Accepts since [dt] last changed. *)
  failures : int;  (** Rejections in a row. *)
  accepted_steps : int;
  rejected_steps : int;
}

let ( let* ) = Result.bind

let solution { at; accepted_steps; rejected_steps; _ } = { t = at.t; y = at.y; accepted_steps; rejected_steps }

(* Largest component of [e], each measured against 1 + |y_i|. *)
let scaled e y = Array.fold_left Float.max 0. (Array.map2 (fun ei yi -> Float.abs ei /. (1. +. Float.abs yi)) e y)

(* The step's result and its local error estimate. Without history the
   estimate is half the gap between backward Euler and its explicit Euler
   predictor, which is the leading term of backward Euler's local error. *)
let attempt rhs { at; prev; _ } h =
  let* be = Bdf1.step rhs h at in
  match prev with
  | None -> Ok (be, Vec.scale 0.5 (Vec.sub be (Vec.axpy h (rhs at.t at.y) at.y)))
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
  let span = p.t_end -. p.t0 in
  let dt_max = Option.value dt_max ~default:(span /. 10.) in
  let dt0 = Float.min dt_max (Option.value dt0 ~default:(1e-6 *. span)) in
  (* A step of length [h] reached [next]. *)
  let accept s h next =
    let streak = s.streak + 1 in
    {
      at = next;
      prev = Some { dt_prev = h; y_prev = s.at.y };
      dt = (if streak = 3 then Float.min (2. *. s.dt) dt_max else s.dt);
      streak = streak mod 3;
      failures = 0;
      accepted_steps = s.accepted_steps + 1;
      rejected_steps = s.rejected_steps;
    }
  in
  (* Halve the step that failed, which is shorter than [s.dt] when it was cut
     to land on [t_end], or give up when halving cannot help. *)
  let reject s h =
    let failures = s.failures + 1 and dt = h /. 2. in
    if failures > max_rejects || dt < dt_min s.at.t then Error (Fail.StepRejected failures)
    else Ok { s with dt; streak = 0; failures; rejected_steps = s.rejected_steps + 1 }
  in
  let rec go s =
    if s.at.t >= p.t_end then Ok (solution s)
    else
      (* The last step is cut to land exactly on [t_end]. *)
      let last = s.dt >= p.t_end -. s.at.t in
      let h = if last then p.t_end -. s.at.t else s.dt in
      match attempt p.rhs s h with
      | Ok (y, est) when scaled est y <= tol -> go (accept s h { t = (if last then p.t_end else s.at.t +. h); y })
      | Ok _ | Error _ -> Result.bind (reject s h) go
  in
  if not (span >= 0.) then invalid_arg "Adaptive.integrate: t_end is before t0"
  else if span > 0. && not (dt0 > 0.) then invalid_arg "Adaptive.integrate: dt0 and dt_max must be positive"
  else if not (Vec.finite p.y0 && Vec.finite (p.rhs p.t0 p.y0)) then Error Fail.Nan
  else
    go
      {
        at = { t = p.t0; y = p.y0 };
        prev = None;
        dt = dt0;
        streak = 0;
        failures = 0;
        accepted_steps = 0;
        rejected_steps = 0;
      }
