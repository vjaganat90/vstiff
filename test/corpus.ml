(* The corpus: a table of cases, each a name and a thunk (a function of ()) returning the text Report prints
   after it; dune diffs the output with corpus.expected. Steps 1 to 6 climb from Newton to the adaptive driver,
   so the first changed line names the lowest layer that broke; regression pins follow. See
   docs/numerics/06-the-corpus.md and docs/glossary.md. *)
open Vstiff

(* Twelve decimals hide last-bit noise, yet show a loose Newton tolerance (1e-6 instead of 1e-10) in the root
   of x^2 - 2. *)
let pp_vec v = "[" ^ String.concat "; " (Array.to_list (Array.map (Printf.sprintf "%.12f") v)) ^ "]"
let failure e = "Error " ^ Fail.to_string e

(* Numerical failures are Error values, not exceptions (docs/ocaml.md). *)
let show = function Ok x -> "Ok " ^ pp_vec x | Error e -> failure e
let max_abs m = Array.fold_left (fun acc row -> Float.max acc (Vec.norm_inf row)) 0. m
let max_error exact y = Vec.norm_inf (Vec.sub y exact)

(* BDF2 with its backward Euler error estimate under the halve/double policy; modules are passed as arguments
   (docs/ocaml.md). *)
let bdf2_halving = Adaptive.integrate (module Bdf2) (module Halving)

(* The error prints to three digits, so a line also moves when the error changes but still meets the bound. *)
let error_below bound = function
  | Ok e -> Printf.sprintf "max error %.2e < %g: %b" e bound (e < bound)
  | Error e -> failure e

(* Step 1: Newton alone, with hand-written Jacobians (the integrators build theirs with Jac.forward). A linear
   system takes one full step; x^2 - 2 from 1 converges to sqrt 2; from 0 the 1x1 Jacobian is singular, so
   Newton returns Error Diverged. docs/numerics/02-newton.md *)
let newton =
  let a = [| [| 3.; 1. |]; [| 1.; 2. |] |] and b = [| 9.; 8. |] in
  let linear x = Vec.sub (Array.map (fun row -> Vec.dot row x) a) b in
  let square x = [| (x.(0) *. x.(0)) -. 2. |] and square' x = [| [| 2. *. x.(0) |] |] in
  [
    ("newton linear 2d", fun () -> show (Newton.solve linear (fun _ -> a) [| 0.; 0. |]));
    ("newton quadratic x0=1", fun () -> show (Newton.solve square square' [| 1. |]));
    ("newton quadratic x0=0", fun () -> show (Newton.solve square square' [| 0. |]));
  ]

(* Step 2: Jac.forward against the canary's exact Jacobian, -Lambda. At the origin nothing cancels, so an absolute
   1e-6 holds; at y0 the 1e4 entry carries rounding noise, so errors are taken relative to 1 + |J_ij| (which the
   label abbreviates). docs/numerics/03-jacobians-and-floating-point.md, sections 5 and 6 *)
let jacobian =
  (* let open M in: the names of M are in scope below (docs/ocaml.md). *)
  let open Problems.Canary in
  let analytic = Array.mapi (fun i l -> Array.mapi (fun j _ -> if i = j then -.l else 0.) lambda) lambda in
  let entry_errors y = Array.map2 Vec.sub (Jac.forward (rhs 0.) y) analytic in
  (* Partial application: relative waits for two matrices (docs/ocaml.md). *)
  let relative = Array.map2 (Array.map2 (fun e a -> e /. (1. +. Float.abs a))) in
  [
    ( "jac canary at origin, max entry error < 1e-6",
      fun () -> string_of_bool (max_abs (entry_errors (Array.map (fun _ -> 0.) y0)) < 1e-6) );
    ( "jac canary at y0, max entry error relative to |J_ij| < 1e-6",
      fun () -> string_of_bool (max_abs (relative (entry_errors y0) analytic) < 1e-6) );
  ]

(* Step 3: backward Euler through the whole stack, Stepper down to Linalg. First order: the slow component's
   error at t = 1 is about 0.184 h, so dt = 2e-6 gives the 3.68e-07 of corpus.expected, under the bound 1e-6, in
   500,000 steps that dominate the test time. docs/numerics/01-odes-and-stiffness.md, section 6 *)
let backward_euler =
  let open Problems.Canary in
  [
    ( "bdf1 canary t=1 dt=2e-6",
      fun () -> error_below 1e-6 (Result.map (max_error (exact 1.)) (Stepper.fixed (module Bdf1) ~dt:2e-6 problem)) );
  ]

(* Step 4: BDF2 is second order, so halving dt divides the error at t = 5 by about 4. The window [3.5, 4.5]
   rejects first order (ratio about 2) and third (about 8) and leaves room for higher-order terms. The backward
   Euler first step adds only O(h^2) to the final error; omega = h / h_prev is 1. docs/numerics/04-bdf.md *)
let bdf2_order =
  let open Problems.Logistic in
  let error dt = Result.map (max_error (exact 5.)) (Stepper.fixed (module Bdf2) ~dt problem) in
  [
    ( "bdf2 logistic [0,5]",
      fun () ->
        match (error 0.01, error 0.005) with
        | Ok coarse, Ok fine ->
            let ratio = coarse /. fine in
            Printf.sprintf "error dt=0.01 %.3e, dt=0.005 %.3e, ratio %.2f in [3.5, 4.5]: %b" coarse fine ratio
              (3.5 <= ratio && ratio <= 4.5)
        | Error e, _ | _, Error e -> failure e );
  ]

(* Step 5: step control on a relaxation oscillation. The run must reach t_end, reject at least one step (so the
   rejection path ran) and end finite; there is no closed form or reference, so no accuracy check. s.t and
   s.stats need no module prefix because the type of s is known (docs/ocaml.md). docs/numerics/05-step-control.md *)
let van_der_pol =
  [
    ( "vdp mu=1000 [0,2000] tol=1e-4",
      fun () ->
        match bdf2_halving ~tol:1e-4 Problems.VanDerPol.problem with
        | Ok s ->
            Printf.sprintf "Ok at t=%g, rejected >= 1: %b, y finite: %b" s.t (s.stats.rejected_steps >= 1) (Vec.finite s.y)
        | Error e -> failure e );
  ]

(* Step 6: Robertson against Refs, a value computed outside vstiff. The 1e-3 bound on y1 is loose (about 1% of
   it); the accuracy pin below is tighter. y1 + y2 + y3 must stay 1 within 1e-8: the right-hand sides sum to zero
   and BDF weights sum to one, so only round-off moves it, and 1e-8 leaves room for it to build up. *)
let robertson =
  [
    ( "robertson t=1e4 tol=1e-6",
      fun () ->
        match bdf2_halving ~tol:1e-6 Problems.Robertson.problem with
        | Ok s ->
            let y1_error = Float.abs (s.y.(0) -. Refs.robertson_y1_at_1e4) in
            let mass_error = Float.abs (Array.fold_left ( +. ) 0. s.y -. 1.) in
            Printf.sprintf "|y1 - ref| < 1e-3: %b, |y1 + y2 + y3 - 1| < 1e-8: %b" (y1_error < 1e-3) (mass_error < 1e-8)
        | Error e -> failure e );
  ]

(* Regression pins: each group fixes rules the steps above leave free; its comment says what moves its lines. *)

(* Stiff canary: dt = 1e-3 is h lambda = 10 on the fast rate, five times past the stability limit of explicit
   Euler (h lambda <= 2), yet the error stays under 1e-6. With Jacobian rows 0 and 2 swapped Newton diverges here,
   while step 3 stays true with other digits: the Jacobian sets Newton's speed, not its root. *)
let stiff_canary =
  let open Problems.Canary in
  [
    ( "bdf2 canary t=1 dt=1e-3 (h lambda = 10)",
      fun () -> error_below 1e-6 (Result.map (max_error (exact 1.)) (Stepper.fixed (module Bdf2) ~dt:1e-3 problem)) );
  ]

(* Orientation: this Jacobian is not symmetric, so Jac.forward returning the transpose fails the line (bound 1e-6
   as in step 2); the diagonal canary cannot show that. The y3^2 term is curved, so a coarse step shows too. The
   transpose also makes the Robertson run so slow that a full run needs a timeout to reach this line. *)
let orientation =
  let f y = [| y.(1); -.y.(0) -. (3. *. y.(1)); (5. *. y.(0)) +. (y.(2) *. y.(2)) |] in
  let y = [| 1.; 2.; 3. |] in
  let analytic = [| [| 0.; 1.; 0. |]; [| -1.; -3.; 0. |]; [| 5.; 0.; 2. *. y.(2) |] |] in
  let relative = Array.map2 (Array.map2 (fun fd a -> (fd -. a) /. (1. +. Float.abs a))) in
  [
    ( "jac non-symmetric at (1, 2, 3), max entry error relative to |J_ij| < 1e-6",
      fun () -> string_of_bool (max_abs (relative (Jac.forward f y) analytic) < 1e-6) );
  ]

(* Partial pivoting: a zero leading entry that only a row exchange gets past (None without it), and a tiny one
   that, without the exchange, makes the first unknown 0 instead of 1. *)
let pivoting =
  let solved = function Some x -> pp_vec x | None -> "None" in
  [
    ( "linalg zero leading pivot",
      fun () -> solved (Linalg.solve [| [| 0.; 1.; 1. |]; [| 2.; 1.; 0. |]; [| 1.; 0.; 3. |] |] [| 5.; 4.; 10. |]) );
    ("linalg tiny leading pivot", fun () -> solved (Linalg.solve [| [| 1e-20; 1. |]; [| 1.; 1. |] |] [| 1.; 2. |]));
  ]

(* Termination and argument checks. Each right-hand side runs on Guard's call budget, so a loop that never ends
   prints a line instead of hanging the run; an invalid argument prints its message. *)
let give_up =
  (* The annotation tells OCaml that s.t is a field of Adaptive.solution. *)
  let reached (s : _ Adaptive.solution) = Printf.sprintf "Ok at t=%g" s.t in
  (* y' = -y from y = 1. *)
  let decay ~t0 ~t_end = { Ode.rhs = Guard.budget (fun _ y -> Vec.scale (-1.) y); t0; t_end; y0 = [| 1. |] } in
  let blow_up () = { Ode.rhs = Guard.budget (fun _ y -> [| y.(0) *. y.(0) |]); t0 = 0.; t_end = 2.; y0 = [| 1. |] } in
  let nan_wall () =
    {
      Ode.rhs = Guard.budget (fun t y -> if t > 0.5 then [| Float.nan |] else Vec.scale (-1.) y);
      t0 = 0.;
      t_end = 1.;
      y0 = [| 1. |];
    }
  in
  [
    (* No way forward: y' = y^2 blows up at t = 1 (exact solution 1 / (1 - t)), and the NaN wall has no values
       past t = 0.5. The step floor 16 eps |t| in Halving ends both early; without it only max_rejects ends
       them, with StepRejected 51. Why the wall takes 46 rejections: docs/numerics/05-step-control.md. *)
    ("adaptive blow-up y' = y^2 from y(0)=1 to t=2", Guard.run reached (fun () -> bdf2_halving ~tol:1e-6 (blow_up ())));
    ( "adaptive rhs NaN past t=0.5",
      Guard.run reached (fun () -> bdf2_halving ~dt0:0.25 ~dt_max:1e6 ~tol:1e-3 (nan_wall ())) );
    (* The first step must be cut to the span, 1, and a rejection must halve that cut step. Without the cut, or
       halving the 1e30 proposal instead, max_rejects runs out first: StepRejected 51. *)
    ( "adaptive dt0 = dt_max = 1e30 on [0, 1]",
      Guard.run reached (fun () -> bdf2_halving ~dt0:1e30 ~dt_max:1e30 ~tol:1e-6 (decay ~t0:0. ~t_end:1.)) );
    (* A zero span is valid: nothing to integrate. *)
    ("adaptive empty span", Guard.run reached (fun () -> bdf2_halving ~tol:1e-6 (decay ~t0:1. ~t_end:1.)));
    (* The last three lines are Check's messages; without the checks they print Ok or Error StepRejected instead. *)
    ("adaptive t_end < t0", Guard.run reached (fun () -> bdf2_halving ~tol:1e-6 (decay ~t0:1. ~t_end:0.)));
    ("adaptive dt0 = 0", Guard.run reached (fun () -> bdf2_halving ~dt0:0. ~tol:1e-6 (decay ~t0:0. ~t_end:1.)));
    ( "bdf1 dt = 0",
      Guard.run (fun y -> "Ok " ^ pp_vec y) (fun () -> Stepper.fixed (module Bdf1) ~dt:0. (decay ~t0:0. ~t_end:1.)) );
  ]

(* Step-size counts. From dt0 = 0.5, also the default dt_max here, the startup estimate rejects the first attempts,
   then the halve and double rules drive the counts; with dt_max = 1e-3 the cap, not accuracy, sets the step.
   Derivation, and what each rule moves: docs/numerics/06-the-corpus.md, sections 5 and 9. *)
let step_control =
  let open Problems.Logistic in
  let counts (s : Halving.stats Adaptive.solution) =
    Printf.sprintf "accepted %d, rejected %d" s.stats.accepted_steps s.stats.rejected_steps
  in
  [
    ( "adaptive logistic [0,5] dt0=0.5 tol=1e-6",
      fun () ->
        match bdf2_halving ~dt0:0.5 ~tol:1e-6 problem with
        | Ok s -> Printf.sprintf "%s, max error %.2e" (counts s) (max_error (exact 5.) s.y)
        | Error e -> failure e );
    ( "adaptive logistic [0,5] dt_max=1e-3 tol=1e-6",
      fun () -> match bdf2_halving ~dt_max:1e-3 ~tol:1e-6 problem with Ok s -> counts s | Error e -> failure e );
  ]

(* Robertson accuracy: the error is 7.2e-07, and 1e-5 sits far below the 1e-3 of step 6. Returning backward Euler's
   result instead of BDF2's (error 1.6e-04), or not cutting the last step to land on t_end, fails this line and
   not step 6. tol limits each step's estimated error, not the final error. *)
let robertson_accuracy =
  [
    ( "robertson t=1e4 tol=1e-6 accuracy",
      fun () ->
        match bdf2_halving ~tol:1e-6 Problems.Robertson.problem with
        | Ok s ->
            let e = Float.abs (s.y.(0) -. Refs.robertson_y1_at_1e4) in
            Printf.sprintf "|y1 - ref| = %.1e < 1e-5: %b" e (e < 1e-5)
        | Error e -> failure e );
  ]

(* Newton must not call an overflow a root: here a too-small Jacobian makes the
   converged step run past max_float. *)
let newton_overflow =
  [
    ( "newton step that converges into an overflow",
      Guard.run
        (fun x -> "Ok " ^ pp_vec x)
        (fun () ->
          Newton.solve (fun x -> [| x.(0) -. Float.max_float |]) (fun _ -> [| [| 0.25 |] |]) [| Float.max_float -. 1e297 |])
    );
  ]

(* The clock: a step below half an ulp of t cannot move t, and a step of a few
   ulps moves t by less than its nominal length. Steps are snapped to the
   floats, a step that cannot move t is rejected, and a remainder below t's
   resolution is taken as the last step. Each run is on a call budget. *)
let clock =
  let y_is_t ~t0 ~t_end = { Ode.rhs = Guard.budget (fun _ _ -> [| 1. |]); t0; t_end; y0 = [| 0. |] } in
  let landed t_end (s : _ Adaptive.solution) = Printf.sprintf "Ok, at t_end: %b, y = %.6g" (s.t = t_end) s.y.(0) in
  let one_ulp = Float.succ 1. in
  [
    ( "adaptive y' = 1 over one ulp from t=1",
      Guard.run (landed one_ulp) (fun () -> bdf2_halving ~tol:1e-6 (y_is_t ~t0:1. ~t_end:one_ulp)) );
    ( "adaptive y' = 1 from t=1e10, dt_max = 5e-7 is below t's resolution",
      Guard.run (landed (1e10 +. 1.)) (fun () ->
          bdf2_halving ~dt_max:5e-7 ~tol:1e-6 (y_is_t ~t0:1e10 ~t_end:(1e10 +. 1.))) );
    ( "adaptive y' = 1 from t=1e15 over 100, dt0 = dt_max = 0.19",
      Guard.run (landed (1e15 +. 100.)) (fun () ->
          bdf2_halving ~dt0:0.19 ~dt_max:0.19 ~tol:1e-6 (y_is_t ~t0:1e15 ~t_end:(1e15 +. 100.))) );
  ]

(* A step that cannot move t is rejected before the method sees it: the only
   call of the right-hand side is the driver's check of the initial state. *)
let too_small =
  [
    ( "adaptive step that cannot move t never reaches the method",
      fun () ->
        let rhs, calls = Instrument.count (Guard.budget (fun _ _ -> [| 1. |])) in
        let outcome =
          Guard.run
            (fun _ -> "Ok")
            (fun () -> bdf2_halving ~dt_max:5e-7 ~tol:1e-6 { Ode.rhs; t0 = 1e10; t_end = 1e10 +. 1.; y0 = [| 0. |] })
            ()
        in
        Printf.sprintf "%s, rhs calls: %d" outcome (calls ()) );
  ]

(* The fixed-step clock: step k ends at t0 + k h and the last one at t_end
   itself, so a right-hand side undefined past t_end is never evaluated there.
   The stage is evaluated at the new time, which only a time-dependent problem
   can see: on y' = 2t (y(1) = 1) BDF2 is second order only if it is. *)
let fixed_clock =
  let sqrt_rest = { Ode.rhs = (fun t _ -> [| Float.sqrt (1. -. t) |]); t0 = 0.; t_end = 1.; y0 = [| 0. |] } in
  let ramp = { Ode.rhs = (fun t _ -> [| 2. *. t |]); t0 = 0.; t_end = 1.; y0 = [| 0. |] } in
  let decay_at_1e10 = { Ode.rhs = (fun _ y -> Vec.scale (-1.) y); t0 = 1e10; t_end = 1e10 +. 1.; y0 = [| 1. |] } in
  [
    ( "bdf1 y' = sqrt(1 - t) on [0, 1], dt = 1/9",
      Guard.run (fun y -> "Ok " ^ pp_vec y) (fun () -> Stepper.fixed (module Bdf1) ~dt:(1. /. 9.) sqrt_rest) );
    ( "bdf2 y' = 2t on [0, 1]",
      fun () ->
        let error dt = Result.map (max_error [| 1. |]) (Stepper.fixed (module Bdf2) ~dt ramp) in
        match (error 0.1, error 0.05) with
        | Ok coarse, Ok fine ->
            let ratio = coarse /. fine in
            Printf.sprintf "error dt=0.1 %.3e, dt=0.05 %.3e, ratio %.2f in [3.5, 4.5]: %b" coarse fine ratio
              (3.5 <= ratio && ratio <= 4.5)
        | Error e, _ | _, Error e -> failure e );
    ( "bdf1 dt = 1e-7 at t=1e10, below t's resolution",
      Guard.run (fun y -> "Ok " ^ pp_vec y) (fun () -> Stepper.fixed (module Bdf1) ~dt:1e-7 decay_at_1e10) );
  ]

(* Table order is output order, the order of the lines in corpus.expected. *)
let corpus =
  List.concat
    [
      newton;
      jacobian;
      backward_euler;
      bdf2_order;
      van_der_pol;
      robertson;
      stiff_canary;
      orientation;
      pivoting;
      give_up;
      step_control;
      robertson_accuracy;
      newton_overflow;
      clock;
      too_small;
      fixed_clock;
    ]

let () = Report.lines corpus
