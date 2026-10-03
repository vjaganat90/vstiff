open Vstiff

(* The corpus is a table of cases: a name and a thunk computing the text
   printed after it. dune diffs this executable's output with corpus.expected. *)

let pp_vec v = "[" ^ String.concat "; " (Array.to_list (Array.map (Printf.sprintf "%.12f") v)) ^ "]"
let failure e = "Error " ^ Fail.to_string e
let show = function Ok x -> "Ok " ^ pp_vec x | Error e -> failure e
let max_abs m = Array.fold_left (fun acc row -> Float.max acc (Vec.norm_inf row)) 0. m
let max_error exact y = Vec.norm_inf (Vec.sub y exact)

(* The adaptive integrator as the corpus runs it: BDF2 with its backward
   Euler error estimate, under the halve/double step-size policy. *)
let bdf2_halving = Adaptive.integrate (module Bdf2) (module Halving)

let error_below bound = function
  | Ok e -> Printf.sprintf "max error %.2e < %g: %b" e bound (e < bound)
  | Error e -> failure e

(* Step 1: damped Newton on a user Jacobian. *)
let newton =
  let a = [| [| 3.; 1. |]; [| 1.; 2. |] |] and b = [| 9.; 8. |] in
  let linear x = Vec.sub (Array.map (fun row -> Vec.dot row x) a) b in
  let square x = [| (x.(0) *. x.(0)) -. 2. |] and square' x = [| [| 2. *. x.(0) |] |] in
  [
    ("newton linear 2d", fun () -> show (Newton.solve linear (fun _ -> a) [| 0.; 0. |]));
    ("newton quadratic x0=1", fun () -> show (Newton.solve square square' [| 1. |]));
    ("newton quadratic x0=0", fun () -> show (Newton.solve square square' [| 0. |]));
  ]

(* Step 2: forward-difference Jacobian against the analytic one of the canary.
   At the origin f(y + h e_j) has no cancellation, so the absolute bound holds
   exactly. At y0 the rounding floor eps |λ_3 y_3| / h is about 3e-5 in
   absolute terms, so there the bound is per entry relative to |J_ij|. *)
let jacobian =
  let open Problems.Canary in
  let analytic = Array.mapi (fun i l -> Array.mapi (fun j _ -> if i = j then -.l else 0.) lambda) lambda in
  let entry_errors y = Array.map2 Vec.sub (Jac.forward (rhs 0.) y) analytic in
  let relative = Array.map2 (Array.map2 (fun e a -> e /. (1. +. Float.abs a))) in
  [
    ( "jac canary at origin, max entry error < 1e-6",
      fun () -> string_of_bool (max_abs (entry_errors (Array.map (fun _ -> 0.) y0)) < 1e-6) );
    ( "jac canary at y0, max entry error relative to |J_ij| < 1e-6",
      fun () -> string_of_bool (max_abs (relative (entry_errors y0) analytic) < 1e-6) );
  ]

(* Step 3: backward Euler, fixed dt, on the canary. *)
let backward_euler =
  let open Problems.Canary in
  [
    ( "bdf1 canary t=1 dt=2e-6",
      fun () -> error_below 1e-6 (Result.map (max_error (exact 1.)) (Stepper.fixed (module Bdf1) ~dt:2e-6 problem)) );
  ]

(* Step 4: BDF2 with a BDF1 startup is second order on the logistic equation. *)
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

(* Step 5: step rejection on van der Pol. *)
let van_der_pol =
  [
    ( "vdp mu=1000 [0,2000] tol=1e-4",
      fun () ->
        match bdf2_halving ~tol:1e-4 Problems.VanDerPol.problem with
        | Ok s ->
            Printf.sprintf "Ok at t=%g, rejected >= 1: %b, y finite: %b" s.t (s.stats.rejected_steps >= 1) (Vec.finite s.y)
        | Error e -> failure e );
  ]

(* Step 6: Robertson to t = 1e4 against the reference, with mass conserved. *)
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

(* Regression pins. *)

(* The canary at a stiff step: h lambda = 10 for the fast rate, so Newton only
   converges with the right Jacobian, and a swapped row makes it diverge. *)
let stiff_canary =
  let open Problems.Canary in
  [
    ( "bdf2 canary t=1 dt=1e-3 (h lambda = 10)",
      fun () -> error_below 1e-6 (Result.map (max_error (exact 1.)) (Stepper.fixed (module Bdf2) ~dt:1e-3 problem)) );
  ]

(* Jacobian orientation: J is not symmetric here, so a transposed J fails. *)
let orientation =
  let f y = [| y.(1); -.y.(0) -. (3. *. y.(1)); (5. *. y.(0)) +. (y.(2) *. y.(2)) |] in
  let y = [| 1.; 2.; 3. |] in
  let analytic = [| [| 0.; 1.; 0. |]; [| -1.; -3.; 0. |]; [| 5.; 0.; 2. *. y.(2) |] |] in
  let relative = Array.map2 (Array.map2 (fun fd a -> (fd -. a) /. (1. +. Float.abs a))) in
  [
    ( "jac non-symmetric at (1, 2, 3), max entry error relative to |J_ij| < 1e-6",
      fun () -> string_of_bool (max_abs (relative (Jac.forward f y) analytic) < 1e-6) );
  ]

(* Partial pivoting: a zero leading entry that needs a row exchange, and a tiny
   one that elimination without the exchange turns into a wrong answer. *)
let pivoting =
  let solved = function Some x -> pp_vec x | None -> "None" in
  [
    ( "linalg zero leading pivot",
      fun () -> solved (Linalg.solve [| [| 0.; 1.; 1. |]; [| 2.; 1.; 0. |]; [| 1.; 0.; 3. |] |] [| 5.; 4.; 10. |]) );
    ("linalg tiny leading pivot", fun () -> solved (Linalg.solve [| [| 1e-20; 1. |]; [| 1.; 1. |] |] [| 1.; 2. |]));
  ]

(* Termination and argument checks. Each right-hand side counts its calls and
   aborts past a budget, so a loop that never ends prints a line instead of
   hanging the run. *)
exception Over_budget

let budgeted rhs =
  let calls = ref 0 in
  fun t y ->
    incr calls;
    if !calls > 5_000_000 then raise Over_budget else rhs t y

let outcome run () =
  match run () with
  | Ok (s : _ Adaptive.solution) -> Printf.sprintf "Ok at t=%g" s.t
  | Error e -> failure e
  | exception Over_budget -> "no answer within 5e6 rhs calls"
  | exception Invalid_argument m -> "Invalid_argument " ^ m

let give_up =
  (* y' = -y from y = 1, on a fresh budget each time. *)
  let decay ~t0 ~t_end = { Ode.rhs = budgeted (fun _ y -> Vec.scale (-1.) y); t0; t_end; y0 = [| 1. |] } in
  let blow_up = { Ode.rhs = budgeted (fun _ y -> [| y.(0) *. y.(0) |]); t0 = 0.; t_end = 2.; y0 = [| 1. |] } in
  let nan_wall () =
    {
      Ode.rhs = budgeted (fun t y -> if t > 0.5 then [| Float.nan |] else Vec.scale (-1.) y);
      t0 = 0.;
      t_end = 1.;
      y0 = [| 1. |];
    }
  in
  [
    ("adaptive blow-up y' = y^2 from y(0)=1 to t=2", outcome (fun () -> bdf2_halving ~tol:1e-6 blow_up));
    ( "adaptive rhs NaN past t=0.5",
      outcome (fun () -> bdf2_halving ~dt0:0.25 ~dt_max:1e6 ~tol:1e-3 (nan_wall ())) );
    ( "adaptive dt0 = dt_max = 1e30 on [0, 1]",
      outcome (fun () -> bdf2_halving ~dt0:1e30 ~dt_max:1e30 ~tol:1e-6 (decay ~t0:0. ~t_end:1.)) );
    ("adaptive empty span", outcome (fun () -> bdf2_halving ~tol:1e-6 (decay ~t0:1. ~t_end:1.)));
    ("adaptive t_end < t0", outcome (fun () -> bdf2_halving ~tol:1e-6 (decay ~t0:1. ~t_end:0.)));
    ("adaptive dt0 = 0", outcome (fun () -> bdf2_halving ~dt0:0. ~tol:1e-6 (decay ~t0:0. ~t_end:1.)));
    ( "bdf1 dt = 0",
      fun () ->
        match Stepper.fixed (module Bdf1) ~dt:0. (decay ~t0:0. ~t_end:1.) with
        | Ok y -> "Ok " ^ pp_vec y
        | Error e -> failure e
        | exception Invalid_argument m -> "Invalid_argument " ^ m );
  ]

(* Step control: from dt0 = 0.5, which is also the default dt_max here, the
   startup estimate rejects the first attempts; after that the halve and double
   rules decide every step, so these counts move if any of them do. With
   dt_max = 1e-3, accuracy alone would allow longer steps, so the cap sets
   every step once dt has grown to it. *)
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

(* Robertson accuracy at tol = 1e-6, far tighter than the 1e-3 acceptance bound. *)
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
    ]

let () = List.iter (fun (name, run) -> Printf.printf "%s: %s\n" name (run ())) corpus
