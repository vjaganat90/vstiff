open Vstiff

let pp_vec v =
  "[" ^ String.concat "; " (Array.to_list (Array.map (Printf.sprintf "%.12f") v)) ^ "]"

let show = function
  | Ok x -> "Ok " ^ pp_vec x
  | Error e -> "Error " ^ Fail.to_string e

let case name result = Printf.printf "%s: %s\n" name result

(* Step 1: damped Newton on a user Jacobian. *)
let () =
  let a = [| [| 3.; 1. |]; [| 1.; 2. |] |] and b = [| 9.; 8. |] in
  let f x = Vec.sub (Array.map (fun row -> Vec.dot row x) a) b in
  case "newton linear 2d" (show (Newton.solve ~f ~jac:(fun _ -> a) [| 0.; 0. |]));
  let f x = [| (x.(0) *. x.(0)) -. 2. |] and jac x = [| [| 2. *. x.(0) |] |] in
  case "newton quadratic x0=1" (show (Newton.solve ~f ~jac [| 1. |]));
  case "newton quadratic x0=0" (show (Newton.solve ~f ~jac [| 0. |]))

let max_abs m = Array.fold_left (fun acc row -> Float.max acc (Vec.norm_inf row)) 0. m

(* Step 2: forward-difference Jacobian against the analytic one of the canary.
   At the origin f(y + h e_j) has no cancellation, so the absolute bound checks
   layout, step and divisor exactly. At y0 the rounding floor eps |λ_3 y_3| / h
   is about 3e-5 in absolute terms, so there the bound is per entry relative
   to |J_ij|. *)
let () =
  let open Problems.Canary in
  let analytic = Array.mapi (fun i l -> Array.mapi (fun j _ -> if i = j then -.l else 0.) lambda) lambda in
  let entry_errors y = Array.map2 Vec.sub (Jac.forward (rhs 0.) y) analytic in
  let relative = Array.map2 (Array.map2 (fun e a -> e /. (1. +. Float.abs a))) in
  let at_origin = max_abs (entry_errors (Array.map (fun _ -> 0.) y0)) in
  let at_y0 = max_abs (relative (entry_errors y0) analytic) in
  case "jac canary at origin, max entry error < 1e-6" (string_of_bool (at_origin < 1e-6));
  case "jac canary at y0, max entry error relative to |J_ij| < 1e-6" (string_of_bool (at_y0 < 1e-6))

let max_error exact y = Vec.norm_inf (Vec.sub y exact)

let show_error name bound = function
  | Ok e -> case name (Printf.sprintf "max error %.2e < %g: %b" e bound (e < bound))
  | Error e -> case name ("Error " ^ Fail.to_string e)

(* Step 3: backward Euler, fixed dt, on the canary. *)
let () =
  let open Problems.Canary in
  Bdf1.integrate ~rhs ~t0:0. ~t_end:1. ~dt:2e-6 y0
  |> Result.map (max_error (exact 1.))
  |> show_error "bdf1 canary t=1 dt=2e-6" 1e-6

(* Step 4: BDF2 with a BDF1 startup is second order on the logistic equation. *)
let () =
  let open Problems.Logistic in
  let error dt = Result.map (max_error (exact 5.)) (Bdf2.integrate ~rhs ~t0:0. ~t_end:5. ~dt y0) in
  match (error 0.01, error 0.005) with
  | Ok coarse, Ok fine ->
      let ratio = coarse /. fine in
      case "bdf2 logistic [0,5]"
        (Printf.sprintf "error dt=0.01 %.3e, dt=0.005 %.3e, ratio %.2f in [3.5, 4.5]: %b" coarse fine
           ratio
           (3.5 <= ratio && ratio <= 4.5))
  | Error e, _ | _, Error e -> case "bdf2 logistic [0,5]" ("Error " ^ Fail.to_string e)

(* Step 5: step rejection on van der Pol. *)
let () =
  let open Problems.Van_der_pol in
  match Adaptive.integrate ~tol:1e-4 ~rhs ~t0:0. ~t_end:2000. y0 with
  | Ok s ->
      case "vdp mu=1000 [0,2000] tol=1e-4"
        (Printf.sprintf "Ok at t=%g, rejected >= 1: %b, y finite: %b" s.t (s.rejected >= 1) (Vec.finite s.y))
  | Error e -> case "vdp mu=1000 [0,2000] tol=1e-4" ("Error " ^ Fail.to_string e)

(* Step 6: Robertson to t = 1e4 against the reference, with mass conserved. *)
let () =
  let open Problems.Robertson in
  match Adaptive.integrate ~tol:1e-6 ~rhs ~t0:0. ~t_end:1e4 y0 with
  | Ok s ->
      let y1_error = Float.abs (s.y.(0) -. Refs.robertson_y1_at_1e4) in
      let mass_error = Float.abs (Array.fold_left ( +. ) 0. s.y -. 1.) in
      case "robertson t=1e4 tol=1e-6"
        (Printf.sprintf "|y1 - ref| < 1e-3: %b, |y1 + y2 + y3 - 1| < 1e-8: %b" (y1_error < 1e-3)
           (mass_error < 1e-8))
  | Error e -> case "robertson t=1e4 tol=1e-6" ("Error " ^ Fail.to_string e)

(* Regression pins. *)

(* The canary at a stiff step: h lambda = 10 for the fast rate, so Newton only
   converges with the right Jacobian, and a swapped row makes it diverge. *)
let () =
  let open Problems.Canary in
  Bdf2.integrate ~rhs ~t0:0. ~t_end:1. ~dt:1e-3 y0
  |> Result.map (max_error (exact 1.))
  |> show_error "bdf2 canary t=1 dt=1e-3 (h lambda = 10)" 1e-6

(* Jacobian orientation: J is not symmetric here, so a transposed J fails. *)
let () =
  let f y = [| y.(1); -.y.(0) -. (3. *. y.(1)); (5. *. y.(0)) +. (y.(2) *. y.(2)) |] in
  let y = [| 1.; 2.; 3. |] in
  let analytic = [| [| 0.; 1.; 0. |]; [| -1.; -3.; 0. |]; [| 5.; 0.; 2. *. y.(2) |] |] in
  let relative = Array.map2 (Array.map2 (fun fd a -> (fd -. a) /. (1. +. Float.abs a))) in
  let err = max_abs (relative (Jac.forward f y) analytic) in
  case "jac non-symmetric at (1, 2, 3), max entry error relative to |J_ij| < 1e-6" (string_of_bool (err < 1e-6))

(* Partial pivoting: a zero leading entry that needs a row exchange, and a tiny
   one that elimination without the exchange turns into a wrong answer. *)
let () =
  let solved = function Some x -> pp_vec x | None -> "None" in
  case "linalg zero leading pivot"
    (solved (Linalg.solve [| [| 0.; 1.; 1. |]; [| 2.; 1.; 0. |]; [| 1.; 0.; 3. |] |] [| 5.; 4.; 10. |]));
  case "linalg tiny leading pivot" (solved (Linalg.solve [| [| 1e-20; 1. |]; [| 1.; 1. |] |] [| 1.; 2. |]))
