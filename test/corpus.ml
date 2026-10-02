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
