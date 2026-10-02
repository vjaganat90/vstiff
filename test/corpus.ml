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
