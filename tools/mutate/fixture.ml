(* The input of the tool's own test: every mutation operator has a site here. It is compiled and never run. *)

let area w h = w *. h /. 2.
let shift k = ((k + 1) * 2) - (k / 3)
let flip x = x *. -1.
let total = Array.fold_left ( +. ) 0.
let double = Array.map (( *. ) 2.)
let tiny ?(eps = 1e-8) x = if Float.abs x <= eps then 0. else x
let clamp lo hi x = if x < lo then lo else if x > hi then hi else x
let inside lo hi x = (lo <= x && not (x >= hi)) || lo = hi
let differs a b = a <> b

let drain n =
  let r = ref n in
  while !r > 0 do
    decr r
  done;
  !r

let sign x = match compare x 0. with 0 -> 0 | c when c > 0 -> 1 | _ -> 0
let head = function [] -> None | x :: _ -> Some x
let quotient a b = try a / b with Division_by_zero -> 0
let finite x = if not (Float.is_nan x) then x else 0.
let pick = function Some x when not (x < 0) -> x | _ -> 0
let gap a b = a -. b
