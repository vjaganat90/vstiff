open Vstiff

let cells = 40
let alpha = 1. /. 50.
let a = 1.
let b = 3.

(* alpha / dx^2 with dx = 1 / (cells + 1), as a product of exact values, so Python computes the same double. *)
let diffusion = alpha *. float_of_int ((cells + 1) * (cells + 1))

(* u_i and v_i are y.(2 i) and y.(2 i + 1); beyond the ends they are the boundary values 1 and 3. A product that is
   added or subtracted is bound by let, so it is rounded before the sum on every platform, as in the Python twin;
   products by 2 and 4 are exact and need no binding. *)
let rhs _t y =
  let u i = if i < 0 || i >= cells then 1. else y.(2 * i) in
  let v i = if i < 0 || i >= cells then 3. else y.((2 * i) + 1) in
  Array.init (2 * cells) (fun k ->
      let i = k / 2 in
      let ui = u i and vi = v i in
      let uuv = ui *. ui *. vi in
      if k mod 2 = 0 then
        let flux = diffusion *. (u (i - 1) -. (2. *. ui) +. u (i + 1)) in
        a +. uuv -. ((b +. 1.) *. ui) +. flux
      else
        let bu = b *. ui in
        let flux = diffusion *. (v (i - 1) -. (2. *. vi) +. v (i + 1)) in
        bu -. uuv +. flux)

let y0 =
  Array.init (2 * cells) (fun k ->
      if k mod 2 = 1 then 3.
      else
        let x = float_of_int ((k / 2) + 1) *. (1. /. float_of_int (cells + 1)) in
        1. +. sin (2. *. Float.pi *. x))

let problem = { Ode.rhs; t0 = 0.; t_end = 10.; y0 }
