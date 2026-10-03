open Vstiff

(* y6 y8 is the one nonlinear term, shared by three components; bound by let, it is rounded the same way on every
   platform and in Python. The sums keep the order of the Python twin, which transcription.py compares. *)
let rhs _t y =
  let y1 = y.(0) and y2 = y.(1) and y3 = y.(2) and y4 = y.(3) in
  let y5 = y.(4) and y6 = y.(5) and y7 = y.(6) and y8 = y.(7) in
  let r = 280. *. y6 *. y8 in
  [|
    (-1.71 *. y1) +. (0.43 *. y2) +. (8.32 *. y3) +. 0.0007;
    (1.71 *. y1) -. (8.75 *. y2);
    (-10.03 *. y3) +. (0.43 *. y4) +. (0.035 *. y5);
    (8.32 *. y2) +. (1.71 *. y3) -. (1.12 *. y4);
    (-1.745 *. y5) +. (0.43 *. y6) +. (0.43 *. y7);
    -.r +. (0.69 *. y4) +. (1.71 *. y5) -. (0.43 *. y6) +. (0.69 *. y7);
    r -. (1.81 *. y7);
    -.r +. (1.81 *. y7);
  |]

let y0 = [| 1.; 0.; 0.; 0.; 0.; 0.; 0.; 0.0057 |]
let problem = { Ode.rhs; t0 = 0.; t_end = 321.8122; y0 }
