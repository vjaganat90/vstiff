let points = 100

(* SplitMix64 as a function of a counter, so nothing mutates and the points are the same on every platform: integer
   arithmetic is exact. The k-th number is the mix of seed + (k + 1) gamma, 53 of its bits as a double in [0, 1). *)
let gamma = 0x9E3779B97F4A7C15L

let mix z =
  let open Int64 in
  let z = mul (logxor z (shift_right_logical z 30)) 0xBF58476D1CE4E5B9L in
  let z = mul (logxor z (shift_right_logical z 27)) 0x94D049BB133111EBL in
  logxor z (shift_right_logical z 31)

let seed = 2026L

let uniform k =
  let z = mix (Int64.add seed (Int64.mul gamma (Int64.of_int (k + 1)))) in
  Int64.to_float (Int64.shift_right_logical z 11) *. 0x1p-53

(* Each (problem, point, component) has a counter of its own; component 127 draws t, and no problem has that many. *)
let draw ~case ~point component = uniform ((((case * points) + point) * 128) + component)

let lines () =
  let line case_index (case : Suite.case) point =
    let between k (lo, hi) = lo +. ((hi -. lo) *. draw ~case:case_index ~point k) in
    let y = Array.mapi between case.box in
    let t = case.problem.t_end *. draw ~case:case_index ~point 127 in
    let f = case.problem.rhs t y in
    let numbers = (t :: Array.to_list y) @ Array.to_list f in
    String.concat "," (case.name :: List.map (Printf.sprintf "%.17g") numbers)
  in
  List.concat (List.mapi (fun i case -> List.init points (line i case)) Suite.cases)
