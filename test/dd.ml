(* Knuth's two-sum and the exact error of a product from one fma, combined as in the QD library of Hida, Li and
   Bailey. The compiler may fuse a product into a sum on arm64; only the low-order
   terms of mul have such a product, where fusing can only help. *)

type t = { hi : float; lo : float }

let of_float x = { hi = x; lo = 0. }
let to_float { hi; lo } = hi +. lo

(* s + e = a + b exactly. *)
let two_sum a b =
  let s = a +. b in
  let bb = s -. a in
  (s, a -. (s -. bb) +. (b -. bb))

(* The same when |a| >= |b|, which normalizes a result so that hi is its nearest float. *)
let quick_two_sum a b =
  let s = a +. b in
  (s, b -. (s -. a))

let add x y =
  let s, e = two_sum x.hi y.hi in
  let t, f = two_sum x.lo y.lo in
  let s, e = quick_two_sum s (e +. t) in
  let hi, lo = quick_two_sum s (e +. f) in
  { hi; lo }

let sub x y = add x { hi = -.y.hi; lo = -.y.lo }

let mul x y =
  let p = x.hi *. y.hi in
  let e = Float.fma x.hi y.hi (-.p) in
  let hi, lo = quick_two_sum p (e +. ((x.hi *. y.lo) +. (x.lo *. y.hi))) in
  { hi; lo }

let dot a b = Array.fold_left add (of_float 0.) (Array.mapi (fun i ai -> mul (of_float ai) (of_float b.(i))) a)
