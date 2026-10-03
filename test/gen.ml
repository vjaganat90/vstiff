(* Integrated shrinking, as in Hedgehog: a value carries its smaller candidates, so map and bind shrink what they
   build. The trees are lazy (Seq), and a run forces only the path that a failing case shrinks along. *)

type 'a tree = Node of 'a * 'a tree Seq.t
type 'a t = Random.State.t -> 'a tree

let root (Node (x, _)) = x
let rec map_tree f (Node (x, smaller)) = Node (f x, Seq.map (map_tree f) smaller)
let return x _ = Node (x, Seq.empty)
let map f g st = map_tree f (g st)

(* f draws from a copy of one split-off state every time, so a shrunk a draws again deterministically. *)
let bind g f st =
  let inner = Random.State.split st in
  let rec go (Node (a, smaller)) =
    let (Node (b, more)) = f a (Random.State.copy inner) in
    Node (b, Seq.append (Seq.map go smaller) more)
  in
  go (g st)

module Syntax = struct
  let ( let* ) = bind
  let ( let+ ) g f = map f g
end

let rec pair_tree (Node (a, xs) as ta) (Node (b, ys) as tb) =
  Node ((a, b), Seq.append (Seq.map (fun ta' -> pair_tree ta' tb) xs) (Seq.map (pair_tree ta) ys))

(* The order of evaluation of function arguments is unspecified: draw a, then b. *)
let pair ga gb st =
  let ta = ga st in
  pair_tree ta (gb st)

let triple ga gb gc = map (fun (a, (b, c)) -> (a, b, c)) (pair ga (pair gb gc))

(* lo first, then halfway back up, then a quarter, ..., x - 1: a binary search for the smallest failing value. *)
let rec int_tree lo x =
  Node (x, Seq.map (int_tree lo) (Seq.unfold (fun d -> if d = 0 then None else Some (x - d, d / 2)) (x - lo)))

let int ~lo ~hi st = int_tree lo (Random.State.int_in_range st ~min:lo ~max:hi)

(* d cut to its first k significant bits, toward 0. *)
let cut k d =
  let m, e = Float.frexp d in
  Float.ldexp (Float.trunc (Float.ldexp m k)) (e - k)

(* o, then o + (x - o) cut to 1, 2, ... bits: a shrunk value has a short binary expansion, and halving toward o
   instead would sink into the subnormals whenever every nonzero value fails. *)
let rec float_tree lo hi o x =
  let candidate k = Float.min hi (Float.max lo (o +. cut k (x -. o))) in
  (* Consecutive cuts are equal while the bits they add are 0. *)
  let fresh k = if candidate k = x || candidate k = candidate (k - 1) then None else Some (candidate k) in
  let smaller = if x = o then Seq.empty else Seq.cons o (Seq.filter_map fresh (Seq.init 53 succ)) in
  Node (x, Seq.map (float_tree lo hi o) smaller)

let float ~lo ~hi st =
  float_tree lo hi (Float.min hi (Float.max lo 0.)) (Float.min hi (lo +. Random.State.float st (hi -. lo)))

(* In base 2, so that the exponents shrunk to integers give powers of two. *)
let log_uniform ~lo ~hi =
  map (fun e -> Float.min hi (Float.max lo (Float.pow 2. e))) (float ~lo:(Float.log2 lo) ~hi:(Float.log2 hi))

let specials =
  [| 0.; -0.; Float.succ 0.; Float.infinity; Float.neg_infinity; Float.nan; Float.max_float; Float.epsilon |]

let special st =
  let x = specials.(Random.State.int st (Array.length specials)) in
  Node (x, if x = 0. then Seq.empty else Seq.return (Node (0., Seq.empty)))

let remove i a = Array.init (Array.length a - 1) (fun k -> if k < i then a.(k) else a.(k + 1))
let replace i t a = Array.mapi (fun k s -> if k = i then t else s) a

(* lens are the lengths the length generator shrinks to from here: one less tries removing each element in turn,
   anything shorter keeps a prefix. *)
let rec array_tree lens ts =
  let n = Array.length ts in
  let drop (Node (m, shorter)) =
    if m = n - 1 then Seq.init n (fun i -> array_tree shorter (remove i ts))
    else if 0 <= m && m < n then Seq.return (array_tree shorter (Array.sub ts 0 m))
    else Seq.empty
  in
  let shrink i =
    let (Node (_, smaller)) = ts.(i) in
    Seq.map (fun t -> array_tree lens (replace i t ts)) smaller
  in
  Node (Array.map root ts, Seq.append (Seq.flat_map drop lens) (Seq.flat_map shrink (Seq.init n Fun.id)))

(* Array.init calls its function on 0, 1, ... in order, so the elements are drawn in order. *)
let array len g st =
  let (Node (n, shorter)) = len st in
  array_tree shorter (Array.init n (fun _ -> g st))

let one_of gs = bind (int ~lo:0 ~hi:(List.length gs - 1)) (List.nth gs)

let sum n f =
  let rec go k acc = if k = n then acc else go (k + 1) (acc +. f k) in
  go 0 0.

let transpose a = Array.init (if a = [||] then 0 else Array.length a.(0)) (fun j -> Array.map (fun row -> row.(j)) a)

(* (I - 2 v v^T / v^T v) a, an orthogonal map of the columns; the zero vector stands for I. *)
let reflect a v =
  let n = Array.length v in
  let vv = sum n (fun k -> v.(k) *. v.(k)) in
  if vv = 0. then a
  else
    let w = Array.init n (fun j -> sum n (fun k -> v.(k) *. a.(k).(j)) /. vv) in
    Array.mapi (fun i row -> Array.mapi (fun j x -> x -. (2. *. v.(i) *. w.(j))) row) a

(* U diag(sigma) V^T is the transpose of V (U diag(sigma))^T, reflections being symmetric. *)
let prescribed sigma =
  let n = Array.length sigma in
  let reflections = array (int ~lo:0 ~hi:n) (array (return n) (float ~lo:(-1.) ~hi:1.)) in
  let diagonal = Array.init n (fun i -> Array.init n (fun j -> if i = j then sigma.(i) else 0.)) in
  let reflect_all vs a = Array.fold_left reflect a vs in
  map (fun (u, v) -> transpose (reflect_all v (transpose (reflect_all u diagonal)))) (pair reflections reflections)

(* Row j copies row i: k skips over i, so j <> i. *)
let singular n =
  let rows = array (return n) (array (return n) (float ~lo:(-1.) ~hi:1.)) in
  map
    (fun (a, (i, k)) ->
      let j = if k >= i then k + 1 else k in
      Array.mapi (fun r row -> if r = j then a.(i) else row) a)
    (pair rows (pair (int ~lo:0 ~hi:(n - 1)) (int ~lo:0 ~hi:(n - 2))))
