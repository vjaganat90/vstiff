(* The property suite: each line runs one property on generated cases and prints ok with the count, or the seed and
   the shrunk case that failed. The kernel properties, then shadows of two theorems about Bdf2.coeffs: T1, the
   formula is exact on quadratics, and T2, the differences of the homogeneous recurrence contract. Each bound is derived in the comment above its property; u is
   the unit roundoff, the largest relative error of one rounding. docs/testing.md, Properties. *)
open Vstiff
open Numerics

let eps = Float.epsilon
let u = eps /. 2.

(* 17 significant digits read back as the same float, so a printed case can be run again. *)
let num = Printf.sprintf "%.17g"
let vector v = "[" ^ String.concat "; " (Array.to_list (Array.map num v)) ^ "]"
let matrix m = "[" ^ String.concat "; " (Array.to_list (Array.map vector m)) ^ "]"

(* An error over its bound, which a property needs at most 1; 0 when both are 0. *)
let ratio error bound = if error = 0. then 0. else error /. bound
let entry = Gen.float ~lo:(-1.) ~hi:1.
let vector_of n g = Gen.array (Gen.return n) g
let signed neg g = Gen.map (fun (s, x) -> if s = 0 then x else neg x) (Gen.pair (Gen.int ~lo:0 ~hi:1) g)
let magnitude ~lo ~hi = signed Float.neg (Gen.log_uniform ~lo ~hi)

(* Singular values from 1 down to 1 / kappa, evenly spaced in the logarithm. *)
let spread n kappa = Array.init n (fun i -> if n = 1 then 1. else kappa ** (-.Float.of_int i /. Float.of_int (n - 1)))

(* A system of size 1 to max_n with 2-norm condition number up to max_kappa, right-hand side in [-1, 1]. *)
let system max_n max_kappa =
  let open Gen.Syntax in
  let* n, kappa = Gen.pair (Gen.int ~lo:1 ~hi:max_n) (Gen.log_uniform ~lo:1. ~hi:max_kappa) in
  Gen.pair (Gen.prescribed (spread n kappa)) (vector_of n entry)

let show_system (a, b) = Printf.sprintf "a = %s, b = %s" (matrix a) (vector b)
let abs_sum v = Array.fold_left (fun s x -> s +. Float.abs x) 0. v
let norm a = Array.fold_left (fun m row -> Float.max m (abs_sum row)) 0. a

(* A x - b in double-double, rounded once: its own rounding is far below what the bounds allow. *)
let residual a x b = Array.mapi (fun i row -> Dd.to_float (Dd.sub (Dd.dot row x) (Dd.of_float b.(i)))) a
let dim a = Float.of_int (Array.length a)

(* Vec: identities that hold exactly, and axpy within its rounding. *)

(* Any finite magnitude, or a special value. *)
let any = Gen.one_of [ magnitude ~lo:1e-300 ~hi:1e300; Gen.special ]

(* Each entry carries a sort key; sorting by the keys permutes the vector, and keys shrunk to 0 leave it in place.
   Float.max keeps a nan wherever it appears, so the norm is the same exactly, nan included (Float.equal). *)
module VecNormPermuted = struct
  type t = (float * int) array

  let name = "vec norm_inf is invariant under permutation"
  let gen = Gen.array (Gen.int ~lo:0 ~hi:8) (Gen.pair any (Gen.int ~lo:0 ~hi:8))

  let show c =
    vector (Array.map fst c)
    ^ ", keys "
    ^ String.concat " " (Array.to_list (Array.map (fun (_, k) -> string_of_int k) c))

  let holds c =
    let sorted = List.stable_sort (fun (_, j) (_, k) -> compare j k) (Array.to_list c) in
    Float.equal (Vec.norm_inf (Array.map fst c)) (Vec.norm_inf (Array.of_list (List.map fst sorted)))
end

(* x - x is 0 exactly when x is finite: infinity - infinity and nan - nan are nan. *)
module VecFinite = struct
  type t = float array

  let name = "vec finite is false exactly when an entry is nan or infinite"
  let gen = Gen.array (Gen.int ~lo:0 ~hi:8) any
  let show = vector
  let holds v = Vec.finite v = Array.for_all (fun x -> x -. x = 0.) v
end

(* Fused or not, a x + y takes at most two roundings, of a x and of the sum, so the result is within
   u (|a x + y| + |a x|) of the exact value, plus 2^-1074 should it be subnormal. The magnitudes keep a x away from
   overflow and underflow, and the exact value comes from double-double. A longer y is cut to x's length. *)
module VecAxpy = struct
  type t = float * (float * float) array * float array

  let name = "vec axpy is a x + y within its rounding"
  let finite = Gen.one_of [ Gen.return 0.; magnitude ~lo:1e-150 ~hi:1e150 ]

  let gen =
    Gen.triple finite
      (Gen.array (Gen.int ~lo:0 ~hi:8) (Gen.pair finite finite))
      (Gen.array (Gen.int ~lo:0 ~hi:2) finite)

  let show (a, xy, extra) =
    Printf.sprintf "a = %s, x = %s, y = %s" (num a)
      (vector (Array.map fst xy))
      (vector (Array.append (Array.map snd xy) extra))

  let error a (x, y) r =
    let exact = Dd.add (Dd.mul (Dd.of_float a) (Dd.of_float x)) (Dd.of_float y) in
    ratio
      (Float.abs (Dd.to_float (Dd.sub (Dd.of_float r) exact)))
      ((u *. (Float.abs (Dd.to_float exact) +. Float.abs (a *. x))) +. Float.succ 0.)

  let worst (a, xy, extra) =
    let r = Vec.axpy a (Array.map fst xy) (Array.append (Array.map snd xy) extra) in
    if Array.length r <> Array.length xy then infinity else Array.fold_left Float.max 0. (Array.map2 (error a) xy r)

  let holds c = worst c <= 1.
end

(* Linalg. Gaussian elimination with partial pivoting is backward stable: (A + dA) x = b with |dA| <= gamma_3n |L| |U|
   (Higham, Accuracy and Stability of Numerical Algorithms, chapter 9), a worst case that grows with n and with the
   growth factor. These matrices have orthogonal factors, which keep the growth factor small, and their roundings do
   not line up: the normwise backward error |A x - b| / (|A| |x| + |b|) stays a small multiple of eps, and n eps
   leaves room for it (infinity norms). *)
let backward_error a x b = ratio (Vec.norm_inf (residual a x b)) ((norm a *. Vec.norm_inf x) +. Vec.norm_inf b)

module LinalgBackward = struct
  type t = Linalg.matrix * Vec.t

  let name = "linalg backward error <= n eps, condition up to 1e8"
  let gen = system 8 1e8
  let show = show_system
  let worst (a, b) = match Linalg.solve a b with Some x -> backward_error a x b /. (dim a *. eps) | None -> infinity
  let holds c = worst c <= 1.
end

(* Both solutions meet the bound above, so |A (x_p - x)| = |r_p - r| <= n eps (|A| (|x| + |x_p|) + 2 |b|). Without
   ties in a pivot column both orders pivot on the same rows and agree exactly; dividing each row by the magnitude
   of its first entry, in half the cases, makes that column all ties, and the orders part ways. *)
module LinalgPermuted = struct
  type t = (Linalg.matrix * Vec.t) * int array * int

  let name = "linalg row permutation agrees within the backward-error bound"

  let gen =
    let open Gen.Syntax in
    let* n, kappa = Gen.pair (Gen.int ~lo:1 ~hi:8) (Gen.log_uniform ~lo:1. ~hi:1e8) in
    Gen.triple
      (Gen.pair (Gen.prescribed (spread n kappa)) (vector_of n entry))
      (vector_of n (Gen.int ~lo:0 ~hi:n))
      (Gen.int ~lo:0 ~hi:1)

  let tied (a, b) =
    let scale i = if a.(i).(0) = 0. then 1. else 1. /. Float.abs a.(i).(0) in
    (Array.mapi (fun i row -> Array.map (( *. ) (scale i)) row) a, Array.mapi (fun i bi -> bi *. scale i) b)

  let systems (ab, keys, ties) =
    let a, b = if ties = 1 then tied ab else ab in
    let order = List.stable_sort (fun i j -> compare keys.(i) keys.(j)) (List.init (Array.length a) Fun.id) in
    let permute v = Array.of_list (List.map (fun i -> v.(i)) order) in
    (a, b, permute a, permute b)

  let show ((_, keys, ties) as c) =
    let a, b, _, _ = systems c in
    Printf.sprintf "%s, row keys %s%s"
      (show_system (a, b))
      (String.concat " " (Array.to_list (Array.map string_of_int keys)))
      (if ties = 1 then ", first column tied" else "")

  let worst c =
    let a, b, pa, pb = systems c in
    match (Linalg.solve a b, Linalg.solve pa pb) with
    | Some x, Some xp ->
        ratio
          (Vec.norm_inf (Vec.sub (residual a xp b) (residual a x b)))
          (dim a *. eps *. ((norm a *. (Vec.norm_inf x +. Vec.norm_inf xp)) +. (2. *. Vec.norm_inf b)))
    | _ -> infinity

  let holds c = worst c <= 1.
end

(* Two equal rows stay equal through every elimination step (same multipliers, same pivot row) until one is the
   pivot; the other then becomes exactly 0 (multiplier 1), and a zero row ends in a zero pivot: None. *)
module LinalgSingular = struct
  type t = Linalg.matrix * Vec.t

  let name = "linalg repeated row gives None"

  let gen =
    let open Gen.Syntax in
    let* n = Gen.int ~lo:2 ~hi:8 in
    Gen.pair (Gen.singular n) (vector_of n entry)

  let show = show_system
  let holds (a, b) = Option.is_none (Linalg.solve a b)
end

(* Scaling by powers of two is exact, or loses digits only far below the smallest singular value when it underflows.
   With A scaled by 2^p and b by 2^q, |p|, |q| <= 960 and |q - p| <= 900, the solution (at most
   kappa sqrt n 2^(q - p)) and every intermediate value stay far from overflow: a non-finite entry can only come from
   the solver. *)
module LinalgFinite = struct
  type t = (Linalg.matrix * Vec.t) * int * int

  let name = "linalg results are finite, scales 2^-960 to 2^960"
  let gen = Gen.triple (system 8 1e8) (signed Int.neg (Gen.int ~lo:0 ~hi:960)) (signed Int.neg (Gen.int ~lo:0 ~hi:900))

  let scaled ((a, b), p, d) =
    let q = Int.max (-960) (Int.min 960 (p + d)) in
    (Array.map (Array.map (fun x -> Float.ldexp x p)) a, Array.map (fun x -> Float.ldexp x q) b)

  let show c = show_system (scaled c)
  let holds c = match Linalg.solve (fst (scaled c)) (snd (scaled c)) with Some x -> Vec.finite x | None -> false
end

(* Jac.forward on f y = M y + c, evaluated in double-double and rounded once, so each f_i is within u |f_i| of the
   exact value, plus the double-double error, at most 4 (n + 1) u^2 S_i with S_i = sum_k |M_ik y_k| + |c_i|. An
   affine map has no truncation error: column j divides f(y + d e_j) - f(y) = M_j d + two roundings of f by the
   step, which is d = 1e-8 (1 + |y_j|) to a relative 1e-7 (the rounding of y_j + d), and the subtraction, the
   reciprocal, the product and the stored step round J_ij four more times. Hence, with f_i(y + d e_j) at most
   |f_i(y)| + |M_ij| d:  |J_ij - M_ij| <= (2 u |f_i(y)| + 8 (n + 1) u^2 S_i) (1 + 1e-7) / d + 5 u |M_ij|. *)
let affine m c y = Array.mapi (fun i row -> Dd.to_float (Dd.add (Dd.dot row y) (Dd.of_float c.(i)))) m

let jac_bound m c y =
  let fy = affine m c y and n = Float.of_int (Array.length y) in
  Array.mapi
    (fun i row ->
      let s = abs_sum (Array.mapi (fun k mik -> mik *. y.(k)) row) +. Float.abs c.(i) in
      Array.mapi
        (fun j mij ->
          let d = 1e-8 *. (1. +. Float.abs y.(j)) in
          (((2. *. u *. Float.abs fy.(i)) +. (8. *. (n +. 1.) *. u *. u *. s)) *. (1. +. 1e-7) /. d)
          +. (5. *. u *. Float.abs mij))
        row)
    m

(* The largest |J_ij - M_ij| over its bound, for J and M of the same shape. *)
let jac_worst j m bound =
  Array.fold_left Float.max 0.
    (Array.mapi
       (fun i row ->
         Array.fold_left Float.max 0. (Array.mapi (fun k jik -> ratio (Float.abs (jik -. m.(i).(k))) bound.(i).(k)) row))
       j)

(* Coordinates of all sizes, so the step is absolute near 0 and relative for large ones. *)
let coordinate = Gen.one_of [ Gen.return 0.; magnitude ~lo:1e-12 ~hi:1e8 ]
let show_affine (m, c, y) = Printf.sprintf "m = %s, c = %s, y = %s" (matrix m) (vector c) (vector y)

module JacAffine = struct
  type t = Linalg.matrix * Vec.t * Vec.t

  let name = "jac affine map within the forward-difference rounding bound"

  let gen =
    let open Gen.Syntax in
    let* n = Gen.int ~lo:1 ~hi:8 in
    Gen.triple (vector_of n (vector_of n entry)) (vector_of n (magnitude ~lo:1e-3 ~hi:1e3)) (vector_of n coordinate)

  let show = show_affine
  let worst (m, c, y) = jac_worst (Jac.forward (affine m c) y) m (jac_bound m c y)
  let holds c = worst c <= 1.
end

(* f from R^n to R^m: J must have a row per output and a column per input. Moving the bottom-left entry of M by 3
   puts a square M of size 2 or more far from symmetric, so its transpose misses the bound by orders of magnitude. *)
module JacOrientation = struct
  type t = Linalg.matrix * Vec.t * Vec.t

  let name = "jac orientation: n inputs and m outputs give m rows, M and not its transpose"

  let gen =
    let open Gen.Syntax in
    let* m, n = Gen.pair (Gen.int ~lo:1 ~hi:6) (Gen.int ~lo:1 ~hi:6) in
    let+ a, c, y = Gen.triple (vector_of m (vector_of n entry)) (vector_of m entry) (vector_of n coordinate) in
    ( Array.mapi (fun i row -> if i = m - 1 then Array.mapi (fun j x -> if j = 0 then x +. 3. else x) row else row) a,
      c,
      y )

  let show = show_affine

  let holds (m, c, y) =
    let j = Jac.forward (affine m c) y and bound = jac_bound m c y in
    let rows = Array.length m and cols = Array.length y in
    Array.length j = rows
    && Array.for_all (fun row -> Array.length row = cols) j
    && jac_worst j m bound <= 1.
    && (rows <> cols || cols = 1
       || jac_worst j (Array.init cols (fun i -> Array.init cols (fun k -> m.(k).(i)))) bound > 1.)
end

(* With M_ij = 0, f_i does not depend on y_j: its double-double sum adds an exact 0 for that term, so f_i(y + d e_j)
   equals f_i(y) and J_ij is exactly 0. A perturbation that reached any other coordinate would show there. *)
module JacColumn = struct
  type t = Linalg.matrix * Vec.t * Vec.t

  let name = "jac perturbing y_j moves only column j"

  let gen =
    let open Gen.Syntax in
    let* n = Gen.int ~lo:1 ~hi:8 in
    Gen.triple
      (vector_of n (vector_of n (Gen.one_of [ Gen.return 0.; entry ])))
      (vector_of n entry) (vector_of n coordinate)

  let show = show_affine

  let holds (m, c, y) =
    let j = Jac.forward (affine m c) y in
    Array.for_all2 (fun jrow mrow -> Array.for_all2 (fun jik mik -> mik <> 0. || jik = 0.) jrow mrow) j m
end

(* Newton on f x = 2^e (M x + a x^3) - c (cubes by entry), with the Jacobian j 2^e (M + 3 a diag(x^2)): j = 1 is
   exact, any other j makes Newton converge at best linearly. *)
type problem = { m : Linalg.matrix; a : float; e : int; c : Vec.t; x0 : Vec.t; j : float }

let newton_f p x =
  Array.mapi (fun i row -> Float.ldexp (Vec.dot row x +. (p.a *. x.(i) *. x.(i) *. x.(i))) p.e -. p.c.(i)) p.m

let newton_jac p x =
  Array.mapi
    (fun i row ->
      Array.mapi (fun k mik -> p.j *. Float.ldexp (if i = k then mik +. (3. *. p.a *. x.(i) *. x.(i)) else mik) p.e) row)
    p.m

let newton p = Newton.solve (newton_f p) (newton_jac p) p.x0

let show_problem p =
  Printf.sprintf "m = %s, a = %s, e = %d, c = %s, x0 = %s, j = %s" (matrix p.m) (num p.a) p.e (vector p.c) (vector p.x0)
    (num p.j)

(* M = A^T A is symmetric positive definite with eigenvalues in [0.01, 1], so for a >= 0 f is strongly monotone and
   its Jacobian is never singular; c and x0 lie in [-reach, reach]. *)
let monotone ~reach factor =
  let open Gen.Syntax in
  let* n, kappa = Gen.pair (Gen.int ~lo:1 ~hi:6) (Gen.log_uniform ~lo:1. ~hi:10.) in
  let start = vector_of n (Gen.float ~lo:(-.reach) ~hi:reach) in
  let+ a, (cubic, j), (c, x0) =
    Gen.triple (Gen.prescribed (spread n kappa)) (Gen.pair (Gen.float ~lo:0. ~hi:1.) factor) (Gen.pair start start)
  in
  let gram =
    Array.init n (fun i ->
        Array.init n (fun k -> Array.fold_left ( +. ) 0. (Array.map (fun row -> row.(i) *. row.(k)) a)))
  in
  { m = gram; a = cubic; e = 0; c; x0; j }

(* Newton's own checks (newton.mli) are all that stands between an overflow and an Ok. Half the cases are wild:
   values anywhere in the range or special, scales from 2^-1000 to 2^1000, Jacobians off by up to a factor 4. The
   other half put a root next to max_float and start within a relative 1e-6 of it, so that a step short enough to
   pass the stopping test can still overshoot past max_float when the Jacobian is too small. *)
module NewtonFinite = struct
  type t = problem

  let name = "newton Ok is finite"
  let near_max = Gen.map (fun d -> Float.max_float *. (1. -. d)) (Gen.log_uniform ~lo:1e-16 ~hi:1e-6)
  let wild = Gen.one_of [ magnitude ~lo:1e-300 ~hi:1e300; Gen.special; signed Float.neg near_max ]
  let factor = Gen.log_uniform ~lo:0.25 ~hi:4.

  let anywhere =
    let open Gen.Syntax in
    let* n, kappa = Gen.pair (Gen.int ~lo:1 ~hi:3) (Gen.log_uniform ~lo:1. ~hi:1e4) in
    let+ (m, a, e), (c, x0, j) =
      Gen.pair
        (Gen.triple
           (Gen.prescribed (spread n kappa))
           (Gen.one_of [ Gen.return 0.; Gen.float ~lo:0. ~hi:1. ])
           (signed Int.neg (Gen.int ~lo:0 ~hi:1000)))
        (Gen.triple (vector_of n wild) (vector_of n wild) factor)
    in
    { m; a; e; c; x0; j }

  let next_to_max =
    let open Gen.Syntax in
    let* n = Gen.int ~lo:1 ~hi:2 in
    let root = Gen.one_of [ magnitude ~lo:1e-3 ~hi:1e3; signed Float.neg near_max ] in
    let+ m, (x, offsets), j =
      Gen.triple
        (Gen.prescribed (spread n 10.))
        (Gen.pair (vector_of n root) (vector_of n (magnitude ~lo:1e-16 ~hi:1e-6)))
        factor
    in
    {
      m;
      a = 0.;
      e = 0;
      c = Array.map (fun row -> Vec.dot row x) m;
      x0 = Array.mapi (fun i xi -> xi *. (1. +. offsets.(i))) x;
      j;
    }

  let gen = Gen.one_of [ anywhere; next_to_max ]
  let show = show_problem
  let holds p = match newton p with Ok x -> Vec.finite x | Error _ -> true
end

(* The stopping test of newton.mli, evaluated again at an Ok: f x = 0, or the next step has
   |dx|_inf <= 1e-10 (1 + |x|_inf). A Jacobian off by a factor j in [0.7, 2] makes Newton converge linearly with
   rate |1 - 1/j| <= 1/2, so the step after the one that stopped it is about half of that one, at most: the test
   passes again within a factor 2, and a looser tolerance inside Newton would fail it. *)
module NewtonStops = struct
  type t = problem

  let name = "newton Ok passes its own stopping test again"
  let gen = monotone ~reach:10. (Gen.float ~lo:0.7 ~hi:2.)
  let show = show_problem

  let worst p =
    match newton p with
    | Error _ -> 0.
    | Ok x -> (
        let fx = newton_f p x in
        if Vec.norm_inf fx = 0. then 0.
        else
          match Linalg.solve (newton_jac p x) (Vec.scale (-1.) fx) with
          | Some dx -> Vec.norm_inf dx /. (1e-10 *. (1. +. Vec.norm_inf x))
          | None -> infinity)

  let holds p = worst p <= 1.
end

(* On f x = M x - c the first step solves the system and the next is about kappa eps |x|, far below the tolerance for
   kappa <= 1e4: Newton returns Ok after a solve and one correction, with the backward error of a solve, below n eps. *)
module NewtonLinear = struct
  type t = problem

  let name = "newton converges on linear problems, condition up to 1e4"

  let gen =
    let open Gen.Syntax in
    let* m, c = system 8 1e4 in
    let+ x0 = vector_of (Array.length m) (Gen.float ~lo:(-10.) ~hi:10.) in
    { m; a = 0.; e = 0; c; x0; j = 1. }

  let show = show_problem
  let worst p = match newton p with Ok x -> backward_error p.m x p.c /. (dim p.m *. eps) | Error _ -> infinity
  let holds p = worst p <= 1.
end

(* With an exact Jacobian the step is a descent direction for every norm of f, yet from starts as far as 1e4 a flat
   direction of M (eigenvalue near 0.01) can make it so long that the cubic defeats every damping down to 1/1024,
   and Newton returns Diverged. Convergence is likely, not certain: a statistical claim. *)
module NewtonConverges = struct
  type t = problem

  let name = "newton converges from random starts on monotone cubics"
  let threshold = 0.99
  let gen = monotone ~reach:1e4 (Gen.return 1.)
  let show = show_problem
  let holds p = Result.is_ok (newton p)
end

(* T1: with h = omega h_prev the formula is exact on quadratics. In units of h_prev the nodes are -1, 0 and omega,
   so P(s) = c0 + c1 s + c2 s^2 must satisfy P(omega) = a1 P(0) + a0 P(-1) + beta omega P'(omega). Double-double
   evaluates the defect so closely that only the rounding of Bdf2.coeffs remains: to first order at most 5
   roundings of a1, 3 of a0 and 3 of beta, so |defect| <= u (5 |a1 P(0)| + 3 |a0 P(-1)| + 3 |beta omega P'(omega)|). *)
let steps =
  Gen.map
    (fun (h_prev, r) -> (h_prev, h_prev *. r))
    (Gen.pair (Gen.log_uniform ~lo:1e-6 ~hi:1.)
       (Gen.one_of [ Gen.float ~lo:0.01 ~hi:2.2; Gen.log_uniform ~lo:1e-6 ~hi:2.2 ]))

let show_steps (h_prev, h) = Printf.sprintf "h_prev = %s, h = %s" (num h_prev) (num h)

module Bdf2Exact = struct
  type t = (float * float) * (float * float * float)

  let name = "bdf2 exact on quadratics up to the rounding of the coefficients (T1)"
  let gen = Gen.pair steps (Gen.triple entry entry entry)
  let show (hs, (c0, c1, c2)) = Printf.sprintf "%s, P = %s + %s s + %s s^2" (show_steps hs) (num c0) (num c1) (num c2)

  let worst ((h_prev, h), (c0, c1, c2)) =
    let omega = h /. h_prev in
    let { Bdf2.a1; a0; beta } = Bdf2.coeffs omega in
    let open Dd in
    let p s = add (of_float c0) (add (mul (of_float c1) s) (mul (of_float c2) (mul s s))) in
    let w = of_float omega in
    let now = mul (of_float a1) (p (of_float 0.)) and before = mul (of_float a0) (p (of_float (-1.))) in
    let slope = mul (of_float beta) (mul w (add (of_float c1) (mul (of_float (2. *. c2)) w))) in
    let mag x = Float.abs (to_float x) in
    ratio
      (mag (sub (p w) (add now (add before slope))))
      (u *. ((5. *. mag now) +. (3. *. mag before) +. (3. *. mag slope)) *. (1. +. 1e-12))

  let holds c = worst c <= 1.
end

(* T1 with p = 1 is a1 + a0 = 1. In Bdf2.coeffs both share one rounded denominator, whose rounding cancels to first
   order in the sum (a1 + a0 = 1); a1 has 4 roundings more and a0 2, so |a1 + a0 - 1| <= u (4 |a1| + 2 |a0| + 1),
   about 5 ulps of 1 at omega = 2.2. The sum is exact in double-double. *)
module Bdf2Sum = struct
  type t = float * float

  let name = "bdf2 |a1 + a0 - 1| within a few ulps (T1)"
  let gen = steps
  let show = show_steps

  let worst (h_prev, h) =
    let { Bdf2.a1; a0; _ } = Bdf2.coeffs (h /. h_prev) in
    let s = Dd.add (Dd.of_float a1) (Dd.of_float a0) in
    ratio (Float.abs (s.hi -. 1. +. s.lo)) (u *. ((4. *. Float.abs a1) +. (2. *. Float.abs a0) +. 1.) *. (1. +. 1e-12))

  let holds c = worst c <= 1.
end

(* T2: since a1 + a0 = 1 the homogeneous recurrence y_{n+2} = a1 y_{n+1} + a0 y_n has d_{n+1} = -a0 d_n for
   d_n = y_{n+1} - y_n, and -a0 grows with omega and stays below 1 for omega < 1 + sqrt 2: with every ratio at most w
   the differences contract by q(w) per step. Rounding adds (a1 + a0 - 1) y_{n+1}, at most
   u (4 |a1| + 2 |a0| + 1) |y_{n+1}| (Bdf2Sum); the three roundings of the recurrence,
   u (|a1 y_{n+1}| + |a0 y_n| + |y_{n+2}|); and those of a0, q(w), the differences and this test, 10 u q(w) |d_n|.
   q(w) is computed in closed form, apart from Bdf2.coeffs. *)
let contraction w = w *. w /. (1. +. (2. *. w))

module Bdf2Contracts = struct
  type t = float * (float * float) * float array

  let name = "bdf2 differences contract by q(w) when every ratio is at most w < 1 + sqrt 2 (T2)"

  (* Each ratio is w times a share in (0, 1], w itself first. *)
  let gen =
    Gen.triple (Gen.float ~lo:0.05 ~hi:2.414) (Gen.pair entry entry)
      (Gen.array (Gen.int ~lo:1 ~hi:40) (Gen.one_of [ Gen.return 1.; Gen.float ~lo:1e-3 ~hi:1. ]))

  let show (w, (y0, y1), shares) =
    Printf.sprintf "w = %s, y0 = %s, y1 = %s, ratios w times %s" (num w) (num y0) (num y1) (vector shares)

  let worst (w, (y0, y1), shares) =
    let q = contraction w in
    let step (worst, y_prev, y) share =
      let { Bdf2.a1; a0; _ } = Bdf2.coeffs (w *. share) in
      let y_next = (a1 *. y) +. (a0 *. y_prev) in
      let d = Float.abs (y -. y_prev) and y_abs = Float.abs y in
      let rounding =
        (((4. *. Float.abs a1) +. (2. *. Float.abs a0) +. 1.) *. y_abs)
        +. Float.abs (a1 *. y)
        +. Float.abs (a0 *. y_prev)
      in
      let slack = u *. ((10. *. q *. d) +. rounding +. Float.abs y_next) in
      (Float.max worst (ratio (Float.abs (y_next -. y)) ((q *. d) +. slack)), y, y_next)
    in
    let worst, _, _ = Array.fold_left step (0., y0, y1) shares in
    worst

  let holds c = worst c <= 1.
end

(* The run seed and the multiplier on the counts: fixed for dune runtest, set from the environment for nightly runs,
   and read once, here. *)
let setting var default valid =
  match Sys.getenv_opt var with
  | None -> Ok default
  | Some s -> (
      match int_of_string_opt s with
      | Some n when valid n -> Ok n
      | _ -> Error (Printf.sprintf "%s=%s is not a valid setting" var s))

let suite ~seed ~scale =
  let prop (module P : Prop.Property) ~count =
    ("prop " ^ P.name, fun () -> Prop.to_string (Prop.check (module P) ~count:(count * scale) ~seed))
  in
  let share (module P : Prop.Statistical) ~count =
    ("prop " ^ P.name, fun () -> Prop.to_string (Prop.check_statistical (module P) ~count:(count * scale) ~seed))
  in
  [
    prop (module VecNormPermuted) ~count:20000;
    prop (module VecFinite) ~count:20000;
    prop (module VecAxpy) ~count:20000;
    prop (module LinalgBackward) ~count:10000;
    prop (module LinalgPermuted) ~count:10000;
    prop (module LinalgSingular) ~count:20000;
    prop (module LinalgFinite) ~count:10000;
    prop (module JacAffine) ~count:10000;
    prop (module JacOrientation) ~count:10000;
    prop (module JacColumn) ~count:10000;
    prop (module NewtonFinite) ~count:20000;
    prop (module NewtonStops) ~count:10000;
    prop (module NewtonLinear) ~count:10000;
    share (module NewtonConverges) ~count:10000;
    prop (module Bdf2Exact) ~count:20000;
    prop (module Bdf2Sum) ~count:20000;
    prop (module Bdf2Contracts) ~count:20000;
  ]

let () =
  match (setting "VSTIFF_PROP_SEED" 2026 (fun _ -> true), setting "VSTIFF_PROP_SCALE" 1 (fun n -> n > 0)) with
  | Ok seed, Ok scale -> Report.lines (suite ~seed ~scale)
  | Error e, _ | _, Error e -> Report.lines [ ("prop settings", fun () -> e) ]
