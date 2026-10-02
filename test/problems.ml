(** The corpus problems. Each exposes only its right-hand side, initial data
    and, where one exists, the exact solution. *)

(** y' = -Λ y with Λ = diag(1, 100, 1e4): a swapped Jacobian row puts a fast
    rate on a slow component and wrecks the answer. *)
module Canary = struct
  let lambda = [| 1.; 100.; 1e4 |]
  let rhs _t y = Array.mapi (fun i yi -> -.lambda.(i) *. yi) y
  let y0 = [| 1.; 1.; 1. |]
  let exact t = Array.map (fun l -> exp (-.l *. t)) lambda
end
