(** Forward-difference Jacobians. *)

(** [forward f y] approximates [J.(i).(j) = d f_i / d y_j] at [y], perturbing
    [y_j] by [1e-8 (1 + |y_j|)]. Costs [n + 1] calls of [f]. *)
val forward : (Vec.t -> Vec.t) -> Vec.t -> Linalg.matrix
