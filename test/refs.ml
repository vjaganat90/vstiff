(** Hardcoded reference values, computed outside vstiff. *)

(** Robertson y1 at t = 1e4 from y(0) = (1, 0, 0): scipy 1.18 solve_ivp with
    Radau (rtol 1e-13), BDF and LSODA (rtol 1e-12), atol 1e-20 and
    finite-difference Jacobians. The three agree to 4e-12. *)
let robertson_y1_at_1e4 = 0.10730042854
