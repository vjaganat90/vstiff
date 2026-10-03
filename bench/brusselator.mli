(** The Brusselator with diffusion on (0, 1), by the method of lines on 40 interior cells (n = 80), to t = 10:
    u' = A + u^2 v - (B + 1) u + alpha u_xx and v' = B u - u^2 v + alpha v_xx, with A = 1, B = 3, alpha = 1/50, u = 1
    and v = 3 at both ends, u(x, 0) = 1 + sin(2 pi x), v(x, 0) = 3. The state interleaves the cells, (u_1, v_1, u_2,
    ...), so the Jacobian is banded; vstiff's forward differences still call the right-hand side 81 times for it. *)

(** The right-hand side, spelled as in bench/compare/problems.py. *)
val rhs : Vstiff.Ode.rhs

val y0 : float array
val problem : Vstiff.Ode.problem
