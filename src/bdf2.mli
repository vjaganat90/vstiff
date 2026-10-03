(** Variable-step BDF2, second order, with backward Euler for the first step
    and as the companion of its error estimate. Zero-stable (earlier errors stay
    bounded) for step ratios [h / h_prev] below [1 + sqrt 2], about 2.414; the
    controller must keep to that. See docs/numerics/04-bdf.md. *)

(** The step solves [y_{n+1} = a1 y_n + a0 y_{n-1} + beta h f(t_{n+1}, y_{n+1})]. *)
type coeffs = { a1 : float; a0 : float; beta : float }

(** [coeffs omega] for the step ratio [omega = h / h_prev]: the weights that
    make the formula exact on quadratics through the last three points.
    Evaluate it to see them; docs/numerics/04-bdf.md shows how to derive them. *)
val coeffs : float -> coeffs

(** [step_with_error] estimates the local error by the gap to backward Euler
    from the same point. The first step is backward Euler, so its estimate is
    half its gap to explicit Euler. See docs/numerics/05-step-control.md. *)
include Ode.Embedded
