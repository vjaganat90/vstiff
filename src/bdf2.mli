(** Variable-step BDF2, second order, with backward Euler for the first step
    and as the companion of its error estimate. *)

(** The step solves [y_{n+1} = a1 y_n + a0 y_{n-1} + beta h f(t_{n+1}, y_{n+1})]. *)
type coeffs = { a1 : float; a0 : float; beta : float }

(** [coeffs omega] for the step ratio [omega = h / h_prev]: the weights that
    make the formula exact on quadratics through the last three points. *)
val coeffs : float -> coeffs

(** [step_with_error] estimates the local error by the gap to backward Euler
    from the same point. *)
include Ode.Embedded
