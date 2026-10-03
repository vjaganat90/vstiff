(** Backward Euler, [y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})]: first order,
    and stable for any step on decaying problems. It needs no history. *)

include Ode.Method
