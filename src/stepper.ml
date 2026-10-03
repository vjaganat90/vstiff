(* Fixed-step driver; the method is passed as a modular explicit (docs/ocaml.md). *)

open Fail.Syntax

let fixed (module M : Ode.Method) ~dt (p : Ode.problem) =
  Check.fixed ~dt p;
  let span = p.t_end -. p.t0 in
  if span = 0. then Ok p.y0
  else
    let n = max 1 (Float.to_int (Float.round (span /. dt))) in
    let h = span /. float_of_int n in
    (* Step k ends at t0 + k h, the last one at t_end itself, and each step is
       the difference of its end times: the clock cannot drift, and the state
       advances by exactly what the clock does. *)
    let time k = if k = n then p.t_end else p.t0 +. (float_of_int k *. h) in
    let rec go k history (at : Ode.point) =
      if k = n then Ok at.y
      else
        let t_next = time (k + 1) in
        let* y, history = M.step p.rhs (t_next -. at.t) history at in
        go (k + 1) history { t = t_next; y }
    in
    go 0 M.start { t = p.t0; y = p.y0 }
