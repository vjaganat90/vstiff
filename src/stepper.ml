(* Fixed-step driver; the method is passed as a modular explicit (docs/ocaml.md). *)

open Fail.Syntax

let fixed (module M : Ode.Method) ~dt (p : Ode.problem) =
  Check.fixed ~dt p;
  let span = p.t_end -. p.t0 in
  if span = 0. then Ok p.y0
  else
    let n = max 1 (Float.to_int (Float.round (span /. dt))) in
    let h = span /. float_of_int n in
    let rec go left history (at : Ode.point) =
      if left = 0 then Ok at.y
      else
        let* y, history = M.step p.rhs h history at in
        go (left - 1) history { t = at.t +. h; y }
    in
    go n M.start { t = p.t0; y = p.y0 }
