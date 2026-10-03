(** Adaptive integration: method [M] proposes each step with an error
    estimate, controller [C] keeps or rejects it and picks the next size. *)

type 'stats solution = { t : float; y : Vec.t; stats : 'stats }

let integrate (module M : Ode.Embedded) (module C : Ode.Controller) ?dt0 ?dt_max ?(max_rejects = 50) ~tol
    (p : Ode.problem) : (C.stats solution, Fail.t) result =
  let span = p.t_end -. p.t0 in
  let dt_max = Option.value dt_max ~default:(span /. 10.) in
  let dt0 = Float.min dt_max (Option.value dt0 ~default:(1e-6 *. span)) in
  Check.adaptive ~dt0 p;
  let rec go c history (at : Ode.point) =
    if at.t >= p.t_end then Ok { t = at.t; y = at.y; stats = C.stats c }
    else
      (* The last step is cut to land exactly on t_end. *)
      let dt = C.proposal c in
      let last = dt >= p.t_end -. at.t in
      let h = if last then p.t_end -. at.t else dt in
      let retry reason = Result.bind (C.rejected c reason ~at:at.t ~h) (fun c -> go c history at) in
      match M.step_with_error p.rhs h history at with
      | Ok (y, err, next) when C.acceptable c ~y ~err ->
          go (C.accepted c) next { Ode.t = (if last then p.t_end else at.t +. h); y }
      | Ok _ -> retry Ode.Too_large
      | Error e -> retry (Ode.Solver e)
  in
  if not (Vec.finite p.y0 && Vec.finite (p.rhs p.t0 p.y0)) then Error Fail.Nan
  else go (C.init ~tol ~dt0 ~dt_max ~max_rejects) M.start { t = p.t0; y = p.y0 }
