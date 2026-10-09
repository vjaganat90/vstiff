(* The loop owns the point and the history; [C] owns the step size and its
   bookkeeping. A rejection leaves the point and the history alone and changes
   only [C]. Theory: docs/numerics/05-step-control.md. *)

module Vec = Numerics.Vec

type 'stats solution = { t : float; y : float array; stats : 'stats }

(* (module C : Ode.Controller) is a modular explicit: the result type names
   C.stats (docs/ocaml.md). *)
let integrate (module M : Ode.Embedded) (module C : Ode.Controller) ?dt0 ?dt_max ?(max_rejects = 50) ~tol
    (p : Ode.problem) : (C.stats solution, Fail.t) result =
  let span = p.t_end -. p.t0 in
  (* ?dt_max arrives as an option; Option.value supplies the default (docs/ocaml.md). *)
  (* A default that underflows to 0 on a tiny span falls back to the span. *)
  let positive d = if d > 0. then d else span in
  let dt_max = Option.value dt_max ~default:(positive (span /. 10.)) in
  (* Start small: how fast y changes is not known yet. Never above dt_max. *)
  let dt0 = Float.min dt_max (Option.value dt0 ~default:(positive (1e-6 *. span))) in
  Check.adaptive ~dt0 ~tol p;
  let f0 = p.rhs p.t0 p.y0 in
  Check.output ~caller:"Adaptive.integrate" p f0;
  (* go loops: every call to go is a tail call, so the stack does not grow (docs/ocaml.md). *)
  let rec go c history (at : Ode.point) =
    if at.t >= p.t_end then Ok { t = at.t; y = at.y; stats = C.stats c }
    else
      let dt = C.proposal c and remaining = p.t_end -. at.t in
      (* The last step lands on t_end itself, and so does a remainder below t's
         resolution. Other steps are snapped to the floats, so the state moves by
         exactly what the clock does, and one below the resolution is rejected:
         above it snapping changes a step by at most 1/32, so BDF2's step ratio
         stays near 2, inside its limit 1 + sqrt 2. *)
      let last = dt >= remaining || remaining <= Clock.resolution at.t in
      let t_next = if last then p.t_end else at.t +. dt in
      let h = if last then remaining else t_next -. at.t in
      let retry reason = Result.bind (C.rejected c reason ~at:at.t ~h) (fun c -> go c history at) in
      if h <= 0. || ((not last) && h < Clock.resolution at.t) then retry Ode.Too_small
      else
        match M.step_with_error p.rhs h history at with
        | Ok (y, err, next) when C.acceptable c ~y ~err -> go (C.accepted c) next { Ode.t = t_next; y }
        | Ok _ -> retry Ode.Too_large
        | Error e -> retry (Ode.Solver e)
  in
  (* Without this check a non-finite start would only end as StepRejected. *)
  if not (Vec.finite p.y0 && Vec.finite f0) then Error Fail.Nan
  else go (C.init ~tol ~dt0 ~dt_max ~max_rejects) M.start { t = p.t0; y = p.y0 }
