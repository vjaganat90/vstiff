open Vstiff

type entry = { row : Measure.t; cpu : float option; verdict : Gate.verdict }

let cpu_cell = function Some s -> Printf.sprintf " %8.2f" s | None -> ""

let title =
  "vstiff: Adaptive.integrate with Bdf2 and Halving, rtol = atol = tol, other arguments at their defaults; rhs calls\n\
   include the finite-difference Jacobians"

let table_header (options : Options.t) =
  Printf.sprintf "%-15s %-6s %8s %9s %12s %10s %6s%s  %s" "problem" "tol" "steps" "rejected" "rhs calls" "error" "scd"
    (if options.cpu then Printf.sprintf " %8s" "cpu s" else "")
    "verdict"

let table_line { row; cpu; verdict } =
  match row.outcome with
  | Measure.Done f ->
      Printf.sprintf "%-15s %-6.0e %8d %9d %12d %10.2e %6.2f%s  %s" row.problem row.tol f.steps f.rejected f.rhs_calls
        f.error f.scd (cpu_cell cpu) (Gate.describe verdict)
  | Failed e ->
      Printf.sprintf "%-15s %-6.0e failed: %s%s  %s" row.problem row.tol (Fail.to_string e) (cpu_cell cpu)
        (Gate.describe verdict)

let csv_header (options : Options.t) =
  "problem,solver,rtol,steps,rejected,rhs,nfev,njev,error,scd" ^ if options.cpu then ",cpu_s" else ""

let csv_line { row; cpu; verdict = _ } =
  let figures =
    match row.outcome with
    | Measure.Done f -> Printf.sprintf "%d,%d,%d,,,%.3e,%.2f" f.steps f.rejected f.rhs_calls f.error f.scd
    | Failed _ -> ",,,,,,"
  in
  Printf.sprintf "%s,bdf2+halving,%g,%s%s" row.problem row.tol figures
    (match cpu with Some s -> Printf.sprintf ",%.2f" s | None -> "")
