open Vstiff

type verdict = Within_band | Regressed of string | Unpinned

let cost_limit_percent = 10
let scd_slack = 0.05

let figures_of (row : Measure.t) = match row.outcome with Measure.Done f -> Some f | Failed _ -> None

let pinned_figures pinned (row : Measure.t) =
  List.find_map
    (fun (p : Measure.t) -> if p.problem = row.problem && p.tol = row.tol then figures_of p else None)
    pinned

(* Integers decide the cost, so a row at exactly +10 % passes and one call more does not, whatever the rounding. *)
let verdict pinned (row : Measure.t) =
  match pinned_figures pinned row with
  | None -> Unpinned
  | Some golden -> (
      match row.outcome with
      | Failed e -> Regressed ("failed with " ^ Fail.to_string e)
      | Done now ->
          let costly = now.rhs_calls * 100 > golden.rhs_calls * (100 + cost_limit_percent) in
          let less_accurate = now.scd < golden.scd -. scd_slack in
          let extra = 100. *. float_of_int (now.rhs_calls - golden.rhs_calls) /. float_of_int golden.rhs_calls in
          let reasons =
            (if costly then [ Printf.sprintf "rhs calls %+.1f%%" extra ] else [])
            @ if less_accurate then [ Printf.sprintf "scd %.2f, pinned %.2f" now.scd golden.scd ] else []
          in
          if reasons = [] then Within_band else Regressed (String.concat "; " reasons))

let describe = function
  | Within_band -> "ok"
  | Regressed what -> "REGRESSED: " ^ what
  | Unpinned -> "UNPINNED: run bench --pin"

let header =
  "(* The golden table: the figures of one run of bench.exe, written by bench.exe --pin and never edited by hand\n\
  \   (bench/README.md). Gate judges later runs against it. *)\n\n\
   open Measure\n\n\
   let rows =\n\
  \  ["

(* An infinite scd, from an exact answer, has no float literal of the %.2f form: such a row is not pinned. *)
let row_source (row : Measure.t) =
  match figures_of row with
  | Some f when Float.is_finite f.error && Float.is_finite f.scd ->
      Some
        (Printf.sprintf
           "    {\n\
           \      problem = %S;\n\
           \      tol = %g;\n\
           \      outcome = Done { steps = %d; rejected = %d; rhs_calls = %d; error = %.3e; scd = %.2f };\n\
           \    };"
           row.problem row.tol f.steps f.rejected f.rhs_calls f.error f.scd)
  | _ -> None

(* No newline at the end: Console.print adds the one. *)
let source rows =
  let lines = List.map row_source rows in
  if List.mem None lines then None else Some (String.concat "\n" ((header :: List.filter_map Fun.id lines) @ [ "  ]" ]))
