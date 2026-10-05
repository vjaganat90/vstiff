(* The pure parts of the bench, checked like the corpus: a table of named cases, one line each printed by Report (the
   corpus's own module, copied by the dune file), which dune diffs with check.expected. It runs in dune runtest, in a
   fraction of a second; the bench itself does not. *)
open Benchlib

let figures ~rhs ~scd = { Measure.steps = 100; rejected = 10; rhs_calls = rhs; error = 1e-6; scd }
let row ~rhs ~scd = { Measure.problem = "p"; tol = 1e-6; outcome = Done (figures ~rhs ~scd) }
let failed fail = { Measure.problem = "p"; tol = 1e-6; outcome = Failed fail }
let pinned = [ row ~rhs:1000 ~scd:4.0 ]
let verdict run = Gate.describe (Gate.verdict pinned run)

(* error and scd are the two measures of the bench; each line has a hand-computable answer. *)
let accuracy =
  let reference = [| 1.; 100. |] in
  [
    ( "error of (1.1, 100) against (1, 100)",
      fun () -> Printf.sprintf "%.3e" (Measure.error ~reference [| 1.1; 100. |]) );
    ( "scd of (1.001, 100) against (1, 100)",
      fun () -> Printf.sprintf "%.2f" (Measure.scd ~reference [| 1.001; 100. |]) );
    ("scd of an exact answer", fun () -> Printf.sprintf "%g" (Measure.scd ~reference [| 1.; 100. |]));
    ( "scd with a zero reference component, off by 1e-3",
      fun () -> Printf.sprintf "%.2f" (Measure.scd ~reference:[| 0.; 1. |] [| 1e-3; 1. |]) );
    ( "error with a nan state is nan",
      fun () -> string_of_bool (Float.is_nan (Measure.error ~reference [| Float.nan; 100. |])) );
  ]

(* The gate: a band of +10 % in rhs calls, compared in integers, and 0.05 digits of scd. *)
let gate =
  [
    ("a row equal to its pin", fun () -> verdict (row ~rhs:1000 ~scd:4.0));
    ("exactly 10 % more rhs calls", fun () -> verdict (row ~rhs:1100 ~scd:4.0));
    ("one call more than that", fun () -> verdict (row ~rhs:1101 ~scd:4.0));
    ("fewer calls and more digits", fun () -> verdict (row ~rhs:600 ~scd:5.0));
    ("0.04 digits lost", fun () -> verdict (row ~rhs:1000 ~scd:3.96));
    ("0.06 digits lost", fun () -> verdict (row ~rhs:1000 ~scd:3.94));
    ("both worse", fun () -> verdict (row ~rhs:1500 ~scd:3.0));
    ("a failed run", fun () -> verdict (failed (Vstiff.Fail.StepRejected 51)));
    ("a tolerance that is not pinned", fun () -> verdict { (row ~rhs:1000 ~scd:4.0) with tol = 1e-8 });
    ("a problem that is not pinned", fun () -> verdict { (row ~rhs:1000 ~scd:4.0) with problem = "q" });
  ]

(* The pin is the source of golden.ml: its rows are what bench --pin prints, and a failure cannot be pinned. *)
let pin =
  let rows_of text =
    String.concat " | " (List.filteri (fun i _ -> 5 <= i && i < 12) (String.split_on_char '\n' text))
  in
  [
    ("source of one row", fun () -> match Gate.source pinned with Some text -> rows_of text | None -> "none");
    ( "source with a failed row",
      fun () -> match Gate.source (failed Vstiff.Fail.Nan :: pinned) with Some _ -> "some" | None -> "none" );
    ( "source with an exact answer, scd infinite",
      fun () -> match Gate.source [ row ~rhs:1000 ~scd:Float.infinity ] with Some _ -> "some" | None -> "none" );
  ]

(* The CSV has the columns of bench/compare/scipy_rows.py: a line has as many fields as the header, with or without the
   CPU time, and a failed run keeps the count. *)
let csv =
  let options cpu = { Options.mode = Csv; cpu; only = None } in
  let entry cpu row = { Render.row; cpu; verdict = Gate.Within_band } in
  let fields line = List.length (String.split_on_char ',' line) in
  let line = entry None (row ~rhs:1000 ~scd:4.0) in
  [
    ("a row", fun () -> Render.csv_line line);
    ("a row with its CPU time", fun () -> Render.csv_line { line with cpu = Some 1.5 });
    ( "fields of the header, a row and a failed run",
      fun () ->
        Printf.sprintf "%d, %d, %d"
          (fields (Render.csv_header (options false)))
          (fields (Render.csv_line line))
          (fields (Render.csv_line (entry None (failed Vstiff.Fail.Nan)))) );
    ( "fields with the CPU time",
      fun () ->
        Printf.sprintf "%d, %d"
          (fields (Render.csv_header (options true)))
          (fields (Render.csv_line { line with cpu = Some 1.5 })) );
  ]

let parse args =
  match Options.parse args with
  | Ok { mode; cpu; only } ->
      let mode =
        match mode with
        | Table -> "table"
        | Csv -> "csv"
        | Pin -> "pin"
        | Transcription -> "transcription"
        | Help -> "help"
      in
      Printf.sprintf "%s, cpu %b, only %s" mode cpu (Option.value only ~default:"none")
  | Error message -> "error: " ^ message

(* The command line: one mode at most, and a pin that is whole. *)
let options =
  [
    ("no flag", fun () -> parse []);
    ("--csv --cpu --only hires", fun () -> parse [ "--csv"; "--cpu"; "--only"; "hires" ]);
    ("--pin", fun () -> parse [ "--pin" ]);
    ("--pin --only hires", fun () -> parse [ "--pin"; "--only"; "hires" ]);
    ("--pin --cpu", fun () -> parse [ "--pin"; "--cpu" ]);
    ("--transcription", fun () -> parse [ "--transcription" ]);
    ("--transcription --only hires", fun () -> parse [ "--transcription"; "--only"; "hires" ]);
    ("--csv --pin", fun () -> parse [ "--csv"; "--pin" ]);
    ("--only without a name", fun () -> parse [ "--only" ]);
    ("--only twice", fun () -> parse [ "--only"; "hires"; "--only"; "robertson" ]);
    ("an unknown flag", fun () -> parse [ "--fast" ]);
  ]

(* The suite and its pins: every case has a reference of its own length, ending where the problem does, and every row
   the bench runs has a pinned row, so that a case added without a pin fails here and not silently in the gate. *)
let suite =
  let consistent (case : Suite.case) =
    case.reference.t_end = case.problem.t_end && Array.length case.reference.y = Array.length case.problem.y0
  in
  let pinned_row (case : Suite.case) tol =
    List.exists (fun (p : Measure.t) -> p.problem = case.name && p.tol = tol) Golden.rows
  in
  let all_pinned (case : Suite.case) = List.for_all (pinned_row case) case.tols in
  [
    ("every case has a reference of its length", fun () -> string_of_bool (List.for_all consistent Suite.cases));
    ("every row is pinned", fun () -> string_of_bool (List.for_all all_pinned Suite.cases));
  ]

(* The problems of bench/: (1, 3) is the Brusselator's rest state, where each cell sees the boundary values and the
   right-hand side is exactly 0; HIRES at its start has only two nonzero derivatives, -1.71 + 0.0007 and 1.71. *)
let problems =
  let rest = Array.init 80 (fun k -> if k mod 2 = 0 then 1. else 3.) in
  let largest = Array.fold_left (fun m v -> Float.max m (Float.abs v)) 0. in
  let hires = Array.map (Printf.sprintf "%.4f") (Hires.rhs 0. Hires.y0) in
  [
    ( "brusselator right-hand side at (1, 3)",
      fun () -> Printf.sprintf "max |f| = %g" (largest (Brusselator.rhs 0. rest)) );
    ("brusselator state has 80 components", fun () -> string_of_int (Array.length Brusselator.y0));
    ("hires right-hand side at y0", fun () -> String.concat "; " (Array.to_list hires));
  ]

(* The export for the transcription check: a hundred points per problem from a pure generator, the same on every
   platform; the digits shown are few, because the last bits of y depend on the platform's rounding. *)
let transcription =
  let lines = Transcription.lines () in
  let count (case : Suite.case) =
    let own line = String.starts_with ~prefix:(case.name ^ ",") line in
    Printf.sprintf "%s %d" case.name (List.length (List.filter own lines))
  in
  let digits field = match float_of_string_opt field with Some v -> Printf.sprintf "%.6g" v | None -> field in
  let first_five line = List.filteri (fun i _ -> i < 5) (String.split_on_char ',' line) in
  [
    ("points per problem", fun () -> String.concat ", " (List.map count Suite.cases));
    ( "the first point: problem, t, y",
      fun () ->
        match lines with line :: _ -> String.concat ", " (List.map digits (first_five line)) | [] -> "no points" );
  ]

let () = Report.lines (List.concat [ accuracy; gate; pin; csv; options; suite; problems; transcription ])
