(* The bench: vstiff's adaptive driver on the problems of Appendix B, one printed row per problem and tolerance. All the
   effects are Console's; this file only composes the pure modules of benchlib. bench/README.md says how to run it. *)
open Benchlib

let names = String.concat ", " (List.map (fun (case : Suite.case) -> case.name) Suite.cases)

let rows (options : Options.t) =
  let wanted (case : Suite.case) = match options.only with None -> true | Some name -> name = case.name in
  let with_tols (case : Suite.case) = List.map (fun tol -> (case, tol)) case.tols in
  List.concat_map with_tols (List.filter wanted Suite.cases)

let entry (options : Options.t) (case, tol) =
  if options.cpu then
    let row, seconds = Console.cpu (fun () -> Measure.run case ~tol) in
    { Render.row; cpu = Some seconds }
  else { Render.row = Measure.run case ~tol; cpu = None }

let failed { Render.row; _ } = match row.Measure.outcome with Measure.Failed _ -> true | Done _ -> false

(* The rows appear as they finish: the longest takes a few seconds. *)
let stream options line selected =
  List.map
    (fun row ->
      let entry = entry options row in
      Console.print (line entry);
      entry)
    selected

let table (options : Options.t) selected =
  Console.print Render.title;
  Console.print (Render.table_header options);
  if List.exists failed (stream options Render.table_line selected) then Console.exit 1

let csv (options : Options.t) selected =
  Console.print (Render.csv_header options);
  if List.exists failed (stream options Render.csv_line selected) then Console.exit 1

let with_rows options run =
  match rows options with
  | [] ->
      Console.warn ("no such problem; the problems are " ^ names);
      Console.exit 2
  | selected -> run selected

let () =
  match Options.parse (Console.args ()) with
  | Error message ->
      Console.warn message;
      Console.warn Options.usage;
      Console.exit 2
  | Ok options -> (
      match options.mode with
      | Help -> Console.print Options.usage
      | Table -> with_rows options (table options)
      | Csv -> with_rows options (csv options))
