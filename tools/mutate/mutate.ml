(* The command line: it reads files and prints, and leaves the rest to the other modules. *)

let usage =
  {|usage: mutate list FILE...        list the mutants, one per line: id, original, mutated
       mutate apply ID          print the module with that mutant
       mutate run [OPTIONS] FILE...
                                run the tests against every mutant, print the report
       mutate report RESULTS    print the report of a results file

options of run:
  -j, --workers N       scratch copies, and mutants run at once (default 4)
  --root DIR            the repository (default .); FILEs are relative to it
  --results FILE        append each outcome to FILE, and skip those already in it
  --equivalent FILE     mutants known to be equivalent: an id and a reason per line; not run
  --only TEXT           only the mutants whose id contains TEXT
  --timeout-factor X    a mutant may take X times the baseline run (default 3)
  --keep                leave the scratch copies in place|}

let fail message =
  prerr_endline message;
  exit 1

let enumerate ~file source =
  match Mutant.enumerate ~file source with Ok mutants -> mutants | Error message -> fail message

let list file = List.iter (fun m -> print_endline (Mutant.to_string m)) (enumerate ~file (Files.read file))

let apply id =
  match Mutant.file_of_id id with
  | None -> fail ("not a mutant id: " ^ id)
  | Some file -> (
      match Mutant.apply ~file (Files.read file) ~id with
      | Ok source -> print_string source
      | Error message -> fail message)

type options = {
  workers : int;
  root : string;
  results : string option;
  equivalent : string option;
  only : string option;
  timeout_factor : float;
  keep : bool;
  files : string list;
}

let defaults =
  {
    workers = 4;
    root = ".";
    results = None;
    equivalent = None;
    only = None;
    timeout_factor = 3.;
    keep = false;
    files = [];
  }

(* A path as the id carries it: ./src/x.ml and src/x.ml are one file. *)
let normalise file = if String.starts_with ~prefix:"./" file then String.sub file 2 (String.length file - 2) else file

let rec parse options = function
  | [] -> Ok { options with files = List.rev options.files }
  | ("-j" | "--workers") :: n :: rest -> (
      match int_of_string_opt n with
      | Some workers when workers > 0 -> parse { options with workers } rest
      | _ -> Error "--workers needs a positive integer")
  | "--root" :: root :: rest -> parse { options with root } rest
  | "--results" :: file :: rest -> parse { options with results = Some file } rest
  | "--equivalent" :: file :: rest -> parse { options with equivalent = Some file } rest
  | "--only" :: text :: rest -> parse { options with only = Some text } rest
  | "--timeout-factor" :: x :: rest -> (
      match float_of_string_opt x with
      | Some timeout_factor when timeout_factor > 0. -> parse { options with timeout_factor } rest
      | _ -> Error "--timeout-factor needs a positive number")
  | "--keep" :: rest -> parse { options with keep = true } rest
  | flag :: _ when String.starts_with ~prefix:"-" flag -> Error ("unknown option " ^ flag)
  | file :: rest -> parse { options with files = normalise file :: options.files } rest

let contains ~text s =
  let n = String.length text in
  let rec at i = i + n <= String.length s && (String.sub s i n = text || at (i + 1)) in
  at 0

let read_results file = List.filter_map Outcome.of_line (String.split_on_char '\n' (Files.read file))

let run options =
  if options.files = [] then fail usage;
  let selected m = match options.only with None -> true | Some text -> contains ~text (Mutant.id m) in
  let mutants =
    List.filter selected
      (List.concat_map (fun file -> enumerate ~file (Files.read (Filename.concat options.root file))) options.files)
  in
  let equivalent = match options.equivalent with None -> [] | Some file -> Outcome.equivalents (Files.read file) in
  let known = Hashtbl.create 64 in
  (match options.results with
  | Some file when Files.exists file ->
      List.iter (fun (m, o) -> Hashtbl.replace known (Mutant.id m) o) (read_results file)
  | _ -> ());
  let total = List.length mutants in
  let record m outcome =
    Hashtbl.replace known (Mutant.id m) outcome;
    Option.iter (fun file -> Files.append file ~contents:(Outcome.to_line (m, outcome) ^ "\n")) options.results
  in
  let progress m outcome =
    Printf.eprintf "[%*d/%d] %-9s %s  %s\n%!"
      (String.length (string_of_int total))
      (Hashtbl.length known) total (Outcome.name outcome) (Mutant.id m) (Outcome.detail outcome)
  in
  List.iter
    (fun m ->
      match (Hashtbl.find_opt known (Mutant.id m), List.assoc_opt (Mutant.id m) equivalent) with
      | None, Some reason -> record m (Outcome.Equivalent reason)
      | _ -> ())
    mutants;
  let todo = List.filter (fun m -> not (Hashtbl.mem known (Mutant.id m))) mutants in
  let entries () =
    List.filter_map (fun m -> Option.map (fun o -> (m, o)) (Hashtbl.find_opt known (Mutant.id m))) mutants
  in
  if todo = [] then print_string (Report.render (entries ()))
  else
    let workers = min options.workers (List.length todo) in
    let config =
      { Runner.root = options.root; workers; timeout_factor = options.timeout_factor; keep = options.keep }
    in
    match
      Runner.run config todo prerr_endline (fun m outcome ->
          record m outcome;
          progress m outcome)
    with
    | Ok timing ->
        let run = { Report.workers; baseline = timing.baseline; timeout = timing.timeout; wall = timing.wall } in
        print_string (Report.render ~run (entries ()))
    | Error message -> fail message

let report file =
  let latest = Hashtbl.create 64 in
  List.iter (fun (m, o) -> Hashtbl.replace latest (Mutant.id m) (m, o)) (read_results file);
  print_string (Report.render (List.of_seq (Hashtbl.to_seq_values latest)))

let main () =
  match List.tl (Array.to_list Sys.argv) with
  | "list" :: (_ :: _ as files) -> List.iter list files
  | [ "apply"; id ] -> apply id
  | "run" :: args -> ( match parse defaults args with Ok options -> run options | Error message -> fail message)
  | [ "report"; file ] -> report file
  | _ ->
      prerr_endline usage;
      exit 2

(* A file that cannot be read is the user's mistake, not a crash. *)
let () = try main () with Sys_error message -> fail message
