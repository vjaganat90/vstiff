(* The command line: it reads files and prints, and leaves the rest to the other modules. *)

let usage =
  {|usage: mutate list FILE...        list the mutants, one per line: id, original, mutated
       mutate apply ID          print the module with that mutant
       mutate report RESULTS    print the report of a results file|}

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

let read_results file = List.filter_map Outcome.of_line (String.split_on_char '\n' (Files.read file))

let report file =
  let latest = Hashtbl.create 64 in
  List.iter (fun (m, o) -> Hashtbl.replace latest (Mutant.id m) (m, o)) (read_results file);
  print_string (Report.render (List.of_seq (Hashtbl.to_seq_values latest)))

let main () =
  match List.tl (Array.to_list Sys.argv) with
  | "list" :: (_ :: _ as files) -> List.iter list files
  | [ "apply"; id ] -> apply id
  | [ "report"; file ] -> report file
  | _ ->
      prerr_endline usage;
      exit 2

(* A file that cannot be read is the user's mistake, not a crash. *)
let () = try main () with Sys_error message -> fail message
