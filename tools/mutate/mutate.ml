(* The command line: it reads files and prints, and leaves the rest to the pure modules. *)

let usage = "usage: mutate list FILE...\n       mutate apply ID"

let fail message =
  prerr_endline message;
  exit 1

let list file =
  match Mutant.enumerate ~file (Files.read file) with
  | Ok mutants -> List.iter (fun m -> print_endline (Mutant.to_string m)) mutants
  | Error message -> fail message

let apply id =
  match Mutant.file_of_id id with
  | None -> fail ("not a mutant id: " ^ id)
  | Some file -> (
      match Mutant.apply ~file (Files.read file) ~id with
      | Ok source -> print_string source
      | Error message -> fail message)

let () =
  match List.tl (Array.to_list Sys.argv) with
  | "list" :: (_ :: _ as files) -> List.iter list files
  | [ "apply"; id ] -> apply id
  | _ ->
      prerr_endline usage;
      exit 2
