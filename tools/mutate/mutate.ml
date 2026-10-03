(* The command line: it reads files and prints, and leaves the rest to the pure modules. *)

let usage = "usage: mutate list FILE..."

let list file =
  match Mutant.enumerate ~file (Files.read file) with
  | Ok mutants -> List.iter (fun m -> print_endline (Mutant.to_string m)) mutants
  | Error message ->
      prerr_endline message;
      exit 1

let () =
  match List.tl (Array.to_list Sys.argv) with
  | "list" :: (_ :: _ as files) -> List.iter list files
  | _ ->
      prerr_endline usage;
      exit 2
