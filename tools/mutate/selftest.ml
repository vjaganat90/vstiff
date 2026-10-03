(* The tool's own test: its output, in selftest.expected, is the mutant list of the fixture. docs/testing.md *)

let () =
  let file = Sys.argv.(1) in
  match Mutant.enumerate ~file (Files.read file) with
  | Ok mutants -> List.iter (fun m -> print_endline (Mutant.to_string m)) mutants
  | Error message -> print_endline message
