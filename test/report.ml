(* The only module that prints. A case computes its text and returns it, so
   cases stay pure; every line is name: text, so a diff against the .expected
   file names the case that changed. Effect rules: docs/architecture.md. *)
let lines cases = List.iter (fun (name, run) -> Printf.printf "%s: %s\n" name (run ())) cases
