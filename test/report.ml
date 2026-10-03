(* The only module that prints: one line per case, its name and its text. *)
let lines cases = List.iter (fun (name, run) -> Printf.printf "%s: %s\n" name (run ())) cases
