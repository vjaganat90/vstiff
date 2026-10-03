(* The command line, the output, the clock and the exit status are effects. The bench keeps them in this one module
   so that every other module stays pure and a run can be rerun from its data. Effect rules: docs/architecture.md. *)

let args () = match Array.to_list Sys.argv with _ :: rest -> rest | [] -> []

(* print_newline flushes, so a long table appears row by row. *)
let print s =
  print_string s;
  print_newline ()

let warn s =
  prerr_string s;
  prerr_newline ()

let cpu f =
  let start = Sys.time () in
  let v = f () in
  (v, Sys.time () -. start)

let exit = Stdlib.exit
