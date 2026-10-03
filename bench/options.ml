type mode = Table | Csv | Help
type t = { mode : mode; cpu : bool; only : string option }

let usage =
  "usage: bench [--cpu] [--only NAME] [--csv | --help]\n\
  \  (no flag)      run every row and print the table\n\
  \  --csv          print the rows as CSV instead\n\
  \  --cpu          add the CPU time of each row; it is not reproducible, so nothing compares it\n\
  \  --only NAME    run only the rows of the problem NAME, as the first column names it\n\
  \  --help         print this text\n"

let modes = [ ("--csv", Csv); ("--help", Help) ]

(* Two mode flags, or two problems, are a mistake: the later one would silently win. The mode keeps its flag for the
   message. *)
let parse args =
  let rec go mode cpu only = function
    | [] -> Ok { mode = (match mode with Some (_, m) -> m | None -> Table); cpu; only }
    | "--cpu" :: rest -> go mode true only rest
    | "--only" :: name :: rest -> if only = None then go mode cpu (Some name) rest else Error "give --only once"
    | [ "--only" ] -> Error "--only needs a problem name"
    | flag :: rest -> (
        match (List.assoc_opt flag modes, mode) with
        | None, _ -> Error (Printf.sprintf "unknown argument %s" flag)
        | Some m, None -> go (Some (flag, m)) cpu only rest
        | Some _, Some _ -> Error (Printf.sprintf "give at most one of %s" (String.concat ", " (List.map fst modes))))
  in
  go None false None args
