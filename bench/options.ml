type mode = Table | Csv | Pin | Help
type t = { mode : mode; cpu : bool; only : string option }

let usage =
  "usage: bench [--cpu] [--only NAME] [--csv | --pin | --help]\n\
  \  (no flag)      run every row and print the table, with a verdict per row against the golden table;\n\
  \                 the exit status is 1 unless every verdict is ok\n\
  \  --csv          print the rows as CSV instead, without verdicts\n\
  \  --pin          print the source of bench/golden.ml for this run (bench/README.md)\n\
  \  --cpu          add the CPU time of each row; it is not reproducible, so nothing compares it\n\
  \  --only NAME    run only the rows of the problem NAME, as the first column names it\n\
  \  --help         print this text\n"

let modes = [ ("--csv", Csv); ("--pin", Pin); ("--help", Help) ]

(* The modes that cover every problem: a partial pin would drop the other rows. *)
let whole = function Pin -> true | Table | Csv | Help -> false

(* Two mode flags, or two problems, are a mistake: the later one would silently win. The mode keeps its flag for the
   message. *)
let parse args =
  let rec go mode cpu only = function
    | [] -> (
        match mode with
        | Some (flag, mode) when whole mode && (cpu || only <> None) ->
            Error (Printf.sprintf "%s covers every problem and takes neither --cpu nor --only" flag)
        | Some (_, mode) -> Ok { mode; cpu; only }
        | None -> Ok { mode = Table; cpu; only })
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
