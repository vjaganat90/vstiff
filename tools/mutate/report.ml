type run = { workers : int; baseline : float; timeout : float; wall : float }

let context = 24
let body = 72
let head n s = if String.length s <= n then s else String.sub s 0 n ^ "..."

let tail n s =
  let l = String.length s in
  if l <= n then s else "..." ^ String.sub s (l - n) n

(* The common prefix and suffix are cut to some context and the differing middle to some length. *)
let window a b =
  let la = String.length a and lb = String.length b in
  let limit = min la lb in
  let rec same_prefix i = if i < limit && a.[i] = b.[i] then same_prefix (i + 1) else i in
  let p = same_prefix 0 in
  let rec same_suffix i = if i < limit - p && a.[la - 1 - i] = b.[lb - 1 - i] then same_suffix (i + 1) else i in
  let s = same_suffix 0 in
  let show str =
    let l = String.length str in
    tail context (String.sub str 0 p)
    ^ head body (String.sub str p (l - p - s))
    ^ head context (String.sub str (l - s) s)
  in
  (show a, show b)

(* A code span; a backtick in the text needs a longer fence, with a space inside it. *)
let code s =
  let longest, _ =
    String.fold_left
      (fun (best, run) c ->
        let run = if c = '`' then run + 1 else 0 in
        (max best run, run))
      (0, 0) s
  in
  if longest = 0 then "`" ^ s ^ "`"
  else
    let fence = String.make (longest + 1) '`' in
    fence ^ " " ^ s ^ " " ^ fence

let duration seconds =
  let s = int_of_float (Float.round seconds) in
  if s < 60 then Printf.sprintf "%d s" s
  else if s < 3600 then Printf.sprintf "%d min %02d s" (s / 60) (s mod 60)
  else Printf.sprintf "%d h %02d min" (s / 3600) (s mod 3600 / 60)

let count pred entries = List.length (List.filter (fun (_, o) -> pred o) entries)
let killed = function Outcome.Killed _ -> true | _ -> false
let survived = function Outcome.Survived -> true | _ -> false
let stillborn = function Outcome.Stillborn _ -> true | _ -> false
let equivalent = function Outcome.Equivalent _ -> true | _ -> false
let timed_out o = o = Outcome.timeout
let row cells = "| " ^ String.concat " | " cells ^ " |\n"

let render ?run entries =
  let position ((m : Mutant.t), _) = (m.file, m.line, m.col, m.operator) in
  let entries = List.sort (fun a b -> compare (position a) (position b)) entries in
  let b = Buffer.create 1024 in
  let add = Buffer.add_string b in
  add "# Mutation report\n\n";
  Option.iter
    (fun r ->
      add
        (Printf.sprintf
           "%d workers, baseline %.1f s, timeout %.1f s, wall time %s. Mutants are built in dune's release profile, so \
            a warning does not make one stillborn.\n\n"
           r.workers r.baseline r.timeout (duration r.wall)))
    run;
  let files = List.sort_uniq compare (List.map (fun ((m : Mutant.t), _) -> m.file) entries) in
  let with_equivalent = count equivalent entries > 0 in
  let cells name es =
    [ name; string_of_int (List.length es) ]
    @ List.map (fun pred -> string_of_int (count pred es)) [ killed; survived; stillborn ]
    @ if with_equivalent then [ string_of_int (count equivalent es) ] else []
  in
  add
    (row
       ([ "Module"; "Mutants"; "Killed"; "Survived"; "Stillborn" ] @ if with_equivalent then [ "Equivalent" ] else []));
  add (row ([ "---"; "---:"; "---:"; "---:"; "---:" ] @ if with_equivalent then [ "---:" ] else []));
  List.iter (fun file -> add (row (cells file (List.filter (fun ((m : Mutant.t), _) -> m.file = file) entries)))) files;
  add (row (cells "**Total**" entries));
  let timeouts = count timed_out entries in
  if timeouts > 0 then
    add (Printf.sprintf "\nOf the %d killed mutants, %d were killed by the timeout.\n" (count killed entries) timeouts);
  let survivors = List.filter (fun (_, o) -> survived o) entries in
  add "\n## Survivors\n";
  if survivors = [] then add "\nNo mutant survived.\n"
  else
    List.iter
      (fun file ->
        match List.filter (fun ((m : Mutant.t), _) -> m.file = file) survivors with
        | [] -> ()
        | here ->
            add (Printf.sprintf "\n### %s\n\n" file);
            List.iter
              (fun ((m : Mutant.t), _) ->
                let before, after = window m.original m.mutated in
                add (Printf.sprintf "- %s: %s -> %s\n" (code (Mutant.id m)) (code before) (code after)))
              here)
      files;
  let known = List.filter (fun (_, o) -> equivalent o) entries in
  if known <> [] then (
    add "\n## Equivalent\n\n";
    List.iter
      (fun (m, o) ->
        match o with
        | Outcome.Equivalent "" -> add (Printf.sprintf "- %s\n" (code (Mutant.id m)))
        | Outcome.Equivalent reason -> add (Printf.sprintf "- %s: %s\n" (code (Mutant.id m)) reason)
        | _ -> ())
      known);
  Buffer.contents b
