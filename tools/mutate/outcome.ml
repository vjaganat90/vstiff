type t = Killed of string | Survived | Stillborn of string | Equivalent of string

let timeout = Killed "timeout"
let lines text = String.split_on_char '\n' text
let first prefix text = List.find_opt (String.starts_with ~prefix) (lines text)

let killed ~log =
  match first "Fatal error: " log with
  | Some line -> Killed (String.sub line 13 (String.length line - 13))
  | None -> (
      match first "File \"" log with
      | Some line -> ( match String.split_on_char '"' line with _ :: name :: _ -> Killed name | _ -> Killed "failed")
      | None -> Killed "failed")

let stillborn ~log = Stillborn (Option.value (first "Error" log) ~default:"does not compile")

let name = function
  | Killed _ -> "killed"
  | Survived -> "survived"
  | Stillborn _ -> "stillborn"
  | Equivalent _ -> "equivalent"

let detail = function Killed d | Stillborn d | Equivalent d -> d | Survived -> ""

let escape s =
  let b = Buffer.create (String.length s) in
  String.iter
    (function
      | '\\' -> Buffer.add_string b "\\\\"
      | '\t' -> Buffer.add_string b "\\t"
      | '\n' -> Buffer.add_string b "\\n"
      | c -> Buffer.add_char b c)
    s;
  Buffer.contents b

let unescape s =
  let b = Buffer.create (String.length s) in
  let rec go i =
    if i < String.length s then
      if s.[i] = '\\' && i + 1 < String.length s then (
        Buffer.add_char b (match s.[i + 1] with 't' -> '\t' | 'n' -> '\n' | c -> c);
        go (i + 2))
      else (
        Buffer.add_char b s.[i];
        go (i + 1))
  in
  go 0;
  Buffer.contents b

let to_line ((m : Mutant.t), outcome) =
  String.concat "\t"
    (List.map escape
       [
         m.file;
         string_of_int m.line;
         string_of_int m.col;
         m.operator;
         name outcome;
         detail outcome;
         m.original;
         m.mutated;
       ])

let of_line text =
  match List.map unescape (String.split_on_char '\t' text) with
  | [ file; line; col; operator; kind; detail; original; mutated ] -> (
      let outcome =
        match kind with
        | "killed" -> Some (Killed detail)
        | "survived" -> Some Survived
        | "stillborn" -> Some (Stillborn detail)
        | "equivalent" -> Some (Equivalent detail)
        | _ -> None
      in
      match (int_of_string_opt line, int_of_string_opt col, outcome) with
      | Some line, Some col, Some outcome -> Some ({ Mutant.file; line; col; operator; original; mutated }, outcome)
      | _ -> None)
  | _ -> None

let equivalents text =
  List.filter_map
    (fun line ->
      let line = String.trim line in
      if line = "" || line.[0] = '#' then None
      else
        match String.index_opt line ' ' with
        | None -> Some (line, "")
        | Some i -> Some (String.sub line 0 i, String.trim (String.sub line i (String.length line - i))))
    (lines text)
