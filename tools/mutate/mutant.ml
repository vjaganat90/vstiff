(* A mutation is a [site]: a node of the parse tree and a copy of it with one change. Every operator edits a single
   node, so one function finds the sites of a node and can serve both to list them and to pick one. *)
open Parsetree

type t = { file : string; line : int; col : int; operator : string; original : string; mutated : string }

let format_id ~file ~line ~col ~operator = Printf.sprintf "%s:%d:%d:%s" file line col operator
let id m = format_id ~file:m.file ~line:m.line ~col:m.col ~operator:m.operator
let to_string m = String.concat "\t" [ id m; m.original; m.mutated ]

(* The operator and the position have no colon, so the file is what is left once the last three fields are gone. *)
let file_of_id id =
  match List.rev (String.split_on_char ':' id) with
  | _operator :: col :: line :: (_ :: _ as file) when int_of_string_opt line <> None && int_of_string_opt col <> None ->
      Some (String.concat ":" (List.rev file))
  | _ -> None

(* [edited] replaces the node, [before] and [after] are the smaller expressions a report shows. The anchor gives the
   position of the id: the operator for a swap and for the operands of a connective, the keyword of an if or a while,
   the start of a literal, a guard or an arm. Two sites of one position differ in their operator name. *)
type site = { anchor : Location.t; operator : string; edited : expression; before : expression; after : expression }

(* Name used in ids, and replacements. A comparison becomes its boundary neighbour and its negation. *)
let swaps =
  [
    ("+.", ("fadd", [ "-." ]));
    ("-.", ("fsub", [ "+." ]));
    ("*.", ("fmul", [ "/." ]));
    ("/.", ("fdiv", [ "*." ]));
    ("+", ("add", [ "-" ]));
    ("-", ("sub", [ "+" ]));
    ("*", ("mul", [ "/" ]));
    ("/", ("div", [ "*" ]));
    ("<", ("lt", [ "<="; ">=" ]));
    ("<=", ("le", [ "<"; ">" ]));
    (">", ("gt", [ ">="; "<=" ]));
    (">=", ("ge", [ ">"; "<" ]));
    ("=", ("eq", [ "<>" ]));
    ("<>", ("ne", [ "=" ]));
    ("&&", ("and", [ "||" ]));
    ("||", ("or", [ "&&" ]));
  ]

let lident e = match e.pexp_desc with Pexp_ident { txt = Longident.Lident s; _ } -> Some s | _ -> None
let ident ~loc name = Ast_helper.Exp.ident ~loc (Location.mkloc (Longident.Lident name) loc)
let negate e = Ast_helper.Exp.apply ~loc:e.pexp_loc (ident ~loc:e.pexp_loc "not") [ (Asttypes.Nolabel, e) ]

(* [not] put on a condition, or taken off one that has it, as (operator name, new condition): double negation would
   only restate the original. *)
let toggle e =
  match e.pexp_desc with
  | Pexp_apply (f, [ (Asttypes.Nolabel, x) ]) when lident f = Some "not" -> ("not_remove", x)
  | _ -> ("not_add", negate e)

(* [head] is the operator, [whole] the node that holds it, and [rebuild] puts a new operator into that node. *)
let swap_sites ~whole head rebuild =
  match Option.bind (lident head) (fun op -> List.assoc_opt op swaps) with
  | None -> []
  | Some (name, targets) ->
      List.map
        (fun target ->
          let edited = rebuild (ident ~loc:head.pexp_loc target) in
          {
            anchor = head.pexp_loc;
            operator = name ^ "_to_" ^ fst (List.assoc target swaps);
            edited;
            before = whole;
            after = edited;
          })
        targets

(* The operands of && and || are conditions too. They are anchored at the connective, because the left operand
   starts where the whole expression does. *)
let operand_sites whole head = function
  | [ (l1, a); (l2, b) ] ->
      let site suffix operand rebuild =
        let name, operand' = toggle operand in
        {
          anchor = head.pexp_loc;
          operator = name ^ suffix;
          edited = rebuild operand';
          before = operand;
          after = operand';
        }
      in
      [
        site "_lhs" a (fun a' -> { whole with pexp_desc = Pexp_apply (head, [ (l1, a'); (l2, b) ]) });
        site "_rhs" b (fun b' -> { whole with pexp_desc = Pexp_apply (head, [ (l1, a); (l2, b') ]) });
      ]
  | _ -> []

let condition_site ~anchor cond rebuild =
  let name, cond' = toggle cond in
  { anchor; operator = name; edited = rebuild cond'; before = cond; after = cond' }

(* A match, a try and a function share their arms. An arm's guard is a condition; its body can be another arm's. *)
let arm_sites cases rebuild =
  let replace i case = List.mapi (fun j c -> if i = j then case else c) cases in
  let guard i case g =
    let name, g' = toggle g in
    {
      anchor = g.pexp_loc;
      operator = "guard_" ^ name;
      edited = rebuild (replace i { case with pc_guard = Some g' });
      before = g;
      after = g';
    }
  in
  let body i case j donor =
    {
      anchor = case.pc_lhs.ppat_loc;
      operator = Printf.sprintf "arm_body_from_%d" (j + 1);
      edited = rebuild (replace i { case with pc_rhs = donor.pc_rhs });
      before = case.pc_rhs;
      after = donor.pc_rhs;
    }
  in
  List.concat
    (List.mapi
       (fun i case ->
         Option.to_list (Option.map (guard i case) case.pc_guard)
         @ List.concat (List.mapi (fun j donor -> if i = j then [] else [ body i case j donor ]) cases))
       cases)

let literal_sites e variants text spell =
  List.map
    (fun (operator, text') ->
      let edited = { e with pexp_desc = Pexp_constant (spell text') } in
      { anchor = e.pexp_loc; operator; edited; before = e; after = edited })
    (variants text)

let sites e =
  let node desc = { e with pexp_desc = desc } in
  match e.pexp_desc with
  | Pexp_ident _ -> swap_sites ~whole:e e Fun.id
  | Pexp_apply (head, args) ->
      swap_sites ~whole:e head (fun head' -> node (Pexp_apply (head', args)))
      @ if lident head = Some "&&" || lident head = Some "||" then operand_sites e head args else []
  | Pexp_ifthenelse (cond, yes, no) ->
      let swapped =
        match no with
        | None -> []
        | Some no ->
            let edited = node (Pexp_ifthenelse (cond, no, Some yes)) in
            [ { anchor = e.pexp_loc; operator = "if_swap"; edited; before = e; after = edited } ]
      in
      condition_site ~anchor:e.pexp_loc cond (fun cond' -> node (Pexp_ifthenelse (cond', yes, no))) :: swapped
  | Pexp_while (cond, body) -> [ condition_site ~anchor:e.pexp_loc cond (fun cond' -> node (Pexp_while (cond', body))) ]
  | Pexp_constant { pconst_desc = Pconst_float (text, None); _ } ->
      literal_sites e Literal.floats text (fun s -> Ast_helper.Const.float s)
  | Pexp_constant { pconst_desc = Pconst_integer (text, None); _ } ->
      literal_sites e Literal.ints text (fun s -> Ast_helper.Const.integer s)
  | Pexp_match (scrutinee, cases) -> arm_sites cases (fun cases' -> node (Pexp_match (scrutinee, cases')))
  | Pexp_try (body, cases) -> arm_sites cases (fun cases' -> node (Pexp_try (body, cases')))
  | Pexp_function (params, constraint_, Pfunction_cases (cases, loc, attrs)) ->
      arm_sites cases (fun cases' -> node (Pexp_function (params, constraint_, Pfunction_cases (cases', loc, attrs))))
  | _ -> []

(* On one line: Pprintast breaks at its margin, which would split a table cell. *)
let show e =
  let buffer = Buffer.create 64 in
  let ppf = Format.formatter_of_buffer buffer in
  Format.pp_set_margin ppf 1_000_000;
  Pprintast.expression ppf e;
  Format.pp_print_flush ppf ();
  String.map (function '\n' | '\t' -> ' ' | c -> c) (Buffer.contents buffer)
  |> String.split_on_char ' '
  |> List.filter (( <> ) "")
  |> String.concat " "

let position (loc : Location.t) =
  let p = loc.loc_start in
  (p.pos_lnum, p.pos_cnum - p.pos_bol + 1)

(* Arms with equal bodies make equal programs: keep the first of each. *)
let distinct ms =
  let seen = Hashtbl.create 64 in
  List.filter
    (fun m ->
      let key = (m.line, m.col, m.mutated) in
      if Hashtbl.mem seen key then false
      else (
        Hashtbl.add seen key ();
        true))
    ms

let mutants ~file sites =
  List.filter_map
    (fun s ->
      let original = show s.before and mutated = show s.after in
      let line, col = position s.anchor in
      if original = mutated then None else Some { file; line; col; operator = s.operator; original; mutated })
    sites
  |> List.sort (fun a b -> compare (a.line, a.col, a.operator) (b.line, b.col, b.operator))
  |> distinct

let parse ~file source =
  let lexbuf = Lexing.from_string source in
  Location.init lexbuf file;
  match Parse.implementation lexbuf with
  | structure -> Ok structure
  | exception exn -> (
      match Location.error_of_exn exn with
      | Some (`Ok report) -> Error (String.trim (Format.asprintf "%a" Location.print_report report))
      | _ -> Error (Printexc.to_string exn))

(* The head of an application is its operator, found with the application itself: visiting it as well would list
   every swap twice. *)
let enumerate ~file source =
  Result.map
    (fun structure ->
      let found = ref [] in
      let expr (self : Ast_iterator.iterator) e =
        found := sites e @ !found;
        match e.pexp_desc with
        | Pexp_apply ({ pexp_desc = Pexp_ident _; _ }, args) -> List.iter (fun (_, a) -> self.expr self a) args
        | _ -> Ast_iterator.default_iterator.expr self e
      in
      let iterator = { Ast_iterator.default_iterator with expr } in
      iterator.structure iterator structure;
      mutants ~file !found)
    (parse ~file source)

(* A trailing newline, which Pprintast leaves out. *)
let print structure = Pprintast.string_of_structure structure ^ "\n"
let reprint ~file source = Result.map print (parse ~file source)

(* The mapper stops at the node it replaces, so the rest of the module is the original's, and it skips the head of
   an application for the same reason as [enumerate]. *)
let apply ~file source ~id:wanted =
  Result.bind (parse ~file source) (fun structure ->
      let found = ref false in
      let expr (self : Ast_mapper.mapper) e =
        let wanted_here s =
          let line, col = position s.anchor in
          format_id ~file ~line ~col ~operator:s.operator = wanted
        in
        match List.find_opt wanted_here (sites e) with
        | Some s ->
            found := true;
            s.edited
        | None -> (
            match e.pexp_desc with
            | Pexp_apply (({ pexp_desc = Pexp_ident _; _ } as head), args) ->
                { e with pexp_desc = Pexp_apply (head, List.map (fun (label, a) -> (label, self.expr self a)) args) }
            | _ -> Ast_mapper.default_mapper.expr self e)
      in
      let mapper = { Ast_mapper.default_mapper with expr } in
      let mutated = mapper.structure mapper structure in
      if !found then Ok (print mutated) else Error ("no mutant " ^ wanted))
