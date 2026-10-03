(* The parser drops plain comments, but keeps doc comments as attributes: on an
   item, or floating as an item of their own. Dropping those too leaves code. *)
let is_doc (a : Parsetree.attribute) =
  match a.attr_name.txt with "ocaml.doc" | "ocaml.text" -> true | _ -> false

let strip =
  let open Ast_mapper in
  let attributes m attrs = default_mapper.attributes m (List.filter (fun a -> not (is_doc a)) attrs) in
  let structure m items =
    default_mapper.structure m
      (List.filter
         (fun (i : Parsetree.structure_item) ->
           match i.pstr_desc with Pstr_attribute a -> not (is_doc a) | _ -> true)
         items)
  in
  let signature m items =
    default_mapper.signature m
      (List.filter
         (fun (i : Parsetree.signature_item) ->
           match i.psig_desc with Psig_attribute a -> not (is_doc a) | _ -> true)
         items)
  in
  { default_mapper with attributes; structure; signature }

(* Pprintast lays the tree out afresh, so layout does not count either. Its
   layout can change with the compiler release: after an upgrade every canon
   may need promoting again. *)
let () =
  let ast = Pparse.parse_implementation ~tool_name:"canon" Sys.argv.(1) in
  Format.printf "%a@." Pprintast.structure (strip.structure strip ast)
