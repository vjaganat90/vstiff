(* The tool's own test; selftest.expected pins its output: the mutants of the fixture, a check that every one of them
   applies, and a few small applications printed in full. docs/testing.md *)

let print_result = function Ok text -> print_string text | Error message -> print_endline ("error: " ^ message)

(* The fixture's mutants, each one applied: the result must parse, and differ from the module printed unchanged. *)
let fixture file =
  let source = Files.read file in
  match (Mutant.enumerate ~file source, Mutant.reprint ~file source) with
  | Ok mutants, Ok unchanged ->
      List.iter (fun m -> print_endline (Mutant.to_string m)) mutants;
      let applied = List.map (fun m -> Mutant.apply ~file source ~id:(Mutant.id m)) mutants in
      let sound = function
        | Ok text -> text <> unchanged && Result.is_ok (Mutant.enumerate ~file text)
        | Error _ -> false
      in
      let ids = List.sort_uniq compare (List.map Mutant.id mutants) in
      Printf.printf "%d mutants, ids distinct: %b, each applies to a different module that parses: %b\n"
        (List.length mutants)
        (List.length ids = List.length mutants)
        (List.for_all sound applied);
      (* Printing alone must keep every mutant: the same operators on the same expressions, at other positions. *)
      let shape = List.map (fun (m : Mutant.t) -> (m.operator, m.original, m.mutated)) in
      Printf.printf "printing the fixture unchanged keeps its mutants: %b\n"
        (Result.map shape (Mutant.enumerate ~file unchanged) = Ok (shape mutants))
  | Error message, _ | _, Error message -> print_endline message

let small = "let f x = if x > 0 then x else -x\nlet g = function Some y -> y | None -> 0\n"

let () =
  fixture Sys.argv.(1);
  List.iter
    (fun id ->
      print_endline id;
      print_result (Mutant.apply ~file:"small.ml" small ~id))
    [ "small.ml:1:11:if_swap"; "small.ml:1:16:gt_to_le"; "small.ml:2:18:arm_body_from_2"; "small.ml:2:18:if_swap" ];
  Printf.printf "a syntax error is an Error: %b\n" (Result.is_error (Mutant.enumerate ~file:"bad.ml" "let f = (\n"));
  Printf.printf "file_of_id: %s %s\n"
    (Option.value (Mutant.file_of_id "src/a:b.ml:3:14:lt_to_le") ~default:"none")
    (Option.value (Mutant.file_of_id "src/a.ml") ~default:"none")
