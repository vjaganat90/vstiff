(* The tool's own test; selftest.expected pins its output. First the mutants of the fixture, each applied, and a few
   small applications printed in full; then the results file, the build log summaries, the equivalent list, the
   window and the report, on outcomes made up for the fixture's mutants. docs/testing.md *)

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
        (Result.map shape (Mutant.enumerate ~file unchanged) = Ok (shape mutants));
      mutants
  | Error message, _ | _, Error message ->
      print_endline message;
      []

let small = "let f x = if x > 0 then x else -x\nlet g = function Some y -> y | None -> 0\n"

let applications () =
  List.iter
    (fun id ->
      print_endline id;
      print_result (Mutant.apply ~file:"small.ml" small ~id))
    [ "small.ml:1:11:if_swap"; "small.ml:1:16:gt_to_le"; "small.ml:2:18:arm_body_from_2"; "small.ml:2:18:if_swap" ];
  Printf.printf "a syntax error is an Error: %b\n" (Result.is_error (Mutant.enumerate ~file:"bad.ml" "let f = (\n"));
  Printf.printf "file_of_id: %s %s\n"
    (Option.value (Mutant.file_of_id "src/a:b.ml:3:14:lt_to_le") ~default:"none")
    (Option.value (Mutant.file_of_id "src/a.ml") ~default:"none")

let literals () =
  let show = function
    | [] -> "none"
    | variants -> String.concat " " (List.map (fun (operator, text) -> operator ^ "=" ^ text) variants)
  in
  print_endline "\n== float literals ==";
  List.iter
    (fun text -> Printf.printf "%s: %s\n" text (show (Literal.floats text)))
    [ "1e-10"; "2."; "0.5"; "1_000."; "1e-4"; "16."; "123456."; "0."; "-1."; "1e308"; "0x1p3"; "1e-320" ];
  print_endline "\n== int literals ==";
  List.iter (fun text -> Printf.printf "%s: %s\n" text (show (Literal.ints text))) [ "0"; "0x10"; "1_000"; "-5" ]

(* One of each outcome in turn, the files alternating, and one survivor whose text has backticks and one equivalent mutant with no reason. *)
let made_up (mutants : Mutant.t list) =
  let outcome i =
    match i mod 7 with
    | 0 -> Outcome.Survived
    | 1 -> Outcome.Stillborn "Error: Unbound value x"
    | 2 -> Outcome.Killed "timeout"
    | 3 -> Outcome.Equivalent "same result on every input"
    | _ -> Outcome.Killed "test/corpus.expected"
  in
  let entries =
    List.mapi
      (fun i (m : Mutant.t) -> ((if i mod 2 = 0 then m else { m with file = "src/other.ml" }), outcome i))
      mutants
  in
  match mutants with
  | m :: _ ->
      let poly = { m with file = "src/poly.ml" } in
      entries
      @ [
          ({ poly with original = "`A x"; mutated = "`B x" }, Outcome.Survived);
          ({ poly with operator = "int_plus1" }, Outcome.Equivalent "");
        ]
  | [] -> entries

let results (entries : (Mutant.t * Outcome.t) list) =
  print_endline "\n== results file ==";
  print_endline (Outcome.to_line (List.hd entries));
  let odd =
    match entries with
    | (m, _) :: _ -> [ ({ m with original = "\"a\\tb\"\n\\"; mutated = "tab\t" }, Outcome.Killed "x\ty") ]
    | [] -> []
  in
  Printf.printf "every line reads back as the entry it came from: %b\n"
    (List.for_all (fun e -> Outcome.of_line (Outcome.to_line e) = Some e) (odd @ entries));
  Printf.printf "another line is skipped: %b\n" (Outcome.of_line "not a results line" = None)

let logs () =
  print_endline "\n== build logs ==";
  let diff =
    "File \"test/corpus.expected\", line 1, characters 0-0:\n--- a/corpus.expected\n+++ b/corpus.exe.output\n"
  in
  let crash =
    "File \"test/dune\", line 5, characters 8-14:\n\
     5 |  (names corpus soak)\n\
    \            ^^^^^^\n\
     Fatal error: exception Invalid_argument(\"index out of bounds\")\n"
  in
  let error =
    "File \"src/numerics/newton.ml\", line 20, characters 7-8:\n20 |  x\n       ^\nError: Unbound value n\n"
  in
  List.iter
    (fun log -> Printf.printf "%s: %s\n" (Outcome.name (Outcome.killed ~log)) (Outcome.detail (Outcome.killed ~log)))
    [ diff; crash; "dune: no such command\n" ];
  List.iter
    (fun log ->
      Printf.printf "%s: %s\n" (Outcome.name (Outcome.stillborn ~log)) (Outcome.detail (Outcome.stillborn ~log)))
    [ error; "" ]

let equivalents () =
  print_endline "\n== equivalent list ==";
  List.iter
    (fun (id, reason) -> Printf.printf "%s -> %S\n" id reason)
    (Outcome.equivalents
       "# a comment\n\nfixture.ml:3:18:fmul_to_fdiv  w is 1 in every test\nfixture.ml:4:21:int_plus1\n")

let windows () =
  print_endline "\n== windows ==";
  let run c = String.make 40 c in
  List.iter
    (fun (a, b) ->
      let a, b = Report.window a b in
      Printf.printf "%s | %s\n" a b)
    [
      ("x +. y", "x -. y");
      (run 'a' ^ " +. " ^ run 'b', run 'a' ^ " -. " ^ run 'b');
      ("if c then " ^ String.make 100 'x' ^ " else y", "if c then y else " ^ String.make 100 'x');
    ]

let report entries =
  print_endline "\n== report ==";
  print_string (Report.render ~run:{ Report.workers = 4; baseline = 6.94; timeout = 20.82; wall = 852. } entries);
  print_endline "\n== report of nothing ==";
  print_string (Report.render [])

let () =
  let mutants = fixture Sys.argv.(1) in
  applications ();
  literals ();
  let entries = made_up mutants in
  results entries;
  logs ();
  equivalents ();
  windows ();
  report entries
