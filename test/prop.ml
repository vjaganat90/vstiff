(* Each property draws from its own state, made from the run seed and its name, so adding a property changes no
   other property's cases. Every evaluation goes through Guard.verdict, the only place that catches. *)

module type Property = sig
  type t

  val name : string
  val gen : t Gen.t
  val show : t -> string
  val holds : t -> bool
end

module type Statistical = sig
  include Property

  val threshold : float
end

type outcome =
  | Passed of int
  | Met of { threshold : float; count : int }
  | Failed of { seed : int; case : string; shrinks : int; why : string }

(* Candidates evaluated while shrinking, at most: a node of a matrix's tree offers thousands. *)
let max_tries = 10_000
let fails holds x = match Guard.verdict holds x with Ok true -> false | Ok false | Error _ -> true

(* The first candidate that still fails, and the tries left after it. *)
let rec first fails tries candidates =
  if tries = 0 then (None, 0)
  else
    match candidates () with
    | Seq.Nil -> (None, tries)
    | Seq.Cons ((Gen.Node (x, _) as t), rest) -> if fails x then (Some t, tries - 1) else first fails (tries - 1) rest

(* Greedy descent: move to the first smaller case that still fails, until none does or the tries run out. *)
let shrink fails tree =
  let rec descend (Gen.Node (x, smaller)) steps tries =
    match first fails tries smaller with Some t, tries -> descend t (steps + 1) tries | None, _ -> (x, steps)
  in
  descend tree 0 max_tries

(* A modular explicit: the type of tree names P.t (docs/ocaml.md). *)
let failed (module P : Property) ~seed (tree : P.t Gen.tree) why =
  let x, shrinks = shrink (fails P.holds) tree in
  let why = match Guard.verdict P.holds x with Error e -> why ^ ", raises " ^ e | Ok _ -> why in
  Failed { seed; case = P.show x; shrinks; why }

let stream seed name = Random.State.make [| seed; Hashtbl.hash name |]

let check (module P : Property) ~count ~seed =
  let st = stream seed P.name in
  let rec run k =
    if k = count then Passed count
    else
      let tree = P.gen st in
      if fails P.holds (Gen.root tree) then
        failed (module P) ~seed tree (Printf.sprintf "at case %d of %d" (k + 1) count)
      else run (k + 1)
  in
  run 0

let check_statistical (module P : Statistical) ~count ~seed =
  let st = stream seed P.name in
  (* Every case runs, since the share needs them all; the first failing one is kept to show. *)
  let rec run k held first =
    if k = count then (held, first)
    else
      let tree = P.gen st in
      if fails P.holds (Gen.root tree) then run (k + 1) held (if Option.is_none first then Some tree else first)
      else run (k + 1) (held + 1) first
  in
  let held, first = run 0 0 None in
  if Float.of_int held >= P.threshold *. Float.of_int count then Met { threshold = P.threshold; count }
  else
    let why = Printf.sprintf "%d of %d held, fewer than %g%%" held count (100. *. P.threshold) in
    match first with
    | Some tree -> failed (module P) ~seed tree why
    | None -> Failed { seed; case = "none failed"; shrinks = 0; why }

let to_string = function
  | Passed count -> Printf.sprintf "ok %d" count
  | Met { threshold; count } -> Printf.sprintf "ok, at least %g%% of %d" (100. *. threshold) count
  | Failed { seed; case; shrinks; why } ->
      Printf.sprintf "FAILED %s, seed %d, shrunk in %d steps: %s" why seed shrinks case
