(* The only test module that raises or catches, so that a case stays pure: an
   exhausted call budget and an Invalid_argument become a line of text, and a
   runaway loop or a bad argument shows in the diff instead of hanging or
   aborting the run. Effect rules: docs/architecture.md. *)
open Vstiff

exception Exhausted

(* run prints this limit as 5e6: keep the two equal. *)
let limit = 5_000_000

(* rhs on a budget: past limit calls it raises instead of evaluating, so a loop
   that never ends stops the case instead of hanging the run. *)
let budget rhs =
  let counted, calls = Instrument.count rhs in
  fun t y -> if calls () >= limit then raise Exhausted else counted t y

(* run ok case is a table entry's function of (): ok formats an Ok; an Error, an
   exhausted budget and an Invalid_argument each print one line. A match case
   written exception catches what its expression raises (docs/ocaml.md). *)
let run ok case () =
  match case () with
  | Ok v -> ok v
  | Error e -> "Error " ^ Fail.to_string e
  | exception Exhausted -> "no answer within 5e6 rhs calls"
  | exception Invalid_argument m -> "Invalid_argument " ^ m
