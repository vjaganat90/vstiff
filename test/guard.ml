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

(* bounded run is Some of what run returns, or None when a budget inside it ran
   out: the soak's way to stop a round without catching anything itself. *)
let bounded run = match run () with v -> Some v | exception Exhausted -> None

(* verdict holds x is Ok of holds x, or Error with the exception it raised: a
   property case that raises fails like one that does not hold, and the run goes
   on (Prop). *)
let verdict holds x = match holds x with b -> Ok b | exception e -> Error (Printexc.to_string e)

(* run ok case is a table entry's function of (): ok formats an Ok; an Error, an
   exhausted budget and an Invalid_argument each print one line. A match case
   written exception catches what its expression raises (docs/ocaml.md). *)
let run ok case () =
  match case () with
  | Ok v -> ok v
  | Error e -> "Error " ^ Fail.to_string e
  | exception Exhausted -> "no answer within 5e6 rhs calls"
  | exception Invalid_argument m -> "Invalid_argument " ^ m
