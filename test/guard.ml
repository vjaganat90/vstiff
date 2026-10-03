(* Runs cases so that whatever escapes them becomes a line of text. The only
   test module that raises or catches. *)
open Vstiff

exception Exhausted

let limit = 5_000_000

(* [rhs] on a budget: past [limit] calls it raises instead of evaluating, so a
   loop that never ends stops the case instead of hanging the run. *)
let budget rhs =
  let counted, calls = Instrument.count rhs in
  fun t y -> if calls () >= limit then raise Exhausted else counted t y

let run ok case () =
  match case () with
  | Ok v -> ok v
  | Error e -> "Error " ^ Fail.to_string e
  | exception Exhausted -> "no answer within 5e6 rhs calls"
  | exception Invalid_argument m -> "Invalid_argument " ^ m
