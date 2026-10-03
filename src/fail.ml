(* Failures are values: docs/architecture.md. Variants, Result and binding operators: docs/ocaml.md. *)

type t = Diverged | StepRejected of int | Nan

let to_string = function
  | Diverged -> "Diverged"
  | StepRejected n -> Printf.sprintf "StepRejected %d" n
  | Nan -> "Nan"

module Syntax = struct
  let ( let* ) = Result.bind
  let ( let+ ) r f = Result.map f r
end
