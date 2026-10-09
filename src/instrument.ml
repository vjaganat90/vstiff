(* Why only this module is mutable: docs/architecture.md. A ref cell is a mutable box: docs/ocaml.md. *)
let count (rhs : Ode.rhs) =
  let calls = ref 0 in
  ( (fun t y ->
      incr calls;
      rhs t y),
    fun () -> !calls )
