let count (rhs : Ode.rhs) =
  let calls = ref 0 in
  ( (fun t y ->
      incr calls;
      rhs t y),
    fun () -> !calls )
