(** Time as binary64 represents it. *)

(** [resolution t] is [16 eps |t|], 16 to 32 ulp of [t] ([eps] is
    [Float.epsilon]): about the shortest step that [t + h] still resolves to a
    few percent. Below it a step barely moves [t], or not at all. *)
val resolution : float -> float
