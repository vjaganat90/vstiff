(** Perturbed spellings of numeric literals, for the operators that edit constants. *)

(** [floats text] pairs an operator name with the perturbed spelling of the float literal [text]: [float_x2],
    [float_div2], [float_zero], [float_x10] and [float_div10]. Ten times and a tenth move the decimal exponent, so
    [1e-6] gives [1e-5]. A perturbation that leaves the value alone or overflows is left out. Each spelling reads back
    as the value it stands for. *)
val floats : string -> (string * string) list

(** [ints text] is [int_plus1] and [int_minus1] of the integer literal [text], written in decimal. *)
val ints : string -> (string * string) list
