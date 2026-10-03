(** Mutation sites of an OCaml source file, found with the compiler's own parser. Pure: sources come in as strings.

    The operator is the last part of an id:
    - [fadd_to_fsub] and the like: [+.] [-.] [*.] [/.] and [+] [-] [*] [/] swap in pairs, [&&] and [||] swap, [=] and
      [<>] swap, and a comparison becomes its boundary neighbour or its negation ([lt_to_le], [lt_to_ge]);
    - [not_add], [not_remove]: [not] added to, or taken off, the condition of an [if] or a [while]; [guard_not_add]
      and [guard_not_remove] do it to a [when] guard, [not_add_lhs], [not_remove_rhs] and so on to an operand of [&&]
      or [||];
    - [if_swap]: the two branches of an [if] trade places;
    - [float_x2], [float_div2], [float_zero], [float_x10], [float_div10] for a float literal, [int_plus1] and
      [int_minus1] for an integer literal;
    - [arm_body_from_k]: the body of a match arm is replaced by the body of arm [k], counting from 1. *)

(** One mutant. [line] and [col] count from 1 and locate the construct the operator edits; [original] and [mutated]
    print the expression that changes, on one line. *)
type t = { file : string; line : int; col : int; operator : string; original : string; mutated : string }

(** [file:line:col:operator]: unique within a file, and stable while the source is unchanged. *)
val id : t -> string

(** [id], [original] and [mutated], tab-separated. *)
val to_string : t -> string

(** [enumerate ~file source] lists the mutants of [source] in order of position, or the syntax error. [file] only
    names the mutants. A mutant whose print equals the original's is not listed, nor is a second one that makes the
    same change at the same position. *)
val enumerate : file:string -> string -> (t list, string) result
