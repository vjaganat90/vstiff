(* The shortest %g spelling that reads back as v, made into a float literal: it needs a dot or an exponent, and the
   exponent loses its plus sign and padding zeros, so 1e-9 stays 1e-9 and not 1e-09. Small whole numbers are written
   out, 20. and not 2e1; round large ones keep the exponent, as the sources do. *)
let spell v =
  let rec shortest digits =
    let s = Printf.sprintf "%.*g" digits v in
    if digits >= 17 || float_of_string s = v then s else shortest (digits + 1)
  in
  if Float.is_integer v && Float.abs v < 1e5 then Printf.sprintf "%.0f." v
  else
    let s = shortest 1 in
    match String.index_opt s 'e' with
    | Some i ->
        let exponent = int_of_string (String.sub s (i + 1) (String.length s - i - 1)) in
        Printf.sprintf "%se%d" (String.sub s 0 i) exponent
    | None -> if String.contains s '.' then s else s ^ "."

let floats text =
  match float_of_string_opt text with
  | None -> []
  | Some v ->
      [
        ("float_x2", v *. 2.);
        ("float_div2", v /. 2.);
        ("float_zero", 0.);
        ("float_x10", v *. 10.);
        ("float_div10", v /. 10.);
      ]
      |> List.filter (fun (_, w) -> Float.is_finite w && w <> v)
      |> List.map (fun (operator, w) -> (operator, spell w))

let ints text =
  match int_of_string_opt text with
  | None -> []
  | Some n -> [ ("int_plus1", string_of_int (n + 1)); ("int_minus1", string_of_int (n - 1)) ]
