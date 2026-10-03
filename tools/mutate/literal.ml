(* v as the shortest digits that read back as v, in scientific form: 1.5e-3 is ("1.5", -3). *)
let decimal v =
  let rec shortest precision =
    let s = Printf.sprintf "%.*e" precision v in
    if precision >= 16 || float_of_string s = v then s else shortest (precision + 1)
  in
  let s = shortest 0 in
  let e = String.index s 'e' in
  (String.sub s 0 e, int_of_string (String.sub s (e + 1) (String.length s - e - 1)))

(* Plain notation from 0.01 to 99999 and scientific outside, as the sources write numbers: 20. and 0.05, 1e-8 and 5e6.
   Both are float literals: a plain one keeps its dot, a scientific one needs no dot. *)
let spell (mantissa, exponent) =
  let negative = mantissa.[0] = '-' in
  let unsigned = if negative then String.sub mantissa 1 (String.length mantissa - 1) else mantissa in
  let digits = String.concat "" (String.split_on_char '.' unsigned) in
  let n = String.length digits in
  let body =
    if exponent < -2 || exponent > 4 then Printf.sprintf "%se%d" unsigned exponent
    else if exponent >= n - 1 then digits ^ String.make (exponent - n + 1) '0' ^ "."
    else if exponent >= 0 then
      String.sub digits 0 (exponent + 1) ^ "." ^ String.sub digits (exponent + 1) (n - exponent - 1)
    else "0." ^ String.make (-exponent - 1) '0' ^ digits
  in
  if negative then "-" ^ body else body

(* Ten times and a tenth move the decimal exponent, so 1e-6 becomes 1e-5 and not 9.999999999999999e-6. *)
let floats text =
  match float_of_string_opt text with
  | Some v when Float.is_finite v ->
      let mantissa, exponent = decimal v in
      let shifted k = float_of_string (Printf.sprintf "%se%d" mantissa (exponent + k)) in
      [
        ("float_x2", v *. 2.);
        ("float_div2", v /. 2.);
        ("float_zero", 0.);
        ("float_x10", shifted 1);
        ("float_div10", shifted (-1));
      ]
      |> List.filter (fun (_, w) -> Float.is_finite w && w <> v)
      |> List.map (fun (operator, w) -> (operator, spell (decimal w)))
  | _ -> []

let ints text =
  match int_of_string_opt text with
  | None -> []
  | Some n -> [ ("int_plus1", string_of_int (n + 1)); ("int_minus1", string_of_int (n - 1)) ]
