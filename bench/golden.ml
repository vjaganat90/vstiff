(* The golden table: the figures of one run of bench.exe, written by bench.exe --pin and never edited by hand
   (bench/README.md). Gate judges later runs against it. *)

open Measure

let rows =
  [
    {
      problem = "robertson";
      tol = 0.0001;
      outcome = Done { steps = 261; rejected = 71; rhs_calls = 11195; error = 6.882e-05; scd = 3.10 };
    };
    {
      problem = "robertson";
      tol = 1e-06;
      outcome = Done { steps = 2531; rejected = 831; rhs_calls = 84267; error = 6.485e-07; scd = 5.13 };
    };
    {
      problem = "robertson";
      tol = 1e-08;
      outcome = Done { steps = 25298; rejected = 8423; rhs_calls = 684580; error = 6.111e-09; scd = 7.15 };
    };
    {
      problem = "hires";
      tol = 0.0001;
      outcome = Done { steps = 328; rejected = 95; rhs_calls = 22082; error = 9.519e-05; scd = 1.81 };
    };
    {
      problem = "hires";
      tol = 1e-06;
      outcome = Done { steps = 3212; rejected = 1059; rhs_calls = 187292; error = 1.028e-06; scd = 3.78 };
    };
    {
      problem = "hires";
      tol = 1e-08;
      outcome = Done { steps = 32587; rejected = 10854; rhs_calls = 1773376; error = 9.874e-09; scd = 5.80 };
    };
    {
      problem = "van_der_pol";
      tol = 0.0001;
      outcome = Done { steps = 4145; rejected = 1349; rhs_calls = 116202; error = 2.623e-04; scd = 3.07 };
    };
    {
      problem = "van_der_pol";
      tol = 1e-06;
      outcome = Done { steps = 39038; rejected = 12984; rhs_calls = 978885; error = 2.823e-06; scd = 5.04 };
    };
    {
      problem = "van_der_pol";
      tol = 1e-08;
      outcome = Done { steps = 420058; rejected = 139997; rhs_calls = 9015713; error = 3.126e-08; scd = 6.99 };
    };
    {
      problem = "brusselator_80";
      tol = 0.0001;
      outcome = Done { steps = 737; rejected = 230; rhs_calls = 453298; error = 5.992e-05; scd = 3.73 };
    };
  ]
