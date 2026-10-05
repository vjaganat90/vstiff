type case = {
  name : string;
  problem : Vstiff.Ode.problem;
  reference : Reference.t;
  tols : float list;
  box : (float * float) array;
}

(* Brusselator-80 stops at 1e-4: 1e-6 takes some 24 s of CPU and 4.3 million calls. The boxes are where the
   transcription check draws states: generic ones, positive where a problem needs it. *)
let cases =
  [
    {
      name = "robertson";
      problem = Problems.Robertson.problem;
      reference = Reference.robertson;
      tols = [ 1e-4; 1e-6; 1e-8 ];
      box = Array.make 3 (0., 1.);
    };
    {
      name = "hires";
      problem = Hires.problem;
      reference = Reference.hires;
      tols = [ 1e-4; 1e-6; 1e-8 ];
      box = Array.make 8 (0., 1.);
    };
    {
      name = "van_der_pol";
      problem = Problems.VanDerPol.problem;
      reference = Reference.van_der_pol;
      tols = [ 1e-4; 1e-6; 1e-8 ];
      box = [| (-2.5, 2.5); (-1., 1.) |];
    };
    {
      name = "brusselator_80";
      problem = Brusselator.problem;
      reference = Reference.brusselator_80;
      tols = [ 1e-4 ];
      box = Array.init 80 (fun k -> if k mod 2 = 0 then (0.5, 4.5) else (1., 5.));
    };
  ]
