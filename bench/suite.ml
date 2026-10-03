type case = { name : string; problem : Vstiff.Ode.problem; reference : Reference.t; tols : float list }

(* The tolerances are Appendix B's. Its vstiff row for Brusselator-80 stops at 1e-4, as this one does: 1e-6 takes some
   24 s of CPU and 4.3 million calls. *)
let cases =
  [
    {
      name = "robertson";
      problem = Problems.Robertson.problem;
      reference = Reference.robertson;
      tols = [ 1e-4; 1e-6; 1e-8 ];
    };
    { name = "hires"; problem = Hires.problem; reference = Reference.hires; tols = [ 1e-4; 1e-6; 1e-8 ] };
    {
      name = "van_der_pol";
      problem = Problems.VanDerPol.problem;
      reference = Reference.van_der_pol;
      tols = [ 1e-4; 1e-6; 1e-8 ];
    };
    { name = "brusselator_80"; problem = Brusselator.problem; reference = Reference.brusselator_80; tols = [ 1e-4 ] };
  ]
