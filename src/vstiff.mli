(** vstiff: variable-step BDF integration of stiff ordinary differential
    equations. This is the whole public API; the numerical kernel lives in the
    library [Numerics], and [Stage] and [Check] are implementation details. *)

(** Problems, and the contracts methods and step-size controllers meet. *)
module Ode = Ode

(** Named failures: [Diverged], [StepRejected n], [Nan]. *)
module Fail = Fail

(** The resolution of the clock, [16 eps |t|]. *)
module Clock = Clock

(** Counting right-hand-side calls. *)
module Instrument = Instrument

(** Backward Euler, an {!Ode.Method}. *)
module Bdf1 = Bdf1

(** Variable-step BDF2, an {!Ode.Embedded} method. *)
module Bdf2 = Bdf2

(** The halve/double step-size controller, an {!Ode.Controller}. *)
module Halving = Halving

(** Fixed-step integration. *)
module Stepper = Stepper

(** Adaptive integration. *)
module Adaptive = Adaptive
