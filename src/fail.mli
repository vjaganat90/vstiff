(** The kernel's named failures, re-exported so that users of [Vstiff] never
    need the kernel: [Vstiff.Fail.t] is [Numerics.Fail.t]. *)

include module type of struct
  include Numerics.Fail
end
