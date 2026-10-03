(** Generators with integrated shrinking. A generator draws a value from a [Random.State.t] together with a lazy
    tree of smaller candidates, so every combinator shrinks what it builds and no property writes a shrinker.
    docs/testing.md, Properties. *)

(** A value and its shrink candidates, simplest first; each candidate carries its own. *)
type 'a tree = Node of 'a * 'a tree Seq.t

(** Draws one value and its tree; equal states draw equal trees. *)
type 'a t = Random.State.t -> 'a tree

(** The value at the top of a tree. *)
val root : 'a tree -> 'a

(** Always the value, which does not shrink. *)
val return : 'a -> 'a t

(** [map f g] applies [f] to every value in [g]'s tree. *)
val map : ('a -> 'b) -> 'a t -> 'b t

(** [bind g f] draws [a] from [g], then from [f a]. It shrinks [a] first, drawing [f] again from the same state, then
    what [f a] drew. *)
val bind : 'a t -> ('a -> 'b t) -> 'b t

(** [let*] is [bind] and [let+] is [map]. *)
module Syntax : sig
  val ( let* ) : 'a t -> ('a -> 'b t) -> 'b t
  val ( let+ ) : 'a t -> ('a -> 'b) -> 'b t
end

(** Tuples shrink one component at a time, the first one first. *)
val pair : 'a t -> 'b t -> ('a * 'b) t

val triple : 'a t -> 'b t -> 'c t -> ('a * 'b * 'c) t

(** Uniform in [lo, hi], shrinking toward [lo]. *)
val int : lo:int -> hi:int -> int t

(** Uniform in [lo, hi], shrinking toward the point of the range nearest 0 by dropping significant bits. *)
val float : lo:float -> hi:float -> float t

(** A magnitude in [lo, hi], [0 < lo <= hi], whose logarithm is uniform; it shrinks toward the magnitude nearest 1,
    through powers of two. *)
val log_uniform : lo:float -> hi:float -> float t

(** One of +0, -0, the smallest subnormal, +infinity, -infinity, nan, [max_float] and [epsilon], shrinking to 0. *)
val special : float t

(** [array len g] is as many values of [g] as [len] draws. It shrinks by dropping elements, to the lengths that
    [len] shrinks to, and by shrinking the elements. *)
val array : int t -> 'a t -> 'a array t

(** Draws from one generator of a non-empty list, shrinking toward the earlier ones. *)
val one_of : 'a t list -> 'a t

(** [prescribed sigma] is [U diag(sigma) V^T], [U] and [V] products of random Householder reflections, so its
    singular values are the [|sigma_i|]. Shrinks toward fewer and simpler reflections. *)
val prescribed : float array -> float array array t

(** [singular n], [n >= 2]: an [n] by [n] matrix with entries in [-1, 1] and one row repeated, so it is exactly
    singular. *)
val singular : int -> float array array t
