#!/usr/bin/env python3
"""Check that the OCaml and the Python right-hand sides are the same function (docs/plans/plan.md, section 7).

bench.exe --transcription prints, for each problem, f(t, y) at 100 seeded points as CSV lines
`problem,t,y_1,...,y_n,f_1,...,f_n`, every number at 17 digits. This script evaluates the right-hand side of
problems.py at the same points and requires every component to agree to 1e-15 relative to the magnitude of its terms.
Without this check a mistranscribed problem would look like a solver failure. From the repository root:

    dune build ./bench/bench.exe
    ./_build/default/bench/bench.exe --transcription | python3 bench/compare/transcription.py [FILE]

It reads standard input, or FILE. The exit status is 1 if a problem is missing, has other than 100 points, or differs.

Why relative to the terms and not to the result: a sum of near-equal terms, or 1 - y^2 near y = 1, loses digits in any
evaluation, so two correct evaluations differ by more than 1e-15 of such a result. The OCaml compiler on arm64 fuses a
product into the sum after it and Python does not, so they do differ: on these points by up to 1.5e-14 of the result
(printed as `plain`) but by 2e-16 of the terms, which is a rounding or two. Built without the fusion, every product
through a function that is never inlined, the two agree bit for bit at these points; not at every point, since
Python's `y ** 2` is a libm `pow`, which on macOS differs from the product in the last bit for some y. A coefficient
that is off by 1e-10 of itself, or more, is caught.
"""

import math
import sys
from collections import defaultdict

import numpy as np

from problems import PROBLEMS

LIMIT = 1e-15
POINTS = 100


class Magnitude:
    """A number that carries, instead of its value, the sum of the absolute values of the terms that make it.

    That sum is what the rounding errors of an evaluation scale with. Evaluating a right-hand side on Magnitudes
    gives, per component, the size that a difference between two evaluations is measured against. Only the operations
    the right-hand sides use are defined: sums, differences and products with numbers or Magnitudes, and integer powers.
    """

    __slots__ = ("m",)

    def __init__(self, m):
        self.m = abs(m)

    @staticmethod
    def of(x):
        return x.m if isinstance(x, Magnitude) else abs(x)

    def __add__(self, other):
        return Magnitude(self.m + Magnitude.of(other))

    __radd__ = __sub__ = __rsub__ = __add__

    def __mul__(self, other):
        return Magnitude(self.m * Magnitude.of(other))

    __rmul__ = __mul__

    def __pow__(self, n):
        return Magnitude(self.m**n)


def terms(rhs, t, y):
    """The magnitude of the terms of each component of rhs(t, y)."""
    return np.array([x.m for x in rhs(t, np.array([Magnitude(v) for v in y], dtype=object))])


def difference(a, b):
    """The absolute difference of two floats, and the relative one: 0 when they are equal, infinite when one is not
    finite."""
    d = abs(a - b)
    if d == 0:
        return 0.0, 0.0
    if not math.isfinite(d):
        return math.inf, math.inf
    return d, d / max(abs(a), abs(b))


def relative_to_terms(d, scale):
    """d over the magnitude of the terms: 0 for no difference, infinite when it cannot be told."""
    if d == 0:
        return 0.0
    return d / scale if math.isfinite(d) and scale > 0 else math.inf


def read(lines):
    """The points of each problem: (t, y, f), the numbers of a line split in half after t."""
    points = defaultdict(list)
    for line in lines:
        fields = line.strip().split(",")
        if len(fields) < 4 or len(fields) % 2 != 0:
            sys.exit(f"cannot read the line {line.strip()[:60]!r}: a name, t, and as many y as f are expected")
        numbers = [float(x) for x in fields[1:]]
        n = (len(numbers) - 1) // 2
        points[fields[0]].append((numbers[0], np.array(numbers[1 : 1 + n]), np.array(numbers[1 + n :])))
    return points


def compare(problem, points):
    """The worst difference of the two right-hand sides, relative to the terms and to the result."""
    worst_terms = worst_plain = 0.0
    for t, y, f_ocaml in points:
        f_python = problem.rhs(t, y)
        scale = terms(problem.rhs, t, y)
        for a, b, s in zip(f_ocaml, f_python, scale):
            d, plain = difference(a, b)
            worst_terms = max(worst_terms, relative_to_terms(d, s))
            worst_plain = max(worst_plain, plain)
    return worst_terms, worst_plain


def main():
    points = read(open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin)
    failures = [f"{name}: no points" for name in PROBLEMS if name not in points]
    failures += [f"{name}: not a problem of problems.py" for name in points if name not in PROBLEMS]
    for name, own in points.items():
        if name not in PROBLEMS:
            continue
        problem = PROBLEMS[name]
        if len(own) != POINTS:
            failures.append(f"{name}: {len(own)} points, expected {POINTS}")
        if any(len(y) != len(problem.y0) for _, y, _ in own):
            failures.append(f"{name}: OCaml and Python differ in the number of components")
            continue
        worst_terms, worst_plain = compare(problem, own)
        verdict = "ok" if worst_terms <= LIMIT else "DIFFERS"
        print(
            f"{name}: {len(own)} points, worst difference {worst_terms:.1e} of the terms"
            f" (limit {LIMIT:.0e}), {worst_plain:.1e} of the result: {verdict}"
        )
        if worst_terms > LIMIT:
            failures.append(f"{name}: the right-hand sides differ by {worst_terms:.1e} of their terms")
    for failure in failures:
        print("FAILED:", failure)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
