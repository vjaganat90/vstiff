#!/usr/bin/env python3
"""Regenerate the scipy rows of the work-precision table (docs/plans/roadmap.md, Appendix B) as CSV.

One row per problem, solver (scipy's BDF, Radau and LSODA through solve_ivp) and tolerance, with rtol = atol and the
default finite-difference Jacobian. The right-hand side is wrapped to count its calls: scipy's reported nfev for BDF and
Radau leaves out the evaluations of the finite-difference Jacobian (n per Jacobian and a few more), LSODA's includes
them, and vstiff's Instrument counter includes them all (roadmap, 2.0), so `rhs` here counts every call and `nfev` is
what scipy reports. The error is measured against bench/compare/references.json. The counts are those of one platform,
macOS on arm64: scipy's BDF and LSODA pick their steps from the last bits of the right-hand side, so another libm can
move a row by a step or more (problems.py says where).

    python3 scipy_rows.py [--out FILE]

Columns: problem, solver, rtol, steps, rejected, rhs, nfev, njev, error, scd. `error` is the roadmap's mixed error
max_i |y_i - ref_i| / (1 + |ref_i|) and `scd` is -log10 max_i |y_i - ref_i| / |ref_i|, both at t_end. `rejected` stays
empty: scipy does not report rejections. The bench prints the same columns for vstiff, with nfev and njev empty.
"""

import argparse
import csv
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

from problems import PROBLEMS

REFERENCES = Path(__file__).resolve().parent / "references.json"
SOLVERS = ("BDF", "Radau", "LSODA")
COLUMNS = ("problem", "solver", "rtol", "steps", "rejected", "rhs", "nfev", "njev", "error", "scd")


class Counted:
    """The right-hand side, counting every call it gets, whoever makes it."""

    def __init__(self, rhs):
        self.rhs = rhs
        self.calls = 0

    def __call__(self, t, y):
        self.calls += 1
        return self.rhs(t, y)


def accuracy(y, reference):
    delta = np.abs(y - reference)
    error = float(np.max(delta / (1.0 + np.abs(reference))))
    scd = -math.log10(float(np.max(delta / np.abs(reference))))
    return error, scd


def row(problem, reference, solver, rtol):
    rhs = Counted(problem.rhs)
    sol = solve_ivp(rhs, (0.0, problem.t_end), problem.y0, method=solver, rtol=rtol, atol=rtol)
    if not sol.success:
        sys.exit(f"{problem.name} {solver} rtol {rtol:g} failed: {sol.message}")
    error, scd = accuracy(sol.y[:, -1], reference)
    return {
        "problem": problem.name,
        "solver": solver,
        "rtol": f"{rtol:g}",
        "steps": len(sol.t) - 1,
        "rejected": "",
        "rhs": rhs.calls,
        "nfev": sol.nfev,
        "njev": sol.njev,
        "error": f"{error:.3e}",
        "scd": f"{scd:.2f}",
    }


def rows():
    stored = json.loads(REFERENCES.read_text())["problems"]
    for problem in PROBLEMS.values():
        reference = np.array(stored[problem.name]["stored"]["y"])
        for solver in SOLVERS:
            for rtol in problem.rtols:
                yield row(problem, reference, solver, rtol)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", help="write the CSV to this file instead of standard output")
    args = parser.parse_args()
    # Where the counts come from, for whoever regenerates them elsewhere: another libm can move a row.
    print(f"scipy {scipy.__version__}, numpy {np.__version__}, Python {platform.python_version()}, "
          f"{platform.system()} {platform.machine()}", file=sys.stderr)
    table = list(rows())  # all of them first, so that a solver that fails leaves no half-written file
    out = open(args.out, "w", newline="") if args.out else sys.stdout
    writer = csv.DictWriter(out, fieldnames=COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(table)
    if args.out:
        out.close()


if __name__ == "__main__":
    main()
