"""The bench problems as scipy sees them.

Each right-hand side repeats the expression order of its OCaml twin (test/problems/problems.ml, bench/hires.ml,
bench/brusselator.ml), so the two agree to a rounding or two, and transcription.py checks that. The twins are written
apart on purpose: a typo on one side must show as a difference, not hide in a copy.

The step counts of scipy's BDF and LSODA follow the last bits of the right-hand side, so two spellings are the ones that
reproduce the roadmap's Appendix B on the machine that measured it, macOS on arm64: a product in `robertson` and a
power in `van_der_pol` (the second differs from its OCaml twin's product). Where numpy's power is libm's `pow`, as
there, it differs from the product in the last bit for some arguments; a libm whose `pow` is correctly rounded gives
the product, and another row.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class Problem:
    name: str  # as the OCaml bench names it in its CSV
    rhs: Callable[[float, np.ndarray], np.ndarray]
    y0: np.ndarray
    t_end: float
    rtols: tuple  # the tolerances of the work-precision table (roadmap, Appendix B); rtol = atol


def robertson(t, y):
    # y2 * y2, not y2 ** 2: with the power scipy's LSODA takes 241 and 443 steps at rtol 1e-6 and 1e-8 instead of the
    # 245 and 467 of Appendix B (the two spellings differ in the last bit for some y2).
    a = 0.04 * y[0]
    b = 1e4 * y[1] * y[2]
    c = 3e7 * y[1] * y[1]
    return np.array([b - a, a - b - c, c])


def hires(t, y):
    y1, y2, y3, y4, y5, y6, y7, y8 = y
    return np.array(
        [
            -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007,
            1.71 * y1 - 8.75 * y2,
            -10.03 * y3 + 0.43 * y4 + 0.035 * y5,
            8.32 * y2 + 1.71 * y3 - 1.12 * y4,
            -1.745 * y5 + 0.43 * y6 + 0.43 * y7,
            -280.0 * y6 * y8 + 0.69 * y4 + 1.71 * y5 - 0.43 * y6 + 0.69 * y7,
            280.0 * y6 * y8 - 1.81 * y7,
            -280.0 * y6 * y8 + 1.81 * y7,
        ]
    )


def van_der_pol(t, y):
    # y1 ** 2, not y1 * y1: with the product scipy's BDF takes 1819 steps at rtol 1e-8, not the 1846 of Appendix B.
    # The twin in test/problems/problems.ml multiplies; it belongs to the corpus and stays as it is.
    mu = 1000.0
    return np.array([y[1], mu * (1.0 - y[0] ** 2) * y[1] - y[0]])


# The Brusselator with diffusion on (0, 1), by the method of lines: u' = A + u^2 v - (B + 1) u + alpha u_xx and
# v' = B u - u^2 v + alpha v_xx, with u = 1 and v = 3 at both ends and u(x, 0) = 1 + sin(2 pi x), v(x, 0) = 3, on 40
# interior cells. The state interleaves the cells: (u_1, v_1, u_2, v_2, ...).
CELLS = 40
ALPHA = 1.0 / 50.0
A = 1.0
B = 3.0
DIFFUSION = ALPHA * float((CELLS + 1) * (CELLS + 1))  # alpha / dx^2 with dx = 1 / (CELLS + 1)


def brusselator(t, y):
    u = y[0::2]
    v = y[1::2]
    u_ext = np.concatenate(([1.0], u, [1.0]))
    v_ext = np.concatenate(([3.0], v, [3.0]))
    du = A + u * u * v - (B + 1.0) * u + DIFFUSION * (u_ext[:-2] - 2.0 * u + u_ext[2:])
    dv = B * u - u * u * v + DIFFUSION * (v_ext[:-2] - 2.0 * v + v_ext[2:])
    out = np.empty_like(y)
    out[0::2] = du
    out[1::2] = dv
    return out


def brusselator_y0():
    x = (np.arange(CELLS) + 1) * (1.0 / (CELLS + 1))
    y0 = np.empty(2 * CELLS)
    y0[0::2] = 1.0 + np.sin(2.0 * np.pi * x)
    y0[1::2] = 3.0
    return y0


PROBLEMS = {
    p.name: p
    for p in (
        Problem("robertson", robertson, np.array([1.0, 0.0, 0.0]), 1e4, (1e-4, 1e-6, 1e-8)),
        Problem("hires", hires, np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057]), 321.8122, (1e-4, 1e-6, 1e-8)),
        Problem("van_der_pol", van_der_pol, np.array([2.0, 0.0]), 2000.0, (1e-4, 1e-6, 1e-8)),
        Problem("brusselator_80", brusselator, brusselator_y0(), 10.0, (1e-4, 1e-6)),
    )
}
