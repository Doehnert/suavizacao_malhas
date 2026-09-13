"""Gauss-Seidel solver for structured 2D grids (reference implementation).

Port of the original ``seidel.py``. Not used by the main pipeline.

The neighbor-value guards below were tightened from the legacy's ``j > 0`` /
``i > 0`` to ``j > 1`` / ``i > 1``: the excluded rows have structurally zero
coefficients (``a_se`` is never set on j == 1, ``a_nw``/``a_sw`` never at
i == 1 / j == 1), so the change is bit-identical to the legacy solver.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def gauss_seidel(
    coefficients: NDArray[np.floating],
    rhs: NDArray[np.floating],
    nx: int,
    ny: int,
    iterations: int = 100,
) -> NDArray[np.floating]:
    """Solve ``coefficients @ T = rhs`` on a structured nx by ny grid.

    Nodes are numbered row-major with the 9-point stencil offsets used by the
    original dissertation solver.
    """
    dim = coefficients.shape[1]
    temperatures = np.zeros(dim)
    for _ in range(iterations):
        for j in range(1, ny + 1):
            for i in range(1, nx + 1):
                p = (i + (j - 1) * nx) - 1
                a_e = a_w = a_n = a_s = a_ne = a_nw = a_se = a_sw = 0.0

                # Boundary conditions (by grid position).
                if j == 1 and 1 < i < nx:
                    a_n = -coefficients[p, p + nx]
                    a_e = -coefficients[p, p + 1]
                    a_w = -coefficients[p, p - 1]
                    a_ne = -coefficients[p, p + nx + 1]
                    a_nw = -coefficients[p, p + nx - 1]
                if j == ny and 1 < i < nx:
                    a_s = -coefficients[p, p - (nx - 1) - 1]
                    a_e = -coefficients[p, p + 1]
                    a_w = -coefficients[p, p - 1]
                    a_se = -coefficients[p, p - (nx - 1)]
                    a_sw = -coefficients[p, p - (nx - 1) - 2]
                if i == 1 and 1 < j < ny:
                    a_e = -coefficients[p, p + 1]
                    a_s = -coefficients[p, p - (nx - 1) - 1]
                    a_n = -coefficients[p, p + nx]
                    a_ne = -coefficients[p, p + nx + 1]
                    a_se = -coefficients[p, p - (nx - 1)]
                if i == nx and 1 < j < ny:
                    a_w = -coefficients[p, p - 1]
                    a_s = -coefficients[p, p - (nx - 1) - 1]
                    a_n = -coefficients[p, p + nx]
                    a_nw = -coefficients[p, p + nx - 1]
                    a_sw = -coefficients[p, p - (nx - 1) - 2]
                if 1 < i < nx and 1 < j < ny:
                    a_w = -coefficients[p, p - 1]
                    a_e = -coefficients[p, p + 1]
                    a_s = -coefficients[p, p - (nx - 1) - 1]
                    a_n = -coefficients[p, p + nx]
                    a_ne = -coefficients[p, p + nx + 1]
                    a_se = -coefficients[p, p - (nx - 1)]
                    a_nw = -coefficients[p, p + nx - 1]
                    a_sw = -coefficients[p, p - (nx - 1) - 2]

                neighbors = {
                    "e": temperatures[p + 1] if i < nx else 0.0,
                    "w": temperatures[p - 1] if i > 1 else 0.0,
                    "n": temperatures[p + nx] if j < ny else 0.0,
                    "s": temperatures[p - (nx - 1) - 1] if j > 1 else 0.0,
                    "ne": temperatures[p + nx + 1] if i < nx and j < ny else 0.0,
                    "nw": temperatures[p + nx - 1] if i > 1 and j < ny else 0.0,
                    "se": temperatures[p - (nx - 1)] if i < nx and j > 1 else 0.0,
                    "sw": temperatures[p - (nx - 1) - 2] if i > 1 and j > 1 else 0.0,
                }

                b_p = rhs[p]
                a_p = coefficients[p, p]
                temperatures[p] = (
                    a_w * neighbors["w"]
                    + a_e * neighbors["e"]
                    + a_s * neighbors["s"]
                    + a_n * neighbors["n"]
                    + a_ne * neighbors["ne"]
                    + a_se * neighbors["se"]
                    + a_nw * neighbors["nw"]
                    + a_sw * neighbors["sw"]
                    + b_p
                ) / a_p
    return temperatures
