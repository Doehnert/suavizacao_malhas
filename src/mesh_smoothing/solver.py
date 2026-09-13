"""Finite volume solver for the manufactured diffusion problem."""

from __future__ import annotations

import itertools

import numpy as np
from numpy.typing import NDArray

from mesh_smoothing.mesh import Mesh
from mesh_smoothing.volume import Volume

PI = np.pi


def manufactured_solution(x: float, y: float) -> float:
    """Analytical temperature T(x, y) = sin(pi x / 2) * sin(pi y / 2)."""
    return float(np.sin(PI * x / 2.0) * np.sin(PI * y / 2.0))


def source_term(x: float, y: float) -> float:
    """Legacy ``Sp`` convention source (malha.py): returns +laplacian(T).

    The physical manufactured source satisfies -laplacian(T) = S with
    S = +(pi^2)/2 * sin(pi x / 2) * sin(pi y / 2); this function returns
    -S (the legacy sign convention), so the solver assembles
    ``rhs = -source_term(...) * area + cross_sum``.
    """
    return float(-(PI**2) / 2.0 * np.sin(PI * x / 2.0) * np.sin(PI * y / 2.0))


def _boundary_value(mesh: Mesh, volume: Volume, x: float, y: float) -> float:
    """Dirichlet value imposed on a fictitious volume from its boundary edge.

    Faithful to legacy ``malha.py``: every matching permutation overwrites
    the value (last match wins), it does not short-circuit on the first
    match. Identical results on the ell domain (at most one condition ever
    matches per ghost here), kept verbatim for port fidelity.
    """
    endpoints = (volume.p1.label, volume.p2.label, volume.p3.label)
    value = 0.0
    for n1, n2 in itertools.permutations(endpoints, 2):
        v1 = mesh.vertices[n1]
        v2 = mesh.vertices[n2]
        if v1[1] == 1.0 and v2[1] == 1.0:  # top edge: T = sin(pi x / 2)
            value = float(np.sin(PI * x / 2.0))
        if (
            v1[0] == 0.5 and v2[0] == 0.5 and v1[1] >= 0.5 and v2[1] >= 0.5
        ):  # inner vertical wall
            value = float(np.sqrt(2.0) / 2.0 * np.sin(PI * x / 2.0))
        if (
            v1[1] == 0.5 and v2[1] == 0.5 and v1[0] >= 0.5 and v2[0] >= 0.5
        ):  # inner horizontal wall
            value = float(np.sin(PI * y / 2.0))
    return value


def solve_diffusion(mesh: Mesh, iterations: int = 10) -> NDArray[np.floating]:
    """Solve the diffusion problem with Picard iteration.

    Each iteration re-assembles the dense linear system because the
    cross-diffusion terms depend on the current temperatures. Returns the
    temperatures for every volume in ``mesh.volumes`` (also stored on each
    ``volume.temperature``). Temperatures are reset to zero on entry so that
    repeated calls are deterministic (the legacy run always started from a
    freshly built mesh, where temperatures default to 0.0).
    """
    for volume in mesh.volumes:
        volume.temperature = 0.0
    count = len(mesh.volumes)
    for _ in range(iterations):
        matrix = np.zeros((count, count))
        rhs = np.zeros(count)
        for index, volume in enumerate(mesh.volumes):
            x, y = volume.centroid()
            if volume.fictitious:
                # Mirror equation: T_ghost + T_neighbor = 2 * T_boundary.
                matrix[index, index] = 1.0
                neighbor = volume.neighbors[0]
                matrix[index, neighbor.label] = 1.0
                rhs[index] = 2.0 * _boundary_value(mesh, volume, x, y)
            else:
                direct_sum = 0.0
                cross_sum = 0.0
                for neighbor in volume.neighbors:
                    cross_sum += volume.cross_diffusion(neighbor, mesh.volumes)
                    direct = volume.direct_diffusion(neighbor)
                    matrix[index, neighbor.label] = -direct
                    direct_sum += direct
                rhs[index] = -source_term(x, y) * volume.area() + cross_sum
                matrix[index, index] = direct_sum

        temperatures = np.linalg.solve(matrix, rhs)
        for volume, temperature in zip(mesh.volumes, temperatures, strict=True):
            volume.temperature = float(temperature)

    return np.array([volume.temperature for volume in mesh.volumes], dtype=float)


def maximum_differences(mesh: Mesh, top: int = 5) -> list[tuple[int, float]]:
    """Return the ``top`` largest |numerical - analytical| differences.

    Evaluated at every volume centroid (including fictitious volumes, as the
    original script did) and sorted descending.
    """
    differences: list[tuple[int, float]] = []
    for volume in mesh.volumes:
        x, y = volume.centroid()
        error = abs(manufactured_solution(x, y) - volume.temperature)
        differences.append((volume.label, float(error)))
    differences.sort(key=lambda item: item[1], reverse=True)
    return differences[:top]
