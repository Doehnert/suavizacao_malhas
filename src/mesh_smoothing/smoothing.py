"""ODT (optimal Delaunay triangulation) mesh smoothing.

In-package port of the ODT fixed-point iteration that the original
dissertation invoked through ``optimesh.odt.fixed_point_uniform(X, cells,
Infinity, 10000, 1)``. Because ``tol=inf`` converged after the first update,
that single call executed exactly one fixed-point step with ``omega=1``; this
module reproduces that step. The exact formulas (circumcentre barycentric
weights, Heron area in dot-product form, inradius, ``ce``-ratio Delaunay
flips) are taken from the open-source optimesh 0.7.x / meshplex 0.15.x era
and reproduce the published ODT test-suite anchors to 1e-12.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

EPS = 1.0e-12


def _cell_geometry(
    points: NDArray[np.floating], cells: NDArray[np.integer]
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """Per-cell area, circumcentre, barycentre, and inradius.

    Half-edges: e0 = P2-P1, e1 = P0-P2, e2 = P1-P0.
    Area (Heron in dot-product form): A = sqrt(1/4 (d01*d12+d12*d20+d20*d01)),
    with the squared value clamped at 0.
    Circumcentre via barycentric weights beta_k = |e_k|^2 <e_{k+1},e_{k+2}>
    / sum_m |e_m|^2 <e_{m+1},e_{m+2}>.
    """
    p0 = points[cells[:, 0]]
    p1 = points[cells[:, 1]]
    p2 = points[cells[:, 2]]
    e0 = p2 - p1
    e1 = p0 - p2
    e2 = p1 - p0
    d01 = np.einsum("ij,ij->i", e0, e1)
    d12 = np.einsum("ij,ij->i", e1, e2)
    d20 = np.einsum("ij,ij->i", e2, e0)
    n0 = np.einsum("ij,ij->i", e0, e0)
    n1 = np.einsum("ij,ij->i", e1, e1)
    n2 = np.einsum("ij,ij->i", e2, e2)
    vol2 = 0.25 * (d01 * d12 + d12 * d20 + d20 * d01)
    vol2[vol2 < 0] = 0.0
    areas = np.sqrt(vol2)
    beta0 = n0 * d12
    beta1 = n1 * d20
    beta2 = n2 * d01
    denom = beta0 + beta1 + beta2
    circumcentres = (
        (beta0 / denom)[:, None] * p0
        + (beta1 / denom)[:, None] * p1
        + (beta2 / denom)[:, None] * p2
    )
    barycentres = (p0 + p1 + p2) / 3.0
    inradii = 2.0 * areas / (np.sqrt(n0) + np.sqrt(n1) + np.sqrt(n2))
    return areas, circumcentres, barycentres, inradii


def _edge_ce_ratios(
    points: NDArray[np.floating], cells: NDArray[np.integer]
) -> tuple[NDArray[np.floating], NDArray[np.integer]]:
    """Per-row (block-major) ce-ratios and edge endpoint keys.

    Row layout matches _flip_until_delaunay: blocks of local edge 0
    (nodes 1,2), local edge 1 (nodes 2,0), local edge 2 (nodes 0,1).
    ce_k = -<e_{k+1},e_{k+2}>/(4A) for the edge opposite local node k.
    """
    ncells = len(cells)
    edges = np.vstack(
        [cells[:, [1, 2]], cells[:, [2, 0]], cells[:, [0, 1]]]
    )
    ea = np.minimum(edges[:, 0], edges[:, 1])
    eb = np.maximum(edges[:, 0], edges[:, 1])
    keys = ea * len(points) + eb
    ce = _ce_ratios(points, cells)
    return ce.T.ravel(), keys


def _ce_ratios(
    points: NDArray[np.floating], cells: NDArray[np.integer]
) -> NDArray[np.floating]:
    """Per-cell, per-local-edge ce ratios (shape (ncells, 3))."""
    p0 = points[cells[:, 0]]
    p1 = points[cells[:, 1]]
    p2 = points[cells[:, 2]]
    e0 = p2 - p1
    e1 = p0 - p2
    e2 = p1 - p0
    d01 = np.einsum("ij,ij->i", e0, e1)
    d12 = np.einsum("ij,ij->i", e1, e2)
    d20 = np.einsum("ij,ij->i", e2, e0)
    vol2 = 0.25 * (d01 * d12 + d12 * d20 + d20 * d01)
    vol2[vol2 < 0] = 0.0
    areas = np.sqrt(vol2)
    dq = 4.0 * np.maximum(areas, EPS)
    return np.stack([-d12 / dq, -d20 / dq, -d01 / dq], axis=1)


def _flip_until_delaunay(
    points: NDArray[np.floating], cells: NDArray[np.integer]
) -> NDArray[np.integer]:
    """Flip interior edges until every interior edge is locally Delaunay.

    An interior edge (shared by exactly two cells) violates Delaunay iff the
    sum of the two adjacent cells' ce-ratios is negative. Boundary edges
    (shared by one cell) are never flipped.
    """
    cells = np.array(cells, dtype=int)
    num_points = len(points)
    while True:
        ce_flat, keys = _edge_ce_ratios(points, cells)
        ukeys, inverse, counts = np.unique(
            keys, return_inverse=True, return_counts=True
        )
        sums = np.zeros(len(ukeys))
        np.add.at(sums, inverse, ce_flat)
        interior = np.flatnonzero(counts == 2)
        flipped = False
        for u in interior:
            if sums[u] < 0:
                rows = np.flatnonzero(keys == ukeys[u])
                ncells = len(cells)
                j0, k0 = rows[0] % ncells, rows[0] // ncells
                j1, k1 = rows[1] % ncells, rows[1] // ncells
                va, vb = ea_of(keys, rows[0], ukeys[u], num_points)
                c0 = cells[j0][k0]  # opposite vertex in cell j0
                c1 = cells[j1][k1]  # opposite vertex in cell j1
                cells[j0] = [c0, c1, vb]
                cells[j1] = [c0, c1, va]
                flipped = True
                break
        if not flipped:
            return cells


def ea_of(keys, row, ukey, num_points):
    """Recover the sorted endpoint pair (a, b) of an edge from its key.

    Keys are encoded as ``min * num_points + max``, so the min endpoint is
    ``ukey // num_points`` and the max endpoint is ``ukey % num_points``.
    """
    a = ukey // num_points
    b = ukey % num_points
    return a, b


def _boundary_masks(
    cells: NDArray[np.integer], num_points: int
) -> tuple[NDArray[np.integer], NDArray[np.bool_]]:
    """Topological boundary detection (meshplex mark_boundary semantics).

    Returns (is_boundary_cell, is_boundary_point). An edge shared by exactly
    one cell is a boundary edge; a point is on the boundary iff it is an
    endpoint of a boundary edge; a cell is a boundary cell iff it contains at
    least one boundary edge.

    Edge rows are block-major (blocks of local edge 0, 1, 2 over all cells,
    matching _flip_until_delaunay), so the per-cell mask is
    reshape(3, -1).T — NOT reshape(-1, 3).
    """
    edges = np.vstack(
        [cells[:, [0, 1]], cells[:, [1, 2]], cells[:, [2, 0]]]
    )
    ea = np.minimum(edges[:, 0], edges[:, 1])
    eb = np.maximum(edges[:, 0], edges[:, 1])
    keys = ea * num_points + eb
    _, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
    is_boundary_edge = counts[inverse] == 1
    is_boundary_cell = np.any(is_boundary_edge.reshape(3, -1).T, axis=1)
    boundary_edges = edges[is_boundary_edge]
    is_boundary_point = np.zeros(num_points, dtype=bool)
    if len(boundary_edges):
        is_boundary_point[boundary_edges[:, 0]] = True
        is_boundary_point[boundary_edges[:, 1]] = True
    return is_boundary_cell, is_boundary_point


def smooth_odt(
    points: NDArray[np.floating],
    cells: NDArray[np.integer],
    num_steps: int = 1,
) -> tuple[NDArray[np.floating], NDArray[np.integer]]:
    """Run ``num_steps`` ODT fixed-point iterations.

    Each step, in order:
    1. flip_until_delaunay (current connectivity),
    2. compute circumcentres/barycentres/areas/inradii on the current state,
    3. boundary cells use barycentres, interior cells use circumcentres,
    4. new point = area-weighted average of its cells' centres,
    5. reset boundary points to their start-of-step positions,
    6. diff = 1.0 * (new - old), clamped to 0.5 * min adjacent inradius,
    7. move, then flip_until_delaunay again on the moved mesh.

    Returns the final (points, cells); cells reflect post-flip connectivity.
    """
    points = np.array(points, dtype=float)
    cells = np.array(cells, dtype=int)

    for _ in range(num_steps):
        cells = _flip_until_delaunay(points, cells)
        areas, circumcentres, barycentres, inradii = _cell_geometry(
            points, cells
        )
        is_boundary_cell, is_boundary_point = _boundary_masks(
            cells, len(points)
        )
        centres = np.where(is_boundary_cell[:, None], barycentres, circumcentres)
        numerator = np.zeros_like(points)
        denominator = np.zeros(len(points))
        weighted = areas[:, None] * centres
        np.add.at(numerator, cells[:, 0], weighted)
        np.add.at(numerator, cells[:, 1], weighted)
        np.add.at(numerator, cells[:, 2], weighted)
        np.add.at(denominator, cells[:, 0], areas)
        np.add.at(denominator, cells[:, 1], areas)
        np.add.at(denominator, cells[:, 2], areas)
        new_points = numerator / denominator[:, None]
        new_points[is_boundary_point] = points[is_boundary_point]

        diff = new_points - points
        min_inradius = np.full(len(points), np.inf)
        for k in range(3):
            np.minimum.at(min_inradius, cells[:, k], inradii)
        max_step = 0.5 * min_inradius
        magnitude = np.linalg.norm(diff, axis=1)
        too_big = magnitude > max_step
        if np.any(too_big):
            scale = np.ones(len(points))
            scale[too_big] = max_step[too_big] / magnitude[too_big]
            diff = diff * scale[:, None]

        points = points + diff
        cells = _flip_until_delaunay(points, cells)

    return points, cells
