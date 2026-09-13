"""Mesh I/O helpers (VTK format via meshio)."""

from __future__ import annotations

import meshio
import numpy as np
from numpy.typing import NDArray


def save_mesh(
    vertices: NDArray[np.floating],
    triangles: NDArray[np.integer],
    path: str,
) -> None:
    """Write vertices and triangles to a VTK file."""
    mesh = meshio.Mesh(
        points=np.asarray(vertices, dtype=float),
        cells=[("triangle", np.asarray(triangles))],
    )
    meshio.write(path, mesh)


def load_mesh(path: str) -> tuple[NDArray[np.floating], NDArray[np.integer]]:
    """Read (vertices, triangles) from a VTK file."""
    mesh = meshio.read(path)
    triangles = mesh.cells_dict["triangle"]
    return np.asarray(mesh.points, dtype=float), np.asarray(triangles)
