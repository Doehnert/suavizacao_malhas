"""Mesh node (vertex) representation."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


class Node:
    """A vertex of the triangular mesh.

    Attributes:
        label: Unique index of the node (row index in
            :attr:`~mesh_smoothing.mesh.Mesh.vertices`).
        position: (x, y, z) coordinates of the node.
        boundary: Whether the node lies on the domain boundary.
        neighbors: Adjacent nodes connected by a triangle edge.
    """

    def __init__(self, label: int, position: NDArray[np.floating]) -> None:
        self.label = label
        self.position = np.asarray(position, dtype=float)
        self.boundary: bool = False
        self.neighbors: list[Node] = []

    def __repr__(self) -> str:
        return f"Node(label={self.label}, position={self.position.tolist()!r})"