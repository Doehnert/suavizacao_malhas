"""Mesh representation and construction for triangular domains."""

from __future__ import annotations

import itertools

import numpy as np
import triangle as tr
from numpy.typing import NDArray

from mesh_smoothing.node import Node
from mesh_smoothing.volume import Volume

# Deliberate vertex displacements applied to the refined "ell" mesh to
# disorder it (original dissertation setup). Net effect of the original
# perturbation block:
#   verts[7][1] -= 0.11
#   verts[16][0] -= 0.1; verts[16][1] += 0.1
#   verts[6] += (0.1, 0.1); verts[8] += (0.1, 0.1)
#   verts[16][0] += 0.1; verts[16][1] += 0.1      # x-offset cancels out
DEFORMED_ELL_PERTURBATIONS: dict[int, tuple[float, float]] = {
    7: (0.0, -0.11),
    6: (0.10, 0.10),
    8: (0.10, 0.10),
    16: (0.0, 0.20),
}


class Mesh:
    """A 2D triangular mesh with nodes, real/fictitious volumes and neighbors.

    Attributes:
        vertices: (n, 3) float array of node coordinates; ghost vertices are
            appended during volume construction.
        triangles: (m, 3) int array of vertex indices (counter-clockwise).
        nodes: List of :class:`Node` objects, one per real vertex.
        volumes: List of :class:`Volume` objects (real + fictitious).
        real_vertices, real_triangles: Geometry snapshot taken before ghost
            vertices were appended; used by :meth:`improve`.
    """

    def __init__(
        self,
        vertices: NDArray[np.floating],
        triangles: NDArray[np.integer],
        vertex_markers: NDArray[np.integer] | None = None,
    ) -> None:
        self.vertices = np.asarray(vertices, dtype=float)
        self.triangles = np.asarray(triangles)
        self.vertex_markers = (
            None if vertex_markers is None else np.asarray(vertex_markers)
        )
        self.nodes: list[Node] = []
        self.volumes: list[Volume] = []
        self.real_vertices = np.array(self.vertices, copy=True)
        self.real_triangles = np.array(self.triangles, copy=True)
        self._rebuild()

    @classmethod
    def create_ell(
        cls,
        area_constraint: float = 0.3,
        perturbations: dict[int, tuple[float, float]] | None = None,
    ) -> Mesh:
        """Build the L-shaped ("ell") domain used in the dissertation.

        The domain is triangulated with the given maximum area constraint,
        scaled to the unit square, deliberately perturbed to disorder it,
        and padded to 3D coordinates (z=0).
        """
        domain = tr.get_data("ell")
        refined = tr.triangulate(domain, f"ra{area_constraint}")
        vertices = np.asarray(refined["vertices"], dtype=float) / 4.0
        triangles = np.asarray(refined["triangles"])
        markers = np.asarray(refined["vertex_markers"])

        for index, (dx, dy) in (perturbations or DEFORMED_ELL_PERTURBATIONS).items():
            vertices[index, 0] += dx
            vertices[index, 1] += dy

        vertices = np.column_stack((vertices, np.zeros(len(vertices))))
        return cls(vertices, triangles, markers)

    @classmethod
    def load(cls, path: str) -> Mesh:  # implemented in Task 5
        """Load a mesh from a VTK file written by :meth:`save`."""
        raise NotImplementedError("meshio I/O lands in Task 5")

    def save(self, path: str) -> None:  # implemented in Task 5
        """Write the current mesh to a VTK file."""
        raise NotImplementedError("meshio I/O lands in Task 5")

    def improve(self, num_steps: int = 1) -> None:
        """Smooth the mesh with the in-package ODT fixed-point iteration.

        [Implemented in Task 5 via mesh_smoothing.smoothing. Temporary body
        keeps the class importable for Task 4 tests.]
        """
        raise NotImplementedError("ODT smoothing lands in Task 5")

    def _rebuild(self) -> None:
        """Recreate the node/volume graph from the current geometry."""
        self.nodes = []
        self.volumes = []
        self._build_nodes()
        self._build_volumes()
        self._build_neighbors()

    def _build_nodes(self) -> None:
        for index, position in enumerate(self.vertices):
            node = Node(index, position)
            if self.vertex_markers is not None and self.vertex_markers[index] == 1:
                node.boundary = True
            self.nodes.append(node)

        for node in self.nodes:
            adjacent: list[Node] = []
            for triangle in self.triangles:
                if node.label in triangle:
                    for index in triangle:
                        neighbor = self.nodes[index]
                        if neighbor.label != node.label and neighbor not in adjacent:
                            adjacent.append(neighbor)
            node.neighbors = adjacent

    def _build_volumes(self) -> None:
        label = 0
        for triangle in self.triangles.tolist():
            n1, n2, n3 = (self.nodes[index] for index in triangle)
            real = Volume(n1, n2, n3, label=label)
            label += 1
            self.volumes.append(real)

            for pair in itertools.combinations(triangle, 2):
                v1 = self.nodes[pair[0]]
                v2 = self.nodes[pair[1]]
                remaining = next(index for index in triangle if index not in pair)
                v3 = self.nodes[remaining]
                pos1, pos2, pos3 = (v1.position, v2.position, v3.position)

                # A fictitious volume mirrors the third vertex across the
                # boundary edge so the ghost volume lies outside the domain.
                if pos1[1] == 1.0 and pos2[1] == 1.0:  # top boundary (y == 1)
                    ghost_position = [pos3[0], pos1[1] + (pos1[1] - pos3[1]), 0.0]
                    self._add_fictitious(v1, v2, ghost_position, label)
                    label += 1
                    continue
                if pos1[0] == 0.0 and pos2[0] == 0.0:  # left boundary (x == 0)
                    ghost_position = [-pos3[0], pos3[1], 0.0]
                    self._add_fictitious(v1, v2, ghost_position, label)
                    label += 1
                    continue
                if pos1[1] == 0.0 and pos2[1] == 0.0:  # bottom boundary (y == 0)
                    ghost_position = [pos3[0], -pos3[1], 0.0]
                    self._add_fictitious(v1, v2, ghost_position, label)
                    label += 1
                    continue
                if pos1[0] == 1.0 and pos2[0] == 1.0:  # right boundary (x == 1)
                    ghost_position = [pos1[0] + (pos1[0] - pos3[0]), pos3[1], 0.0]
                    self._add_fictitious(v1, v2, ghost_position, label)
                    label += 1
                    continue
                if (
                    pos1[0] == 0.5
                    and pos1[1] >= 0.5
                    and pos2[0] == 0.5
                    and pos2[1] >= 0.5
                ):  # inner vertical wall (x == 0.5, y >= 0.5)
                    ghost_position = [pos1[0] + (pos1[0] - pos3[0]), pos3[1], 0.0]
                    if ghost_position[0] > 0.5 and ghost_position[1] > 0.5:
                        self._add_fictitious(v1, v2, ghost_position, label)
                        label += 1
                    continue
                if (
                    pos1[1] == 0.5
                    and pos1[0] >= 0.5
                    and pos2[1] == 0.5
                    and pos2[0] >= 0.5
                ):  # inner horizontal wall (y == 0.5, x >= 0.5)
                    ghost_position = [pos3[0], pos1[1] + (pos1[1] - pos3[1]), 0.0]
                    if ghost_position[1] > 0.5 and ghost_position[0] > 0.5:
                        self._add_fictitious(v1, v2, ghost_position, label)
                        label += 1
                    continue

    def _add_fictitious(
        self,
        v1: Node,
        v2: Node,
        position: list[float],
        label: int,
    ) -> Volume:
        """Create a ghost node + fictitious volume mirroring a boundary edge.

        The ghost node is appended to ``vertices`` only (matching the original
        behavior where ghost nodes were not part of ``nodes``).
        """
        ghost_label = len(self.vertices)
        ghost_node = Node(ghost_label, position)
        self.vertices = np.vstack((self.vertices, np.asarray(position)))
        ghost = Volume(v1, v2, ghost_node, label=label)
        ghost.fictitious = True
        self.volumes.append(ghost)
        return ghost

    def _build_neighbors(self) -> None:
        """Link volumes that share a face (edge-map acceleration of the
        original O(n^2) permutation scan; same neighbor sets)."""
        edge_map: dict[frozenset[int], list[Volume]] = {}
        for volume in self.volumes:
            for a, b in (
                (volume.p1, volume.p2),
                (volume.p2, volume.p3),
                (volume.p3, volume.p1),
            ):
                key = frozenset((a.label, b.label))
                edge_map.setdefault(key, []).append(volume)

        for volume in self.volumes:
            volume.neighbors = []
            for a, b in (
                (volume.p1, volume.p2),
                (volume.p2, volume.p3),
                (volume.p3, volume.p1),
            ):
                key = frozenset((a.label, b.label))
                for other in edge_map[key]:
                    if other is not volume and other not in volume.neighbors:
                        volume.neighbors.append(other)