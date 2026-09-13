"""Control volume (triangle) of the finite volume method."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from mesh_smoothing.node import Node


def _perp(vector: NDArray[np.floating]) -> NDArray[np.floating]:
    """Return a 2D vector perpendicular to ``vector`` (90-degree rotation)."""
    return np.array([-vector[1], vector[0]])


def _segment_intersection(
    a1: NDArray[np.floating],
    a2: NDArray[np.floating],
    b1: NDArray[np.floating],
    b2: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Intersection point of segments a1-a2 and b1-b2."""
    delta_a = a2 - a1
    delta_b = b2 - b1
    delta_p = a1 - b1
    denom = np.dot(_perp(delta_a), delta_b)
    num = np.dot(_perp(delta_a), delta_p)
    return (num / float(denom)) * delta_b + b1


class Volume:
    """A triangular control volume defined by three nodes.

    Attributes:
        p1, p2, p3: The :class:`Node` objects at the triangle corners
            (counter-clockwise order).
        faces: The three oriented edges ``(p1, p2), (p2, p3), (p3, p1)``.
        label: Unique index of the volume in the mesh.
        fictitious: True for ghost volumes introduced at the boundary.
        temperature: Solution value stored at the volume centroid.
        neighbors: Adjacent volumes sharing one face.
        nonorthogonality: Accumulated face non-orthogonality angles.
        skewness: Accumulated face skewness ratios.
        fitness: Genetic-algorithm fitness (not used by the solver).
        P: Centroid coordinates as a 2D numpy array.
    """

    def __init__(self, p1: Node, p2: Node, p3: Node, label: int = 0) -> None:
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.faces = [(self.p1, self.p2), (self.p2, self.p3), (self.p3, self.p1)]
        self.label = label
        self.fictitious: bool = False
        self.temperature: float = 0.0
        self.neighbors: list[Volume] = []
        self.nonorthogonality = 0.0
        self.skewness = 0.0
        self.fitness = 0.0
        x, y = self.centroid()
        self.P = np.array([x, y])

    def __repr__(self) -> str:
        return f"Volume(label={self.label}, fictitious={self.fictitious})"

    def centroid(self) -> tuple[float, float]:
        """Return the centroid coordinates ``(x, y)`` of the triangle."""
        x = (self.p1.position[0] + self.p2.position[0] + self.p3.position[0]) / 3
        y = (self.p1.position[1] + self.p2.position[1] + self.p3.position[1]) / 3
        return float(x), float(y)

    def area(self) -> float:
        """Triangle area via the signed determinant formula (as in the original)."""
        xa, xb, xc = (
            self.p1.position[0],
            self.p2.position[0],
            self.p3.position[0],
        )
        ya, yb, yc = (
            self.p1.position[1],
            self.p2.position[1],
            self.p3.position[1],
        )
        matrix = [[xa, xb, xc], [ya, yb, yc], [1, 1, 1]]
        return 0.5 * np.linalg.det(matrix)

    def face_normal(self, face: tuple[Node, Node]) -> NDArray[np.floating]:
        """Unit normal vector of a face (pointing left of its direction)."""
        xa, ya = face[0].position[0], face[0].position[1]
        xb, yb = face[1].position[0], face[1].position[1]
        delta = np.hypot(xb - xa, yb - ya)
        return np.array([(yb - ya) / delta, -(xb - xa) / delta])

    def face_direction(self, face: tuple[Node, Node]) -> NDArray[np.floating]:
        """Unit vector along a face (from face[0] toward face[1])."""
        xa, ya = face[0].position[0], face[0].position[1]
        xb, yb = face[1].position[0], face[1].position[1]
        delta = np.hypot(xb - xa, yb - ya)
        return np.array([(xb - xa) / delta, (yb - ya) / delta])

    def direction_to(self, other: Volume) -> NDArray[np.floating]:
        """Unit vector from this volume's centroid toward ``other``'s centroid."""
        xp, yp = self.centroid()
        xa, ya = other.centroid()
        delta = np.hypot(xa - xp, ya - yp)
        return np.array([(xa - xp) / delta, (ya - yp) / delta])

    def common_face_with(self, other: Volume) -> tuple[Node, Node]:
        """Return the oriented face shared with ``other``.

        Raises:
            ValueError: If the volumes do not share a face.
        """
        for face in self.faces:
            for face_other in other.faces:
                if (
                    face[0].label == face_other[0].label
                    and face[1].label == face_other[1].label
                ) or (
                    face[1].label == face_other[0].label
                    and face[0].label == face_other[1].label
                ):
                    return face
        raise ValueError(f"Volume {self.label} and {other.label} share no face")

    def direct_diffusion(self, other: Volume) -> float:
        """Direct-diffusion coefficient across the face shared with ``other``."""
        face = self.common_face_with(other)
        normal = self.face_normal(face)
        direction = self.direction_to(other)

        a = face[0].position
        b = face[1].position
        face_length = np.hypot(b[0] - a[0], b[1] - a[1])

        xp, yp = self.centroid()
        xa, ya = other.centroid()
        centroid_distance = np.hypot(xa - xp, ya - yp)

        coefficient = (normal.dot(normal) / normal.dot(direction)) * (
            face_length / centroid_distance
        )
        return float(coefficient)

    def cross_diffusion(self, other: Volume, volumes: list[Volume]) -> float:
        """Cross-diffusion term across the face shared with ``other``.

        Accumulates ``self.nonorthogonality`` (angle between the face normal
        and the centroid direction) and ``self.skewness`` (ratio of the face
        midpoint deviation to the centroid distance) as side effects, exactly
        like the original implementation.
        """
        face = self.common_face_with(other)
        n = self.face_normal(face)
        direction = self.direction_to(other)
        face_dir = self.face_direction(face)

        self.nonorthogonality += self._angle_between(n, direction)

        a = face[0].position[:2]
        b = face[1].position[:2]

        intersection = _segment_intersection(self.P, other.P, a, b)
        midpoint = (a + b) / 2
        deviation = midpoint - intersection
        centroid_vector = other.P - self.P
        self.skewness += np.linalg.norm(deviation) / np.linalg.norm(centroid_vector)

        a_3d = np.append(a, 0.0)
        b_3d = np.append(b, 0.0)
        t_side_b: list[float] = []
        t_side_a: list[float] = []
        for volume in volumes:
            node_positions = (
                volume.p1.position,
                volume.p2.position,
                volume.p3.position,
            )
            if any(np.array_equal(b_3d, pos) for pos in node_positions):
                t_side_b.append(volume.temperature)
            if any(np.array_equal(a_3d, pos) for pos in node_positions):
                t_side_a.append(volume.temperature)

        t_b_mean = sum(t_side_b) / len(t_side_b)
        t_a_mean = sum(t_side_a) / len(t_side_a)

        cross = -(direction.dot(face_dir) / n.dot(direction)) * (t_b_mean - t_a_mean)
        return float(cross)

    @staticmethod
    def _angle_between(v1: NDArray[np.floating], v2: NDArray[np.floating]) -> float:
        """Angle in radians between two vectors."""
        unit = np.clip(np.dot(v1, v2), -1.0, 1.0)
        return float(np.arccos(unit))
