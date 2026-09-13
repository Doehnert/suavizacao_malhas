"""Plotting helpers for mesh visualization."""

from __future__ import annotations

import matplotlib.pyplot as plt

from mesh_smoothing.mesh import Mesh
from mesh_smoothing.solver import manufactured_solution


def plot_mesh(mesh: Mesh, color: str = "r", show: bool = True) -> None:
    """Plot triangle edges and temperature annotations at volume centroids.

    Fictitious volumes are drawn with blue markers, real volumes with green
    markers, preserving the original plot behavior.
    """
    plt.clf()
    for volume in mesh.volumes:
        x, y = volume.centroid()
        annotation = (
            f"{round(manufactured_solution(x, y), 2)} "
            f"{round(volume.temperature, 2)}:{volume.label}"
        )
        plt.plot(x, y, "bo" if volume.fictitious else "go")
        plt.text(x, y, annotation)

        corners = (
            volume.p1.position,
            volume.p2.position,
            volume.p3.position,
            volume.p1.position,
        )
        plt.plot([p[0] for p in corners], [p[1] for p in corners], color)
    if show:
        plt.show()
