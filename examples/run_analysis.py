"""End-to-end example: build the L-shaped domain, smooth it, solve, report.

Run with:  uv run python examples/run_analysis.py
"""

from __future__ import annotations

from mesh_smoothing.mesh import Mesh
from mesh_smoothing.solver import maximum_differences, solve_diffusion

TOP = 5


def main() -> None:
    mesh = Mesh.create_ell()
    mesh.improve()
    solve_diffusion(mesh, iterations=10)

    print(f"Top {TOP} largest differences:")
    for label, difference in maximum_differences(mesh, top=TOP):
        print(f"  volume {label}: {difference:.6e}")

    mesh.save("optimized_mesh.vtk")


if __name__ == "__main__":
    main()