"""Command-line interface for the mesh-smoothing project."""

from __future__ import annotations

import argparse
from pathlib import Path

from mesh_smoothing.mesh import Mesh
from mesh_smoothing.solver import maximum_differences, solve_diffusion


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mesh-smoothing",
        description=(
            "Finite volume diffusion solver on a smoothed triangular mesh "
            "(L-shaped domain, manufactured-solution validation)."
        ),
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=10,
        help="Picard iterations (default: 10)",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=1,
        help="ODT smoothing steps (default: 1 — the dissertation ran a single "
        "fixed-point step, reproduced exactly by smooth_odt)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="number of largest differences to print (default: 5)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="show the mesh plot after solving",
    )
    parser.add_argument(
        "--save",
        type=Path,
        help="write the optimized mesh to a .vtk file",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the dissertation workflow: build, smooth, solve, report."""
    args = build_parser().parse_args(argv)

    mesh = Mesh.create_ell()
    mesh.improve(num_steps=args.steps)
    solve_diffusion(mesh, iterations=args.iterations)

    differences = maximum_differences(mesh, top=args.top)
    printed = [f"{value:.6f}" for _, value in differences]
    print(f"Top {args.top} largest differences: {printed}")

    if args.save is not None:
        mesh.save(str(args.save))
    if args.plot:
        from mesh_smoothing.plotting import plot_mesh

        plot_mesh(mesh)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
