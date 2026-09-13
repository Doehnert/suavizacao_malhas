# mesh-smoothing

Finite volume method (FVM) diffusion solver with triangular mesh smoothing.
This package is the refactored companion code of a master's dissertation:
it solves a manufactured diffusion problem on a smoothed unstructured
triangular mesh, with finite volume discretization, ghost (fictitious)
volumes for boundary conditions, and validation against the analytical
solution.

## Features

- L-shaped domain meshed with the `triangle` library and smoothed with the
  built-in ODT fixed-point smoother (in-package port of the algorithm the
  dissertation ran via optimesh)
- Finite volume discretization with direct and cross-diffusion terms
- Fictitious (ghost) volumes for Dirichlet boundary treatment
- Picard iteration over the cross-diffusion coupling
- Manufactured-solution validation: `T(x, y) = sin(πx/2) sin(πy/2)`
- Reports the largest numerical-vs-analytical differences
- VTK export of the optimized mesh (`meshio`)
- Optional `experimental/` reference algorithms (genetic algorithm and a
  Gauss-Seidel solver for structured grids)

## Requirements

- Python ≥ 3.10
- [uv](https://docs.astral.sh/uv/) (package manager)

## Installation

```bash
uv sync --all-extras
```

## Quick start

Run the full workflow (build domain → smooth mesh → solve → report):

```bash
uv run mesh-smoothing
```

You should see the top 5 largest differences between the numerical and
analytical solutions, e.g.:

```
Top 5 largest differences: ['0.697114', '0.416746', '0.282053', '0.239332', '0.234416']
```

(Measured baseline on the current mesh, `--iterations 10 --top 5`. The worst
errors sit near the right edge x=1 where the legacy boundary treatment
imposes T=0 but the manufactured solution is sin(pi y/2) != 0 — an O(1)
boundary layer inherited verbatim from the dissertation code.)

Or run the example script directly:

```bash
uv run python examples/run_analysis.py
```

Options for the CLI:

```bash
uv run mesh-smoothing --iterations 10 --top 5 --plot --save optimized_mesh.vtk
```

## Project structure

```
src/mesh_smoothing/
├── __init__.py      # public API
├── node.py          # mesh nodes
├── volume.py        # control volumes and FVM diffusion terms
├── mesh.py          # mesh construction, ghost volumes, smoothing
├── smoothing.py     # in-package ODT fixed-point smoother
├── solver.py        # FVM solver and manufactured-solution validation
├── io.py            # VTK read/write
├── plotting.py      # matplotlib visualization
├── cli.py           # command-line entry point
└── experimental/    # reference algorithms (not used by the pipeline)
    ├── genetic.py
    └── seidel.py
```

## How it works

1. The L-shaped domain is triangulated with a maximum area constraint and
   scaled to the unit square; a few vertices are deliberately perturbed to
   disorder the mesh (as in the dissertation setup).
2. The built-in ODT fixed-point smoother (`smoothing.py`) improves element
   quality — an in-package port of the algorithm the dissertation ran through
   optimesh.
3. Boundary edges get fictitious volumes that mirror the interior across the
   boundary, providing the Dirichlet data through a mirror equation.
4. The diffusion equation is discretized with direct and cross-diffusion
   terms and solved with a dense linear solver in a Picard loop.
5. The numerical solution is compared with the manufactured analytical
   solution at every volume centroid.

## Testing

```bash
uv run pytest
```

## Linting and formatting

```bash
uv run ruff check .
uv run ruff format --check .
```

## Pre-commit

```bash
uv run pre-commit install
```

## License

Proprietary — this code was written for a master's dissertation.
