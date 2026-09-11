# Mesh Smoothing — Total Refactor Design

**Date:** 2026-09-11
**Status:** Approved by user
**Scope:** Full modernization of the master's dissertation mesh smoothing project (`suavizacao_malhas`)

## Purpose

The existing project is a flat, un-packaged set of 5 Python files implementing a Finite Volume Method (FVM) solver for a diffusion equation on 2D triangular meshes, with mesh optimization via `optimesh`. It was written for a master's dissertation and carries pervasive Portuguese nomenclature, 2019-era pinned dependencies, dead code, wildcard imports, no tests, no README, and no package tooling.

The refactor transforms it into a modern, installable Python package named **`mesh-smoothing`** with:

- English naming throughout (files, classes, methods, variables, comments, strings)
- Proper `src/` layout and package structure
- `uv` as the package manager (`pyproject.toml` build)
- Type hints and docstrings
- Dead code isolated to an optional `experimental/` subpackage (DNA genetic algorithm, Gauss–Seidel solver)
- A pytest test suite for core computation
- Ruff linting + pre-commit dev tooling
- A professional README

## Decisions Locked With User

| Decision | Choice |
|----------|--------|
| Dead code modules (DNA.py, seidel.py) | Keep, moved to `experimental/` subpackage |
| Repository/package name | `mesh-smoothing` (package: `mesh_smoothing`) |
| Polish areas | Full code modernization + README + tests + dev tooling (ruff, pre-commit) |
| Package structure approach | Approach 1: flat `src/mesh_smoothing/` with `experimental/` subpackage |

## Target Package Structure

```
mesh-smoothing/
├── pyproject.toml            # uv-managed, hatchling build, ruff/pytest config
├── README.md
├── .pre-commit-config.yaml
├── .gitignore                # updated for uv, python, tooling
├── src/
│   └── mesh_smoothing/
│       ├── __init__.py       # Public API: Mesh, Node, Volume, FvmSolver
│       ├── node.py           # from no.py — Node class
│       ├── mesh.py           # from malha.py — Mesh class / MeshFactory
│       ├── volume.py         # from volume.py — Volume class
│       ├── solver.py         # FVM solver logic extracted from Mesh.solve()
│       ├── io.py             # VTK write/read helpers extracted from Mesh
│       ├── plotting.py       # plotting helper extracted from Mesh (matplotlib)
│       ├── cli.py            # argparse CLI → runs the dissertation workflow
│       └── experimental/
│           ├── __init__.py
│           ├── genetic.py    # from DNA.py — GeneticAlgorithm
│           └── seidel.py     # from seidel.py — GaussSeidelSolver
├── tests/
│   ├── conftest.py
│   ├── test_node.py
│   ├── test_volume.py
│   ├── test_mesh.py
│   └── test_solver.py
└── examples/
    └── run_analysis.py       # from main.py — end-to-end example workflow
```

## Naming Translation (Portuguese → English)

### Files

| Original | Renamed | Notes |
|----------|---------|-------|
| `malha.py` | `mesh.py` | |
| `no.py` | `node.py` | |
| `volume.py` | `volume.py` | unchanged |
| `DNA.py` | `experimental/genetic.py` | |
| `seidel.py` | `experimental/seidel.py` | |
| `main.py` | `examples/run_analysis.py` | becomes an example, not a CLI |

### Classes / Methods

| Original | Renamed |
|----------|---------|
| `Malha` | `Mesh` |
| `No` | `Node` |
| `Volume2` | `Volume` |
| `DNA` | `GeneticAlgorithm` |
| `melhora_malha()` | `improve()` |
| `carrega_malha()` | `load()` (classmethod) |
| `plotar()` | `plot()` |
| `resolve()` | `solve()` |
| `salva_malha()` | `save()` |
| `centroide()` | `centroid()` |
| `difusao_direta()` | `direct_diffusion()` |
| `difusao_cruzada()` | `cross_diffusion()` |
| `calculaFitness()` | `calculate_fitness()` |

### Variables / Attributes

| Original | Renamed |
|----------|---------|
| `fronteira` | `boundary` |
| `vizinhos` / `vizinho` | `neighbors` / `neighbor` |
| `vertices_reais` | `real_vertices` |
| `volumes_reais` | `real_volumes` |
| `ficticio` | `fictitious` |
| `face_comum` | `common_face` |
| `cont_vertices` | `vertex_count` |
| `novoNo` | `new_node` |
| `nome` | `label` |
| `valor` | `value` |
| `f_pt1..f_pt3` | `fictitious_point_1..3` |

All Portuguese comments, docstrings, and user-facing strings translated to English.

## Code Modernization Rules

1. Remove all `from pylab import *` → explicit `import numpy as np`, `matplotlib.pyplot as plt`, `scipy.spatial` etc. per module need.
2. Remove all wildcard inter-module imports (`from volume import *` → `from mesh_smoothing.volume import Volume`).
3. Add type hints to all public methods, function signatures, and class attributes.
4. Add docstrings (Google style) on all classes and public methods.
5. Remove dead code: ~80 lines commented Gurobi block, ~35 lines commented angle-smoothing block, debug `pass` blocks (`vol.nome == 15`), unused `dna` parameter, dead `calculaFitness` block in malha.py.
6. Replace bare `except: pass` with specific exception handling (Volume.area computation).
7. Use context managers for file I/O.
8. Remove unused imports: pandas, sympy, termplotlib.
9. Fix copy-paste comment errors (all boundary sections said "esquerda").
10. Replace `vol.ficticio == True` with `if vol.fictitious:`.
11. Extract solver logic from `Mesh.solve()` into `solver.py` (`FvmSolver` or function).
12. Extract VTK I/O into `io.py`; `Mesh.save()` writes to a configurable path (no hardcoded `out.vtk`/`foo.vtk`).
13. Extract plotting into `plotting.py`; `main.py`'s matplotlib usage moves there.
14. Move heavy computation out of `__init__`: `Mesh.load()` classmethod for loading, `Mesh.create_ell()` alternative factory building the dissertation domain.
15. Remove `LA`/unused numpy imports (malha.py).
16. Replace hardcoded magic numbers with named constants / parameters (vertex displacement indices 7,16,6,8; area constraint 0.3; iteration count 10000) — surface as constructor/method parameters with documented defaults matching the dissertation result.
17. Optimize the O(n²) neighbor construction loop where practical (keep correctness; simple dict-based acceleration preferred).

## pyproject.toml

- Build backend: `hatchling`
- Project name: `mesh-smoothing`
- Requires-Python: `>=3.10,<4.0`
- Dependencies (direct, relaxed bounds):
  - `numpy>=2.0`
  - `scipy>=1.12`
  - `matplotlib>=3.8`
  - `meshio>=5.3`
  - `optimesh>=3.0`
  - `triangle>=20240115.4`
- `[project.optional-dependencies]` dev: `pytest>=8`, `ruff>=0.5`, `pre-commit>=3`
- `[tool.pytest.ini_options]` testpaths, `[tool.ruff]` line-length 88, target-version py310, select E/F/I/UP/B (a pragmatic ruff default set)
- Entry point: `mesh-smoothing = "mesh_smoothing.cli:main"` — a small CLI that runs the dissertation workflow (generate domain → improve mesh → solve → report max differences → optional plot). `run_analysis.py` example calls the same library API.

## Tests

- `test_node.py`: construction, defaults (boundary False), label/value assignment.
- `test_volume.py`: centroid of a known triangle; area computation; direct/cross diffusion terms on a small known configuration; fictitious-volume single-neighbor behavior.
- `test_mesh.py`: `create_ell()` produces expected vertex/volume counts; `load()` + neighbor connectivity; `improve()` terminates and does not change counts; VTK round-trip save/load preserves coordinates.
- `test_solver.py`: solve on the manufactured solution `T(x,y) = sin(πx/2)·sin(πy/2)` on the L-shaped domain; test asserts (a) solve completes, (b) numerical vs analytical maximum absolute difference is finite, (c) the top-5 largest-difference listing is deterministic across two identical runs.

## README.md Sections

- Title + one-paragraph description (master's dissertation project)
- Features (FVM diffusion solver, triangular mesh optimization, ghost/boundary volumes, manufactured-solution validation)
- Requirements (Python >=3.10, uv)
- Installation (`uv sync`)
- Quick start (`uv run mesh-smoothing` and/or `uv run python examples/run_analysis.py`)
- Project structure tree
- How it works (short algorithm explanation)
- Testing (`uv run pytest`)
- Linting (`uv run ruff check .`)
- License note

## Dev Tooling

- Ruff: lint + format, configured in pyproject.toml
- Pre-commit: `ruff check`, `ruff format --check`, `pytest` hooks (or suitable subset)
- `.gitignore`: extend/generate for uv, `.venv`, `__pycache__`, build artifacts, `.ruff_cache`, `.pytest_cache`, `*.vtk` outputs

## Verification

After implementation:

1. `uv sync` installs cleanly (resolves modern dependency versions).
2. `uv run ruff check .` passes (or documented known-warnings).
3. `uv run pytest` passes (tests above).
4. `uv run mesh-smoothing` runs end-to-end and produces the reported top-5 largest differences; max difference is finite and comparable in magnitude to the original (dissertation numbers preserved).
5. `uv run pre-commit run --all-files` passes.
6. README instructions can be followed verbatim by a fresh user.

## Out of Scope

- Changing the numerical scheme or results fidelity (refactor preserves behavior).
- Migrating to another mesh library (PyVista, GMSH) — `triangle` + `optimesh` are preserved as core deps.
- GUI, notebooks, or cloud deployment.
- Rewriting the dissertation itself; this is the accompanying codebase.