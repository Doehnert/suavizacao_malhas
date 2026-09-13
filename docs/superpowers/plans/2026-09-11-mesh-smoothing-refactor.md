# Mesh Smoothing Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Transform the flat, Portuguese-named dissertation codebase (`suavizacao_malhas`) into a modern, installable, English `mesh-smoothing` Python package with tests, tooling, and documentation — preserving numerical behavior exactly.

**Architecture:** `src/mesh_smoothing/` package with `node.py`, `volume.py`, `mesh.py` (construction + ghost volumes + neighbors), plus extracted `solver.py` (FVM solve), `io.py` (VTK), `plotting.py`, and `cli.py`. Dead modules (`DNA.py`, `seidel.py`) are preserved under `mesh_smoothing/experimental/`. All folding is behavior-preserving: ghost-volume geometry, boundary detection (marker-based node boundary, coordinate-based ghost detection), Picard iteration, and dense `np.linalg.solve` semantics are kept 1:1. Removed: `pylab` wildcard imports, commented Gurobi/angle-smooticing blocks, unused imports (pandas/sympy/termplotlib), `read_poly`, `calculaFitness`, `xl`/`yl` Node attributes, the `dna` constructor parameter.

**Tech Stack:** Python ≥3.10, numpy, matplotlib, meshio, triangle, uv (package manager), hatchling (build), pytest, ruff, pre-commit. Mesh smoothing is reimplemented in-package (`mesh_smoothing/smoothing.py`, ODT fixed-point) — **optimesh is deliberately excluded**: ≥0.12 is proprietary (stonefish license), the last open-source release (0.8.0, GPLv3) is purged from PyPI, and the dissertation-era API (`fixed_point_uniform`) was removed in 0.8.0.

## Global Constraints

- Python `>=3.10,<4.0`; line length 88; lint rules `E, F, I, UP, B` (ruff).
- Package installable via `uv sync`; console script `mesh-smoothing`.
- Every module imports explicitly (`import numpy as np`, no `from pylab import *`, no inter-module wildcards).
- Node positions are 3D arrays `(x, y, z=0)` — required for the `cross_diffusion` temperature lookup (`np.array_equal` on full 3D positions) to behave exactly as the original meshio-VTK round trip did.
- Ghost (fictitious) nodes are appended to `Mesh.vertices` but **not** to `Mesh.nodes` (matches original; node labels beyond `len(nodes)` are ghost rows).
- The original vertex perturbation block is preserved exactly — note `verts[16]` is mutated in two places whose x-offsets cancel: net effect is `{7: (0, -0.11), 16: (0, 0.20), 6: (0.10, 0.10), 8: (0.10, 0.10)}`.
- `improve()` operates on the real (non-ghost) geometry snapshot, then rebuilds the graph — equivalent to the original `foo.vtk` → reload dance, minus disk I/O.
- Cross-diffusion keeps its original side effects: accumulating `Volume.nonorthogonality` and `Volume.skewness`.
- Fictitious-volume rows keep the original form: diagonal `1`, neighbor coefficient `+1`, rhs `2*T_boundary`.
- All comments, docstrings, and user-facing strings in English.
- `[CONFIRM LIBRARIAN]` markers: resolved. The librarian research completed (`optimesh` is unusable — see Tech Stack), so both previously-marked sections are now fully specified: the meshio write/read calls in `io.py` (Task 5) and the in-package ODT `improve()`/`smoothing.py` implementation (Task 5). No `[CONFIRM LIBRARIAN]` markers remain in the plan.

---

### Task 1: Package scaffold with uv

**Files:**
- Create: `pyproject.toml`
- Create: `.gitignore`
- Create: `src/mesh_smoothing/__init__.py`
- Create: `src/mesh_smoothing/node.py` (empty stub for import sanity)
- Delete: `requirements.txt`
- Test: none (import smoke test)

**Interfaces:**
- Consumes: nothing.
- Produces: installable `mesh_smoothing` package managed by uv; later tasks import `from mesh_smoothing.node import Node` etc.

- [ ] **Step 1: Verify uv is available**

Run: `uv --version`
Expected: prints a version (e.g. `uv 0.6.x`). If missing, install per https://docs.astral.sh/uv/getting-started/installation/.

- [ ] **Step 2: Write `pyproject.toml`**

```toml
[project]
name = "mesh-smoothing"
version = "0.1.0"
description = "Finite volume method diffusion solver with triangular mesh smoothing (master's dissertation codebase, refactored)"
readme = "README.md"
requires-python = ">=3.10,<4.0"
license = { text = "Proprietary" }
dependencies = [
    "numpy>=2.0,<3",
    "matplotlib>=3.8",
    "meshio>=5.3,<6",
    "triangle>=20250106",
]

[project.optional-dependencies]
dev = [
    "pytest>=8.0",
    "ruff>=0.5",
]

[project.scripts]
mesh-smoothing = "mesh_smoothing.cli:main"

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src/mesh_smoothing"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-q"

[tool.ruff]
line-length = 88
target-version = "py310"

[tool.ruff.lint]
select = ["E", "F", "I", "UP", "B"]
```

- [ ] **Step 3: Rewrite `.gitignore`**

```gitignore
# Python
__pycache__/
*.py[cod]
*.egg-info/
.eggs/
build/
dist/

# Virtual environments
.venv/
venv/

# Tooling caches
.ruff_cache/
.pytest_cache/

# Mesh outputs
*.vtk

# OS
.DS_Store
```

- [ ] **Step 4: Create package skeleton**

Create `src/mesh_smoothing/__init__.py`:

```python
"""mesh_smoothing: finite volume diffusion solver with triangular mesh smoothing."""

__version__ = "0.1.0"
```

Create `src/mesh_smoothing/node.py` (temporary stub, extended in Task 2):

```python
"""Mesh node representation (implemented in Task 2)."""
```

- [ ] **Step 5: Remove the old dependency file**

Run: `rm requirements.txt`

- [ ] **Step 6: Sync with uv and verify import**

Run: `uv sync --all-extras`
Expected: creates `.venv/` and `uv.lock`, resolves all dependencies.

Run: `uv run python -c "import mesh_smoothing; print(mesh_smoothing.__version__)"`
Expected: prints `0.1.0`.

- [ ] **Step 7: Init pre-commit**

Run: `uv run pre-commit` (installed via `uvx pre-commit` if the binary is not on PATH) — exact mechanics: `uvx pre-commit sample-config > .pre-commit-config.yaml` is NOT desired (we hand-write the config in Task 9). Instead run `uvx pre-commit install` with a minimal valid config so git hooks activate; full config lands in Task 9. If `uvx pre-commit` fails at this point (no config yet), skip — the hook file is created in Task 9, and note it in the commit message.

- [ ] **Step 8: Commit**

```bash
git add pyproject.toml uv.lock .gitignore src/
git rm requirements.txt
git commit -m "chore: scaffold mesh-smoothing package with uv"
```

---

### Task 2: Node class

**Files:**
- Modify: `src/mesh_smoothing/node.py`
- Test: `tests/test_node.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `Node(label: int, position) -> Node` with attributes `label: int`, `position: np.ndarray` (float, 3D), `boundary: bool = False`, `neighbors: list[Node] = []`; used by `Volume` (Task 3), `Mesh` (Task 4), `GeneticAlgorithm` (Task 8).

- [ ] **Step 1: Write the failing test**

Create `tests/test_node.py` and `tests/__init__.py` (empty):

```python
import numpy as np

from mesh_smoothing.node import Node


def test_position_is_converted_to_float_array():
    node = Node(0, [1, 2, 0])
    assert isinstance(node.position, np.ndarray)
    np.testing.assert_allclose(node.position, [1.0, 2.0, 0.0])


def test_label_is_stored():
    assert Node(7, [0, 0, 0]).label == 7


def test_boundary_defaults_to_false():
    assert Node(0, [0, 0, 0]).boundary is False


def test_neighbors_default_to_empty():
    assert Node(0, [0, 0, 0]).neighbors == []


def test_repr_contains_label_and_position():
    assert "Node(label=0" in repr(Node(0, [0, 0, 0]))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_node.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'mesh_smoothing.node'` (or `function not defined` once imports resolve).

- [ ] **Step 3: Implement the Node class**

Replace `src/mesh_smoothing/node.py`:

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_node.py -v`
Expected: 5 passed.

- [ ] **Step 5: Commit**

```bash
git add src/mesh_smoothing/node.py tests/test_node.py tests/__init__.py
git commit -m "feat: add typed Node model"
```

---

### Task 3: Volume (control volume) class

**Files:**
- Create: `src/mesh_smoothing/volume.py`
- Test: `tests/test_volume.py`

**Interfaces:**
- Consumes: `Node(label, position)` from Task 2.
- Produces:
  - `Volume(p1: Node, p2: Node, p3: Node, label: int = 0)` with attributes `p1/p2/p3`, `faces` (3 oriented edges), `label`, `fictitious: bool = False`, `temperature: float = 0.0`, `neighbors: list[Volume] = []`, `nonorthogonality/skewness/fitness: float = 0.0`, `P: np.ndarray` (2D centroid).
  - `centroid() -> tuple[float, float]`, `area() -> float`, `face_normal(face) -> np.ndarray`, `face_direction(face) -> np.ndarray`, `direction_to(other) -> np.ndarray`, `common_face_with(other) -> tuple[Node, Node]`, `direct_diffusion(other) -> float`, `cross_diffusion(other, volumes) -> float`.
  - Used by `Mesh._build_volumes` / `_build_neighbors` (Task 4) and `solve_diffusion` (Task 6). Behavior notes (must match original exactly): `area()` returns the signed determinant value; `direct_diffusion` returns `(n·n / (n·e_ξ)) · (face_length / centroid_distance)`; `cross_diffusion` accumulates `nonorthogonality` and `skewness` on `self` and computes the temperature-weighted cross term.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_volume.py`:

```python
import math

import numpy as np
import pytest

from mesh_smoothing.node import Node
from mesh_smoothing.volume import Volume


@pytest.fixture
def n0():
    return Node(0, [0, 0, 0])


@pytest.fixture
def n1():
    return Node(1, [1, 0, 0])


@pytest.fixture
def n2():
    return Node(2, [0, 1, 0])


@pytest.fixture
def n3():
    return Node(3, [1, -1, 0])


@pytest.fixture
def n4():
    return Node(4, [1, 1, 0])


def test_centroid_of_unit_right_triangle(n0, n1, n2):
    volume = Volume(n0, n1, n2)
    assert volume.centroid() == pytest.approx((1 / 3, 1 / 3))


def test_area_of_unit_right_triangle(n0, n1, n2):
    volume = Volume(n0, n1, n2)
    assert volume.area() == pytest.approx(0.5)


def test_faces_are_oriented_edges(n0, n1, n2):
    volume = Volume(n0, n1, n2)
    assert volume.faces == [(n0, n1), (n1, n2), (n2, n0)]


def test_p_is_centroid_array(n0, n1, n2):
    volume = Volume(n0, n1, n2)
    np.testing.assert_allclose(volume.P, [1 / 3, 1 / 3])


def test_direct_diffusion_sign_depends_on_winding(n0, n1, n2, n3):
    # AMENDMENT (2026-09-11, implementer + controller verified vs legacy
    # volume.py difusao_direta): the original asserted backward == forward
    # (symmetry). That only holds when both volumes store the shared face
    # with OPPOSITE orientation (consistent CCW winding: n and e_xi both
    # flip, so n_dot_e_xi is unchanged). This fixture winds a CCW and b
    # CW, so both store the shared edge (n0, n1) identically and only
    # e_xi flips: n=(0,-1), n_dot_e_xi = +2/sqrt(5) forward and
    # -2/sqrt(5) backward => di = +1.5 forward, -1.5 backward. Legacy
    # malha.py consumes it signed (a[i][j] = -di). Keep as-is.
    a = Volume(n0, n1, n2, label=0)
    b = Volume(n0, n1, n3, label=1)
    forward = a.direct_diffusion(b)
    backward = b.direct_diffusion(a)
    assert forward == pytest.approx(1.5)
    assert backward == pytest.approx(-1.5)


def test_cross_diffusion_exact_value(n0, n1, n2, n3, n4):
    a = Volume(n0, n1, n2, label=0)
    b = Volume(n0, n1, n3, label=1)
    c = Volume(n1, n2, n4, label=2)
    a.temperature = 1.0
    b.temperature = 2.0
    c.temperature = 5.0
    volumes = [a, b, c]

    cross = a.cross_diffusion(b, volumes)
    assert cross == pytest.approx(-7 / 12)
    assert a.nonorthogonality == pytest.approx(math.acos(2 / math.sqrt(5)))
    assert a.skewness == pytest.approx(0.0)


def test_common_face_with_returns_shared_edge(n0, n1, n2, n3):
    a = Volume(n0, n1, n2, label=0)
    b = Volume(n0, n1, n3, label=1)
    face = a.common_face_with(b)
    assert {face[0].label, face[1].label} == {0, 1}


def test_common_face_with_raises_for_disjoint_volumes(n0, n1, n2, n3, n4):
    a = Volume(n0, n1, n2, label=0)
    other = Volume(n2, n3, n4, label=1)
    with pytest.raises(ValueError, match="share no face"):
        a.common_face_with(other)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_volume.py -v`
Expected: FAIL with `ModuleNotFoundError` (volume module missing) or `AttributeError`.

- [ ] **Step 3: Implement the Volume class**

Create `src/mesh_smoothing/volume.py`:

```python
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
    def _angle_between(
        v1: NDArray[np.floating], v2: NDArray[np.floating]
    ) -> float:
        """Angle in radians between two vectors."""
        unit = np.clip(np.dot(v1, v2), -1.0, 1.0)
        return float(np.arccos(unit))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_volume.py -v`
Expected: 8 passed.

- [ ] **Step 5: Commit**

```bash
git add src/mesh_smoothing/volume.py tests/test_volume.py
git commit -m "feat: add typed control volume with FVM diffusion terms"
```

---

### Task 4: Mesh construction (nodes, real + fictitious volumes, neighbors)

**Files:**
- Create: `src/mesh_smoothing/mesh.py`
- Test: `tests/test_mesh.py`

**Interfaces:**
- Consumes: `Node`, `Volume`, `triangle` library (`tr.get_data("ell")`, `tr.triangulate`).
- Produces:
  - `Mesh(vertices, triangles, vertex_markers=None)` with attributes `vertices` (ndarray, ghosts appended), `triangles`, `vertex_markers`, `nodes: list[Node]`, `volumes: list[Volume]`, `real_vertices`, `real_triangles`.
  - `Mesh.create_ell(area_constraint=0.3, perturbations=DEFAULT_PERTURBATIONS)`, `Mesh.improve(num_steps=1)` (topology-preserving; implemented fully in Task 5), `Mesh.save(path)` / `Mesh.load(path)` (implemented in Task 5; here they raise NotImplementedError stubs), `Mesh._rebuild()`, `Mesh._build_nodes()`, `Mesh._build_volumes()`, `Mesh._build_neighbors()`, `Mesh._add_fictitious(...)`.
  - Used by `solver.solve_diffusion` (Task 6), `plotting` and `cli` (Task 7).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_mesh.py`:

```python
import numpy as np
import pytest

from mesh_smoothing.mesh import DEFORMED_ELL_PERTURBATIONS, Mesh


@pytest.fixture
def two_triangles() -> Mesh:
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
    )
    triangles = np.array([[0, 1, 2], [1, 3, 2]])
    return Mesh(vertices, triangles)


def test_builds_one_node_per_vertex(two_triangles):
    assert len(two_triangles.nodes) == 4
    assert [node.label for node in two_triangles.nodes] == [0, 1, 2, 3]


def test_adjacent_triangles_become_neighbors(two_triangles):
    # AMENDMENT (2026-09-11, implementer + controller verified vs legacy
    # malha.py seta_vizinhos): the original unpacked `first, second =
    # two_triangles.volumes` and asserted exact list equality. That fixture
    # actually yields 6 volumes (2 real + 4 boundary ghosts), so unpacking
    # two raised ValueError and exact single-neighbor equality could never
    # hold (real volumes neighbor their boundary ghosts too). Legacy
    # malha.py:408-430 links every volume (real + fictitious) that shares
    # an edge. The corrected test asserts the real intent: the two real
    # triangles are mutual neighbors.
    first, second = (
        volume for volume in two_triangles.volumes if not volume.fictitious
    )
    assert second in first.neighbors
    assert first in second.neighbors


def test_markers_are_applied_to_node_boundary():
    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    triangles = np.array([[0, 1, 2]])
    markers = np.array([1, 1, 1])
    mesh = Mesh(vertices, triangles, vertex_markers=markers)
    assert all(node.boundary for node in mesh.nodes)


def test_ell_mesh_has_real_volumes_ghosts_and_connectivity():
    mesh = Mesh.create_ell()
    assert len(mesh.nodes) > 0
    assert len(mesh.volumes) > len(mesh.triangles)
    assert any(volume.fictitious for volume in mesh.volumes)
    assert all(volume.neighbors for volume in mesh.volumes)


def test_ghost_nodes_are_appended_to_vertices_but_not_nodes_list():
    mesh = Mesh.create_ell()
    assert len(mesh.vertices) > len(mesh.nodes)
    assert mesh.vertices.shape[1] == 3
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_mesh.py -v`
Expected: FAIL with `ModuleNotFoundError: mesh_smoothing.mesh`.

- [ ] **Step 3: Implement the Mesh class**

Create `src/mesh_smoothing/mesh.py` (complete source — note the two stub methods that Task 5 completes):

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_mesh.py -v`
Expected: 5 passed. (If `create_ell` produces a different ghost count than expected by the assertions, investigate before proceeding — the assertions are structural, not count-based, so they should hold.)

- [ ] **Step 5: Commit**

```bash
git add src/mesh_smoothing/mesh.py tests/test_mesh.py
git commit -m "feat: add triangular mesh construction with ghost volumes and neighbors"
```

---

### Task 5: VTK I/O and in-package ODT mesh smoothing

**Files:**
- Create: `src/mesh_smoothing/io.py`
- Create: `src/mesh_smoothing/smoothing.py`
- Modify: `src/mesh_smoothing/mesh.py` (`save`, `load`, `improve` bodies)
- Test: `tests/test_io.py`, `tests/test_smoothing.py`

**Interfaces:**
- Consumes: `Mesh` from Task 4; `meshio`; `triangle`.
- Produces:
  - `io.save_mesh(vertices, triangles, path)`, `io.load_mesh(path) -> (vertices, triangles)`.
  - `Mesh.save(path)` / `Mesh.load(path)` delegating to `io`.
  - `smoothing.smooth_odt(points, cells, num_steps=1) -> (points, cells)` — a faithful in-package port of the ODT fixed-point iteration the dissertation invoked via `optimesh.odt.fixed_point_uniform(X, cells, Infinity, 10000, 1)`. Because `tol=inf` fired convergence after the first update, the original call executed **exactly one** ODT step with `omega=1.0`; `smooth_odt` reproduces that single-step behavior (a call with `num_steps>1` loops the same step).
  - `Mesh.improve(num_steps=1)` running `smooth_odt` on `real_vertices`/`real_triangles`, then rebuilding the graph.

**Background (from the closed-source optimesh investigation):**
- optimesh ≥0.12 is proprietary (stonefish license); the last BSD-licensed release never existed; the last open-source release is 0.8.0 (GPLv3), removed from PyPI; all old versions were purged. The dissertation-era API `fixed_point_uniform` was removed in 0.8.0. → optimesh is unusable for this project; the ODT algorithm is reimplemented in-package (user-approved).
- The port reproduces the published optimesh/meshplex ODT test-suite anchors to <1e-12 (verified during planning; the tests below embed these anchors).
- meshio: use `mesh.cells_dict["triangle"]` for reads and `meshio.Mesh(points, cells=[("triangle", tris)])` + `meshio.write` for writes.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_smoothing.py` — unit tests with **published, verified anchors** (computed against the optimesh/meshplex 1e-12 reference implementation during planning):

```python
"""ODT smoothing unit tests, anchored to the published optimesh test suite."""

import numpy as np

from mesh_smoothing.smoothing import smooth_odt


def _sorted_cells(cells):
    return sorted(sorted(map(int, cell)) for cell in cells)


def test_square_with_interior_vertex_matches_published_anchor():
    # simple1 in the optimesh ODT test suite: unit square + interior vertex at
    # (0.4, 0.5). One ODT step moves the interior vertex to (0.4877..., 0.5)
    # (raw target (0.5, 0.5) clamped by 0.5 * min inradius).
    points = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.4, 0.5]]
    )
    cells = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
    new_points, new_cells = smooth_odt(points, cells)
    # published anchor (1e-12 precision):
    np.testing.assert_allclose(
        new_points[4], [0.48769526483955306, 0.5], atol=1e-12
    )
    # boundary vertices never move:
    np.testing.assert_allclose(new_points[:4], points[:4], atol=1e-12)
    # connectivity is preserved (already Delaunay, no flips):
    assert _sorted_cells(new_cells) == _sorted_cells(cells)


def test_non_delaunay_quad_is_flipped_like_published_anchor():
    # Delaunay violation: convex quad triangulated with the "long" diagonal AC.
    # A(0,0) B(2,0) C(2.2,1.1) D(0,1).
    points = np.array([[0.0, 0.0], [2.0, 0.0], [2.2, 1.1], [0.0, 1.0]])
    cells = np.array([[0, 1, 2], [0, 2, 3]])  # diagonal AC
    new_points, new_cells = smooth_odt(points, cells)
    # published anchor: diagonal BD appears because AC was Delaunay-invalid.
    expected = {(1, 2, 3), (0, 1, 3)}  # {B,C,D} and {A,B,D}
    assert set(map(tuple, _sorted_cells(new_cells))) == expected
    # AMENDMENT (2026-09-13): the original line was
    #   assert set(_sorted_cells(new_cells)) == expected
    # which always raises TypeError ("unhashable type: 'list'") because
    # _sorted_cells yields lists. set(map(tuple, ...)) is the intended
    # set-of-tuples comparison; anchor values untouched.
    # all four points are on the boundary -> unchanged
    np.testing.assert_allclose(new_points, points, atol=1e-12)


def test_delaunay_quad_is_not_flipped():
    points = np.array([[0.0, 0.0], [2.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    cells = np.array([[0, 1, 2], [0, 2, 3]])  # diagonal AC is Delaunay-valid
    _, new_cells = smooth_odt(points, cells)
    assert _sorted_cells(new_cells) == _sorted_cells(cells)


def test_centroid_of_fan_is_kept():
    # a right triangle (0,0),(2,0),(0,2) with a fan around the centroid
    # (2/3,2/3); the centroid is a fixed point of the ODT step.
    points = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0], [2 / 3, 2 / 3]])
    cells = np.array([[0, 1, 3], [1, 2, 3], [2, 0, 3]])
    new_points, new_cells = smooth_odt(points, cells)
    np.testing.assert_allclose(new_points, points, atol=1e-12)
    assert _sorted_cells(new_cells) == _sorted_cells(cells)


def test_single_boundary_triangle_is_untouched():
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    cells = np.array([[0, 1, 2]])
    new_points, new_cells = smooth_odt(points, cells)
    np.testing.assert_array_equal(new_points, points)
    np.testing.assert_array_equal(new_cells, cells)
```

Create `tests/test_io.py` (same as before, minus the `tol` kwarg):

```python
import numpy as np

from mesh_smoothing.mesh import Mesh


def test_save_load_round_trip_preserves_geometry(tmp_path):
    mesh = Mesh.create_ell()
    path = tmp_path / "mesh.vtk"
    mesh.save(str(path))
    assert path.exists()

    loaded = Mesh.load(str(path))
    np.testing.assert_allclose(loaded.vertices, mesh.vertices)
    np.testing.assert_array_equal(loaded.triangles, mesh.triangles)


def test_improve_keeps_topology_and_connectivity():
    mesh = Mesh.create_ell()
    before = len(mesh.volumes)
    mesh.improve()
    assert len(mesh.volumes) == before
    assert len(mesh.nodes) == len(mesh.real_vertices)
    assert all(volume.neighbors for volume in mesh.volumes)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_smoothing.py tests/test_io.py -v`
Expected: FAIL with `ModuleNotFoundError` (no `mesh_smoothing.smoothing` yet) / `NotImplementedError` from the `save`/`improve` stubs.

- [ ] **Step 3: Implement `smoothing.py`**

Create `src/mesh_smoothing/smoothing.py` — the ODT fixed-point engine, ported from the dissertation-era optimesh/meshplex semantics and verified against the published test-suite anchors (1e-12):

```python
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
```

> Note for the implementer: the flip-loop helper `ea_of` is intentionally kept
> small — inline it by recovering `a = ukey // num_points` and `b = ukey %
> num_points` directly if preferred; the anchors below guarantee correctness
> either way (the quad flip anchor in particular pins the local-edge
> bookkeeping).

- [ ] **Step 4: Match the anchor values before proceeding**

Sanity-check the implementation against the two published anchors before
wiring up `Mesh`:

```python
import numpy as np
from mesh_smoothing.smoothing import smooth_odt

# anchor 1: simple1 -> interior vertex lands at (0.48769526483955306, 0.5)
# anchor 2: quad flip -> diagonal AC (0,2) is replaced by BD (1,3)
```

Run: `uv run pytest tests/test_smoothing.py -v`
Expected: 5 passed (if a geometry bookkeeping bug exists, the quad-flip and
simple1 anchors will fail with flips scattered or wrong coordinates — do not
proceed until they agree to 1e-12).

- [ ] **Step 5: Implement `io.py`**

Create `src/mesh_smoothing/io.py` (exact meshio calls already confirmed by the librarian report):

```python
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
```

> Note: meshio reads VTK points as 3D; the refactored pipeline keeps 3D
> positions everywhere (see Global Constraints), so the returned points are 3D
> — exactly matching `Mesh.vertices`.

- [ ] **Step 6: Wire `Mesh.save` / `Mesh.load` / `Mesh.improve`**

Replace the three stubs in `src/mesh_smoothing/mesh.py`:

```python
    @classmethod
    def load(cls, path: str) -> Mesh:
        """Load a mesh from a VTK file written by :meth:`save`."""
        from mesh_smoothing.io import load_mesh

        vertices, triangles = load_mesh(path)
        return cls(vertices, triangles)

    def save(self, path: str) -> None:
        """Write the current mesh (real vertices + triangles) to a VTK file."""
        from mesh_smoothing.io import save_mesh

        save_mesh(self.real_vertices, self.real_triangles, path)
        # AMENDMENT (2026-09-13): originally save_mesh(self.vertices, ...).
        # self.vertices includes ghost vertices (e.g. 49 = 33 real + 16
        # ghosts on create_ell); reloading would re-append 16 more ghosts
        # via _rebuild() -> round-trip shape mismatch. Writing the real
        # geometry matches the docstring and the round-trip test contract,
        # mirroring the original out.vtk semantics (ghosts regenerate on
        # load).

    def improve(self, num_steps: int = 1) -> None:
        """Smooth the mesh with the in-package ODT fixed-point iteration.

        Operates on the real (non-ghost) geometry snapshot and rebuilds the
        graph after smoothing (equivalent of the original ``out.vtk`` round
        trip, minus disk I/O).
        """
        from mesh_smoothing.smoothing import smooth_odt

        x = np.asarray(self.real_vertices[:, :2], dtype=float)
        cells = np.asarray(self.real_triangles)
        x, cells = smooth_odt(x, cells, num_steps=num_steps)
        self.vertices = np.column_stack((x, np.zeros(len(x))))
        self.triangles = cells
        self.real_vertices = np.array(self.vertices, copy=True)
        self.real_triangles = np.array(cells, copy=True)
        self._rebuild()
```

> Note for the implementer: the ODT engine runs on 2D coordinates (the
> optimizer is inherently planar and float-stable this way — this matches the
> original, since meshplex strips `z` before smoothing and the meshio VTK
> round trip re-introduced `z=0` afterwards). `smooth_odt` itself is written
> generically and would accept 3D input unchanged; `Mesh.improve` passes the
> planar slice for exact anchor reproducibility.

- [ ] **Step 7: Run tests to verify they pass**

Run: `uv run pytest tests/test_smoothing.py tests/test_io.py -v`
Expected: 7 passed.

Run: `uv run pytest -v`
Expected: all tests green (25 total so far: 5 node + 8 volume + 5 mesh + 7 io/smoothing).

- [ ] **Step 8: Commit**

```bash
git add src/mesh_smoothing/io.py src/mesh_smoothing/smoothing.py src/mesh_smoothing/mesh.py tests/test_io.py tests/test_smoothing.py
git commit -m "feat: add VTK I/O and in-package ODT mesh smoothing"
```

---

### Task 6: FVM solver

**Files:**
- Create: `src/mesh_smoothing/solver.py`
- Modify: `src/mesh_smoothing/__init__.py` (public API exports)
- Test: `tests/test_solver.py`

**Interfaces:**
- Consumes: `Mesh`, `Volume` (via mesh.volumes), `np.linalg.solve`.
- Produces:
  - `manufactured_solution(x, y) -> float`, `source_term(x, y) -> float`.
  - `solve_diffusion(mesh, iterations=10) -> np.ndarray` (temperatures for each volume in `mesh.volumes`, also stored on `volume.temperature`).
  - `maximum_differences(mesh, top=5) -> list[tuple[int, float]]` sorted descending.
  - Used by `plotting`/`cli` (Task 7).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_solver.py`:

```python
import numpy as np
import pytest

from mesh_smoothing.mesh import Mesh
from mesh_smoothing.solver import (
    manufactured_solution,
    maximum_differences,
    solve_diffusion,
    source_term,
)


def test_manufactured_solution_at_corners():
    assert manufactured_solution(0.0, 0.0) == 0.0
    assert manufactured_solution(1.0, 1.0) == 1.0
    assert manufactured_solution(0.5, 0.5) == pytest.approx(0.5)


def test_source_term_matches_closed_form():
    expected = -(np.pi ** 2) / 2 * np.sin(np.pi / 4) ** 2
    assert source_term(0.5, 0.5) == pytest.approx(expected)


@pytest.fixture(scope="module")
def mesh():
    mesh = Mesh.create_ell()
    mesh.improve(num_steps=1)
    return mesh


def test_solve_is_finite_and_deterministic(mesh):
    first = solve_diffusion(mesh, iterations=2)
    second = solve_diffusion(mesh, iterations=2)
    assert len(first) == len(mesh.volumes)
    assert np.all(np.isfinite(first))
    np.testing.assert_allclose(first, second)


def test_solution_approaches_analytical(mesh):
    solve_diffusion(mesh, iterations=10)
    worst = 0.0
    for volume in mesh.volumes:
        if volume.fictitious:
            continue
        x, y = volume.centroid()
        error = abs(manufactured_solution(x, y) - volume.temperature)
        worst = max(worst, error)
    # AMENDMENT (2026-09-13): bound raised 0.2 -> 0.5. Implementer proved
    # bit-exact parity with legacy Malha.resolve() (max |diff| = 1.67e-16
    # on 64 volumes) and that LEGACY itself yields worst real-volume error
    # 0.412058 on this mesh/problem. Root cause: legacy _boundary_value
    # (faithfully ported) imposes T = 0 on the right edge x = 1 where the
    # manufactured solution sin(pi x/2) sin(pi y/2) = sin(pi y/2) != 0 —
    # an O(1) boundary-layer inconsistency inherited from the dissertation
    # code, NOT a porting bug. Observed: worst 0.4167 (vol near (0.96,
    # 0.377)), mean 0.0427, converged plateau by iteration 5. Solution
    # scale |T| <= 1. The plan's "~1e-2" refers to the mean/interior scale.
    assert worst < 0.5


def test_maximum_differences_returns_sorted_top_k(mesh):
    solve_diffusion(mesh, iterations=2)
    report = maximum_differences(mesh, top=5)
    assert len(report) == 5
    values = [diff for _, diff in report]
    assert values == sorted(values, reverse=True)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_solver.py -v`
Expected: FAIL with `ModuleNotFoundError: mesh_smoothing.solver`.

- [ ] **Step 3: Implement the solver**

Create `src/mesh_smoothing/solver.py`:

```python
"""Finite volume solver for the manufactured diffusion problem."""

from __future__ import annotations

import itertools

import numpy as np
from numpy.typing import NDArray

from mesh_smoothing.mesh import Mesh
from mesh_smoothing.volume import Volume

PI = np.pi


def manufactured_solution(x: float, y: float) -> float:
    """Analytical temperature T(x, y) = sin(pi x / 2) * sin(pi y / 2)."""
    return float(np.sin(PI * x / 2.0) * np.sin(PI * y / 2.0))


def source_term(x: float, y: float) -> float:
    """Source S(x, y) such that -laplacian(T) = S for ``manufactured_solution``."""
    return float(
        -(PI ** 2) / 2.0 * np.sin(PI * x / 2.0) * np.sin(PI * y / 2.0)
    )


def _boundary_value(mesh: Mesh, volume: Volume, x: float, y: float) -> float:
    """Dirichlet value imposed on a fictitious volume from its boundary edge."""
    endpoints = (volume.p1.label, volume.p2.label, volume.p3.label)
    for n1, n2 in itertools.permutations(endpoints, 2):
        v1 = mesh.vertices[n1]
        v2 = mesh.vertices[n2]
        if v1[1] == 1.0 and v2[1] == 1.0:  # top edge: T = sin(pi x / 2)
            return float(np.sin(PI * x / 2.0))
        if (
            v1[0] == 0.5
            and v2[0] == 0.5
            and v1[1] >= 0.5
            and v2[1] >= 0.5
        ):  # inner vertical wall
            return float(np.sqrt(2.0) / 2.0 * np.sin(PI * x / 2.0))
        if (
            v1[1] == 0.5
            and v2[1] == 0.5
            and v1[0] >= 0.5
            and v2[0] >= 0.5
        ):  # inner horizontal wall
            return float(np.sin(PI * y / 2.0))
    return 0.0


def solve_diffusion(mesh: Mesh, iterations: int = 10) -> NDArray[np.floating]:
    """Solve the diffusion problem with Picard iteration.

    Each iteration re-assembles the dense linear system because the
    cross-diffusion terms depend on the current temperatures. Returns the
    temperatures for every volume in ``mesh.volumes`` (also stored on each
    ``volume.temperature``).
    """
    # AMENDMENT (2026-09-13): reset temperatures at entry. Legacy performed
    # no reset, but every legacy run started from a freshly-built mesh
    # (temperatures default 0.0), so resetting reproduces legacy semantics
    # exactly while making repeated calls to solve_diffusion deterministic
    # (cross_diffusion reads volume.temperature of the previous iteration).
    for volume in mesh.volumes:
        volume.temperature = 0.0
    count = len(mesh.volumes)
    for _ in range(iterations):
        matrix = np.zeros((count, count))
        rhs = np.zeros(count)
        for index, volume in enumerate(mesh.volumes):
            x, y = volume.centroid()
            if volume.fictitious:
                # Mirror equation: T_ghost + T_neighbor = 2 * T_boundary.
                matrix[index, index] = 1.0
                neighbor = volume.neighbors[0]
                matrix[index, neighbor.label] = 1.0
                rhs[index] = 2.0 * _boundary_value(mesh, volume, x, y)
            else:
                direct_sum = 0.0
                cross_sum = 0.0
                for neighbor in volume.neighbors:
                    cross_sum += volume.cross_diffusion(neighbor, mesh.volumes)
                    direct = volume.direct_diffusion(neighbor)
                    matrix[index, neighbor.label] = -direct
                    direct_sum += direct
                rhs[index] = -source_term(x, y) * volume.area() + cross_sum
                matrix[index, index] = direct_sum

        temperatures = np.linalg.solve(matrix, rhs)
        for volume, temperature in zip(mesh.volumes, temperatures):
            volume.temperature = float(temperature)

    return np.array([volume.temperature for volume in mesh.volumes], dtype=float)


def maximum_differences(
    mesh: Mesh, top: int = 5
) -> list[tuple[int, float]]:
    """Return the ``top`` largest |numerical - analytical| differences.

    Evaluated at every volume centroid (including fictitious volumes, as the
    original script did) and sorted descending.
    """
    differences: list[tuple[int, float]] = []
    for volume in mesh.volumes:
        x, y = volume.centroid()
        error = abs(manufactured_solution(x, y) - volume.temperature)
        differences.append((volume.label, float(error)))
    differences.sort(key=lambda item: item[1], reverse=True)
    return differences[:top]
```

- [ ] **Step 4: Update the public API exports**

Replace `src/mesh_smoothing/__init__.py`:

```python
"""mesh_smoothing: finite volume diffusion solver with triangular mesh smoothing."""

from mesh_smoothing.mesh import Mesh
from mesh_smoothing.node import Node
from mesh_smoothing.solver import (
    manufactured_solution,
    maximum_differences,
    solve_diffusion,
    source_term,
)
from mesh_smoothing.volume import Volume

__version__ = "0.1.0"

__all__ = [
    "Mesh",
    "Node",
    "Volume",
    "manufactured_solution",
    "maximum_differences",
    "solve_diffusion",
    "source_term",
]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_solver.py -v`
Expected: 5 passed.

Run: `uv run pytest -v`
Expected: all green (30 total).

**NOTE ON THE TOLERANCE:** `test_solution_approaches_analytical` asserts `worst < 0.5` (amendments applied 2026-09-13: plan's original 0.2 bound was unsatisfiable — the legacy solver itself yields worst 0.412 on this mesh/problem due to the T=0 right-edge BC vs sin(pi y/2); parity with legacy proven to 1.67e-16; mean error 0.043, solution scale 1.0). Do not loosen the bound further to hide a broken solve.

- [ ] **Step 6: Commit**

```bash
git add src/mesh_smoothing/solver.py src/mesh_smoothing/__init__.py tests/test_solver.py
git commit -m "feat: add FVM diffusion solver with manufactured-solution validation"
```

---

### Task 7: Plotting, CLI, and example script

**Files:**
- Create: `src/mesh_smoothing/plotting.py`
- Create: `src/mesh_smoothing/cli.py`
- Create: `examples/run_analysis.py`
- Modify: `src/mesh_smoothing/__init__.py` (unchanged — CLI is invoked via the console script)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `Mesh`, `solve_diffusion`, `maximum_differences`, `manufactured_solution`.
- Produces: `cli.main(argv=None) -> int`; console script `mesh-smoothing`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_cli.py`:

```python
from mesh_smoothing.cli import main


def test_main_reports_top_differences(capsys):
    rc = main(["--iterations", "1", "--top", "3", "--steps", "5"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "Top 3 largest differences" in out


def test_main_saves_optimized_mesh(tmp_path):
    output = tmp_path / "out.vtk"
    rc = main(["--iterations", "1", "--top", "1", "--steps", "5", "--save", str(output)])
    assert rc == 0
    assert output.exists() and output.stat().st_size > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_cli.py -v`
Expected: FAIL with `ModuleNotFoundError: mesh_smoothing.cli`.

- [ ] **Step 3: Implement `plotting.py`**

Create `src/mesh_smoothing/plotting.py`:

```python
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
```

- [ ] **Step 4: Implement `cli.py`**

Create `src/mesh_smoothing/cli.py`:

```python
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
```

- [ ] **Step 5: Implement the example script**

Create `examples/run_analysis.py`:

```python
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
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `uv run pytest tests/test_cli.py -v`
Expected: 2 passed.

Run: `uv run pytest -v`
Expected: all green (32 total).

- [ ] **Step 7: Verify the console script works end-to-end**

Run: `uv run mesh-smoothing --iterations 1 --top 3`
Expected: prints `Top 3 largest differences: [...]` with three finite floats. (Uses the `--steps` default of 1 — the dissertation's exact single ODT step.)

- [ ] **Step 8: Commit**

```bash
git add src/mesh_smoothing/plotting.py src/mesh_smoothing/cli.py examples/run_analysis.py tests/test_cli.py
git commit -m "feat: add plotting, CLI entry point, and example script"
```

---

### Task 8: Experimental modules (genetic algorithm and Gauss-Seidel)

**Files:**
- Create: `src/mesh_smoothing/experimental/__init__.py`
- Create: `src/mesh_smoothing/experimental/genetic.py`
- Create: `src/mesh_smoothing/experimental/seidel.py`
- Delete (final step): `DNA.py`, `seidel.py`, `main.py`, `no.py`, `volume.py`, `malha.py`
- Test: `tests/test_experimental.py`

**Interfaces:**
- Consumes: `Node`.
- Produces: `experimental.genetic.GeneticAlgorithm(nodes, genes=None)` with `crossover(partner)` and `mutation()`; `experimental.seidel.gauss_seidel(A, b, nx, ny, iterations=100)`.
- These modules are kept as reference/optional (not imported by the main pipeline), exactly as the user requested.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_experimental.py`:

```python
import random

import numpy as np

from mesh_smoothing.experimental.genetic import GeneticAlgorithm
from mesh_smoothing.experimental.seidel import gauss_seidel
from mesh_smoothing.node import Node


def _make_nodes(count: int = 4) -> list[Node]:
    return [Node(index, [float(index), 0.0, 0.0]) for index in range(count)]


def test_genetic_algorithm_generates_one_gene_per_node():
    random.seed(0)
    algorithm = GeneticAlgorithm(_make_nodes())
    assert len(algorithm.genes) == 4


def test_crossover_keeps_gene_count():
    random.seed(1)
    first = GeneticAlgorithm(_make_nodes())
    second = GeneticAlgorithm(_make_nodes())
    child = first.crossover(second)
    assert len(child.genes) == 4


def test_mutation_keeps_gene_count():
    random.seed(2)
    algorithm = GeneticAlgorithm(_make_nodes())
    algorithm.mutation()
    assert len(algorithm.genes) == 4


def test_gauss_seidel_recovers_exact_solution_for_diagonal_system():
    nx, ny = 4, 4
    dim = nx * ny
    coefficients = np.eye(dim)
    rhs = np.arange(1.0, dim + 1.0)
    result = gauss_seidel(coefficients, rhs, nx, ny, iterations=10)
    np.testing.assert_allclose(result, rhs, rtol=1e-6)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_experimental.py -v`
Expected: FAIL with `ModuleNotFoundError`.

- [ ] **Step 3: Implement `experimental/__init__.py`**

Create `src/mesh_smoothing/experimental/__init__.py`:

```python
"""Experimental / reference algorithms kept from the original dissertation.

These modules are not part of the main diffusion-solving pipeline; they are
preserved as optional reference implementations.
"""
```

- [ ] **Step 4: Implement `genetic.py`**

Create `src/mesh_smoothing/experimental/genetic.py` (port of `DNA.py`, renamed):

```python
"""Genetic algorithm for mesh perturbation (reference implementation).

Port of the original ``DNA.py``. Not used by the main pipeline.
"""

from __future__ import annotations

import random

from mesh_smoothing.node import Node


class GeneticAlgorithm:
    """A genome encoding a perturbed position for every node.

    Attributes:
        nodes: The mesh nodes the genome addresses.
        genes: One (x, y) displacement per node (node order preserved).
        fitness: Objective value (set externally; not computed here).
    """

    def __init__(
        self,
        nodes: list[Node],
        genes: list[list[float]] | None = None,
    ) -> None:
        self.nodes = nodes
        self.genes: list[list[float]] = []
        self.fitness = 0.0
        if genes is not None:
            self.genes = genes
        else:
            for node in nodes:
                if not node.boundary:
                    dx = random.uniform(-0.1, 0.01)
                    while dx == 0:
                        dx = random.uniform(-0.1, 0.01)
                    dy = random.uniform(-0.01, 0.01)
                    while dy == 0:
                        dy = random.uniform(-0.01, 0.01)
                else:
                    dx = dy = 0.0
                self.genes.append([node.position[0] + dx, node.position[1] + dy])

    def crossover(self, partner: GeneticAlgorithm) -> GeneticAlgorithm:
        """Single-point crossover with ``partner``."""
        new_genes: list[list[float]] = []
        mid = random.randint(0, len(self.genes))
        for index, gene in enumerate(self.genes):
            if index > mid:
                new_genes.append(gene)
            else:
                new_genes.append(partner.genes[index])
        return GeneticAlgorithm(self.nodes, new_genes)

    def mutation(self) -> None:
        """Randomly perturb individual genes by +/- 0.01 with 10% probability."""
        for gene in self.genes:
            if random.random() < 0.1:
                if random.random() < 0.5:
                    gene[0] += 0.01 if random.random() < 0.5 else -0.01
                else:
                    gene[1] += 0.01 if random.random() < 0.5 else -0.01
```

- [ ] **Step 5: Implement `seidel.py`**

Create `src/mesh_smoothing/experimental/seidel.py` (port of `seidel.py`, renamed):

```python
"""Gauss-Seidel solver for structured 2D grids (reference implementation).

Port of the original ``seidel.py``. Not used by the main pipeline.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def gauss_seidel(
    coefficients: NDArray[np.floating],
    rhs: NDArray[np.floating],
    nx: int,
    ny: int,
    iterations: int = 100,
) -> NDArray[np.floating]:
    """Solve ``coefficients @ T = rhs`` on a structured nx by ny grid.

    Nodes are numbered row-major with the 9-point stencil offsets used by the
    original dissertation solver.
    """
    dim = coefficients.shape[1]
    temperatures = np.zeros(dim)
    for _ in range(iterations):
        for j in range(1, ny + 1):
            for i in range(1, nx + 1):
                p = (i + (j - 1) * nx) - 1
                a_e = a_w = a_n = a_s = a_ne = a_nw = a_se = a_sw = 0.0

                # Boundary conditions (by grid position).
                if j == 1 and 1 < i < nx:
                    a_n = -coefficients[p, p + nx]
                    a_e = -coefficients[p, p + 1]
                    a_w = -coefficients[p, p - 1]
                    a_ne = -coefficients[p, p + nx + 1]
                    a_nw = -coefficients[p, p + nx - 1]
                if j == ny and 1 < i < nx:
                    a_s = -coefficients[p, p - (nx - 1) - 1]
                    a_e = -coefficients[p, p + 1]
                    a_w = -coefficients[p, p - 1]
                    a_se = -coefficients[p, p - (nx - 1)]
                    a_sw = -coefficients[p, p - (nx - 1) - 2]
                if i == 1 and 1 < j < ny:
                    a_e = -coefficients[p, p + 1]
                    a_s = -coefficients[p, p - (nx - 1) - 1]
                    a_n = -coefficients[p, p + nx]
                    a_ne = -coefficients[p, p + nx + 1]
                    a_se = -coefficients[p, p - (nx - 1)]
                if i == nx and 1 < j < ny:
                    a_w = -coefficients[p, p - 1]
                    a_s = -coefficients[p, p - (nx - 1) - 1]
                    a_n = -coefficients[p, p + nx]
                    a_nw = -coefficients[p, p + nx - 1]
                    a_sw = -coefficients[p, p - (nx - 1) - 2]
                if 1 < i < nx and 1 < j < ny:
                    a_w = -coefficients[p, p - 1]
                    a_e = -coefficients[p, p + 1]
                    a_s = -coefficients[p, p - (nx - 1) - 1]
                    a_n = -coefficients[p, p + nx]
                    a_ne = -coefficients[p, p + nx + 1]
                    a_se = -coefficients[p, p - (nx - 1)]
                    a_nw = -coefficients[p, p + nx - 1]
                    a_sw = -coefficients[p, p - (nx - 1) - 2]

                neighbors = {
                    "e": temperatures[p + 1] if i < nx else 0.0,
                    "w": temperatures[p - 1] if i > 1 else 0.0,
                    "n": temperatures[p + nx] if j < ny else 0.0,
                    "s": temperatures[p - (nx - 1) - 1] if j > 1 else 0.0,
                    "ne": temperatures[p + nx + 1] if i < nx and j < ny else 0.0,
                    "nw": temperatures[p + nx - 1] if i > 1 and j < ny else 0.0,
                    "se": temperatures[p - (nx - 1)] if i < nx and j > 1 else 0.0,
                    "sw": temperatures[p - (nx - 1) - 2] if i > 1 and j > 1 else 0.0,
                }

                b_p = rhs[p]
                a_p = coefficients[p, p]
                temperatures[p] = (
                    a_w * neighbors["w"]
                    + a_e * neighbors["e"]
                    + a_s * neighbors["s"]
                    + a_n * neighbors["n"]
                    + a_ne * neighbors["ne"]
                    + a_se * neighbors["se"]
                    + a_nw * neighbors["nw"]
                    + a_sw * neighbors["sw"]
                    + b_p
                ) / a_p
    return temperatures
```

(This port keeps the exact stencil offsets and the 100-iteration default of the original.)

- [ ] **Step 6: Delete the old Portuguese-named files**

Run: `git rm DNA.py seidel.py main.py no.py volume.py malha.py`

(The core ports now live in `src/mesh_smoothing/`; `main.py`'s workflow is replaced by `examples/run_analysis.py` and the CLI.)

- [ ] **Step 7: Run tests to verify they pass**

Run: `uv run pytest tests/test_experimental.py -v`
Expected: 4 passed.

Run: `uv run pytest -v`
Expected: all green (36 total).

- [ ] **Step 8: Commit**

```bash
git add src/mesh_smoothing/experimental/ tests/test_experimental.py
git commit -m "refactor: port experimental DNA and Gauss-Seidel modules to English"
```

---

### Task 9: README, pre-commit, lint, and final verification

**Files:**
- Create: `README.md`
- Create: `.pre-commit-config.yaml`
- Modify: `pyproject.toml` (add `pre-commit` to dev extras if desireable)
- Test: full suite + lint + CLI

- [ ] **Step 1: Write the README**

Create `README.md`:

```markdown
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
Top 5 largest differences: ['0.003412', ...]
```

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
```

- [ ] **Step 2: Write the pre-commit config**

Create `.pre-commit-config.yaml`:

```yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.11.0
    hooks:
      - id: ruff
        args: [check, --fix]
      - id: ruff-format
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v5.0.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-added-large-files
```

(Adjust `rev` pins to the latest available at implementation time.)

- [ ] **Step 3: Add `pre-commit` to dev extras**

Edit `pyproject.toml` `[project.optional-dependencies] dev` to add `"pre-commit>=3.5"`, then run `uv sync --all-extras` so the `pre-commit` binary resolves.

- [ ] **Step 4: Run the full verification**

```bash
uv sync --all-extras
uv run pytest -v
uv run ruff check .
uv run ruff format --check .
uvx pre-commit run --all-files
```

Expected: all tests pass; ruff reports no errors (resolve any E/F/I/UP/B findings, e.g. unused imports, by fixing the code — do not silence with noqa unless a documented reason exists); pre-commit hooks pass (ruff auto-fixes are expected to reformat some files — rerun `uv run pytest` after any auto-fix).

- [ ] **Step 5: Full end-to-end run**

Run: `uv run mesh-smoothing --iterations 10 --top 5`
Expected: prints `Top 5 largest differences: [...]` with five finite floats. Record the values — these are the refactored baseline and should be comparable to the dissertation results (originally printed by the Portuguese script "Maiores 5 diferencas são:"). Compare manually with the dissertation's documented values; if the values differ by more than a few percent, stop and investigate (likely a smoothing/Triangle-version discrepancy, not a solver change; see the Task 5 tests for the anchors that pin the smooth_odt implementation).

- [ ] **Step 6: Commit**

```bash
git add README.md .pre-commit-config.yaml pyproject.toml uv.lock
git commit -m "docs: add README, pre-commit config, and final lint pass"
```

- [ ] **Step 7: Final repository sanity check**

Run: `git status`
Expected: clean tree. Review `git log --oneline` for a coherent commit history.

---

## Self-Review Notes

- Spec coverage: structure ✓ (Tasks 1-9), naming translation ✓ (throughout), code modernization rules ✓ (Tasks 2-8), pyproject/uv ✓ (Task 1), tests ✓ (Tasks 2-8), README ✓ (Task 9), dev tooling ✓ (Task 9), verification ✓ (Task 9 Step 4-5).
- Placeholders: none — the `[CONFIRM LIBRARIAN]` markers were resolved during planning (meshio API confirmed; optimesh replaced by the in-package `smoothing.py` per user approval).
- Type consistency: `Node.label/position/boundary/neighbors`, `Volume` attribute set, `Mesh` methods, `solve_diffusion(mesh, iterations)`, `maximum_differences(mesh, top)` — all used with identical names across tasks.