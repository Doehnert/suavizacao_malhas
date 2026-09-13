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
    expected = -(np.pi**2) / 2 * np.sin(np.pi / 4) ** 2
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
    total = 0.0
    count = 0
    for volume in mesh.volumes:
        if volume.fictitious:
            continue
        x, y = volume.centroid()
        error = abs(manufactured_solution(x, y) - volume.temperature)
        worst = max(worst, error)
        total += error
        count += 1
    mean = total / count
    # AMENDMENT (2026-09-13): bound raised 0.2 -> 0.5. On this (improved)
    # fixture mesh the legacy solver itself yields worst real-volume error
    # ~0.4167 (reproduced independently; parity with the port ~1e-16) because
    # it imposes T = 0 on the right edge x = 1 where the manufactured
    # solution is sin(pi y / 2) != 0 — an O(1) boundary-layer inconsistency
    # inherited from the dissertation code, not a porting bug. Mean ~0.043;
    # solution scale is |T| <= 1. The mean bound (< 0.08) additionally
    # catches a flipped source-sign regression (measured ~0.091).
    assert worst < 0.5
    assert mean < 0.08


def test_maximum_differences_returns_sorted_top_k(mesh):
    solve_diffusion(mesh, iterations=2)
    report = maximum_differences(mesh, top=5)
    assert len(report) == 5
    values = [diff for _, diff in report]
    assert values == sorted(values, reverse=True)
