import numpy as np
import pytest

from mesh_smoothing.mesh import Mesh


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
