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