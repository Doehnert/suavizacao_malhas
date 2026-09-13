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
