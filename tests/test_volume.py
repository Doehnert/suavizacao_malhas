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
    a = Volume(n0, n1, n2, label=0)
    b = Volume(n0, n1, n3, label=1)
    forward = a.direct_diffusion(b)
    backward = b.direct_diffusion(a)
    assert forward == pytest.approx(1.5)
    # Signed coefficient: reverse direction flips e_xi, so backward == -forward.
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
