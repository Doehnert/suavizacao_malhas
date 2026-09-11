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