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


def test_gene_perturbation_respects_boundary_nodes():
    random.seed(5)
    interior = Node(3, [0.3, 0.3, 0.0])
    boundary = Node(4, [0.5, 0.0, 0.0])
    boundary.boundary = True
    algorithm = GeneticAlgorithm([interior, boundary])
    interior_gene, boundary_gene = algorithm.genes
    # Non-boundary nodes: perturbed by a nonzero dx in (-0.1, 0.01) and a
    # nonzero dy in (-0.01, 0.01) (review strengthening, Task 8 review).
    dx = interior_gene[0] - 0.3
    dy = interior_gene[1] - 0.3
    assert dx != 0.0 and -0.1 <= dx <= 0.01
    assert dy != 0.0 and -0.01 <= dy <= 0.01
    # Boundary nodes: genes stay exactly on the node position.
    assert boundary_gene == [0.5, 0.0]


def test_gauss_seidel_recovers_exact_solution_for_diagonal_system():
    nx, ny = 4, 4
    dim = nx * ny
    coefficients = np.eye(dim)
    rhs = np.arange(1.0, dim + 1.0)
    result = gauss_seidel(coefficients, rhs, nx, ny, iterations=10)
    np.testing.assert_allclose(result, rhs, rtol=1e-6)
