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