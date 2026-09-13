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
