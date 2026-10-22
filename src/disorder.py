"""Random disorder realizations for the SSH model."""

import numpy as np


def sample_hopping_disorder(
    n_cells: int, t1: float, t2: float, strength: float, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """Sample independent uniform bond variations in [-strength, strength]."""
    intracell = t1 + rng.uniform(-strength, strength, size=n_cells)
    intercell = t2 + rng.uniform(-strength, strength, size=n_cells - 1)
    return intracell, intercell


def sample_onsite_disorder(
    n_cells: int, strength: float, rng: np.random.Generator
) -> np.ndarray:
    """Sample independent uniform site energies in [-strength, strength]."""
    return rng.uniform(-strength, strength, size=2 * n_cells)
