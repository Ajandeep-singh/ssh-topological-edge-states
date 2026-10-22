"""Hamiltonian construction and observables for the 1D SSH chain."""

from __future__ import annotations

import numpy as np


def build_ssh_hamiltonian(
    n_cells: int,
    t1: float,
    t2: float,
    *,
    intracell_hoppings: np.ndarray | None = None,
    intercell_hoppings: np.ndarray | None = None,
    onsite_energies: np.ndarray | None = None,
) -> np.ndarray:
    """Build an open SSH-chain Hamiltonian in the site order A1, B1, ... .

    ``intracell_hoppings`` has length N and ``intercell_hoppings`` length
    N-1. ``onsite_energies`` has length 2N. Missing arrays use the clean
    values (or zero onsite energy).
    """
    if not isinstance(n_cells, (int, np.integer)) or n_cells < 1:
        raise ValueError("n_cells must be a positive integer")

    t1_bonds = np.full(n_cells, t1, dtype=float) if intracell_hoppings is None else np.asarray(intracell_hoppings, dtype=float)
    t2_bonds = np.full(n_cells - 1, t2, dtype=float) if intercell_hoppings is None else np.asarray(intercell_hoppings, dtype=float)
    onsite = np.zeros(2 * n_cells, dtype=float) if onsite_energies is None else np.asarray(onsite_energies, dtype=float)
    if t1_bonds.shape != (n_cells,):
        raise ValueError(f"intracell_hoppings must have shape ({n_cells},)")
    if t2_bonds.shape != (n_cells - 1,):
        raise ValueError(f"intercell_hoppings must have shape ({n_cells - 1},)")
    if onsite.shape != (2 * n_cells,):
        raise ValueError(f"onsite_energies must have shape ({2 * n_cells},)")

    hamiltonian = np.diag(onsite)
    for cell in range(n_cells):
        a_site, b_site = 2 * cell, 2 * cell + 1
        hamiltonian[a_site, b_site] = t1_bonds[cell]
        hamiltonian[b_site, a_site] = t1_bonds[cell]
        if cell < n_cells - 1:
            next_a_site = 2 * (cell + 1)
            hamiltonian[b_site, next_a_site] = t2_bonds[cell]
            hamiltonian[next_a_site, b_site] = t2_bonds[cell]
    return hamiltonian


def chiral_symmetry_error(hamiltonian: np.ndarray) -> float:
    """Return max|Gamma H + H Gamma|; zero means chiral symmetry holds."""
    n_sites = hamiltonian.shape[0]
    if n_sites % 2:
        raise ValueError("SSH Hamiltonian must have an even number of sites")
    gamma_diagonal = np.tile([1.0, -1.0], n_sites // 2)
    # Gamma is diagonal, so multiplying by it scales rows/columns elementwise.
    anticommutator = gamma_diagonal[:, None] * hamiltonian + hamiltonian * gamma_diagonal[None, :]
    return float(np.max(np.abs(anticommutator)))


def edge_weights(states: np.ndarray, n_cells: int, edge_cells: int = 4) -> np.ndarray:
    """Probability in the first/last ``edge_cells`` unit cells for each state."""
    if edge_cells < 1 or edge_cells * 2 > n_cells:
        raise ValueError("edge_cells must be positive and no more than half the chain")
    n_edge_sites = 2 * edge_cells
    indices = np.r_[0:n_edge_sites, 2 * n_cells - n_edge_sites : 2 * n_cells]
    return np.sum(np.abs(states[indices, :]) ** 2, axis=0)


def iprs(states: np.ndarray) -> np.ndarray:
    """Inverse participation ratio: sum_i |psi_i|^4 for each eigenstate."""
    return np.sum(np.abs(states) ** 4, axis=0)


def analyze_spectrum(
    energies: np.ndarray,
    states: np.ndarray,
    n_cells: int,
    *,
    edge_cells: int = 4,
    energy_window: float = 0.5,
    minimum_edge_weight: float = 0.5,
) -> dict[str, float]:
    """Select a candidate edge-mode pair inside a central energy window.

    The pair is selected from states satisfying both ``|E| < energy_window``
    and ``P_edge >= minimum_edge_weight``, then ranked by edge weight. The
    energy window is set to the clean half-gap by the sweep. This is a
    transparent finite-chain diagnostic, not a topological invariant; inspect
    probability profiles when interpreting strong-disorder results.
    """
    weights = edge_weights(states, n_cells, edge_cells)
    participation = iprs(states)
    candidates = np.flatnonzero(
        (np.abs(energies) < energy_window) & (weights >= minimum_edge_weight)
    )
    if len(candidates) < 2:
        return {
            "edge_energy_1": float("nan"),
            "edge_energy_2": float("nan"),
            "edge_energy_pair_midpoint": float("nan"),
            "edge_energy_mean_abs": float("nan"),
            "edge_energy_splitting": float("nan"),
            "edge_weight_mean": float("nan"),
            "edge_weight_min": float("nan"),
            "edge_ipr_mean": float("nan"),
            "bulk_gap_estimate": float("nan"),
            "edge_modes_detected": 0.0,
        }
    edge_indices = candidates[np.argsort(weights[candidates])[-2:]]
    edge_indices = edge_indices[np.argsort(energies[edge_indices])]
    bulk_indices = np.setdiff1d(np.arange(len(energies)), edge_indices)
    bulk_energies = energies[bulk_indices]
    below_zero = bulk_energies[bulk_energies < 0]
    above_zero = bulk_energies[bulk_energies > 0]
    if len(below_zero) and len(above_zero):
        bulk_gap = float(np.min(above_zero) - np.max(below_zero))
    else:
        bulk_gap = float("nan")

    selected_energies = energies[edge_indices]
    selected_weights = weights[edge_indices]
    return {
        "edge_energy_1": float(selected_energies[0]),
        "edge_energy_2": float(selected_energies[1]),
        "edge_energy_pair_midpoint": float(np.mean(selected_energies)),
        "edge_energy_mean_abs": float(np.mean(np.abs(selected_energies))),
        "edge_energy_splitting": float(abs(selected_energies[1] - selected_energies[0])),
        "edge_weight_mean": float(np.mean(selected_weights)),
        "edge_weight_min": float(np.min(selected_weights)),
        "edge_ipr_mean": float(np.mean(participation[edge_indices])),
        "bulk_gap_estimate": bulk_gap,
        "edge_modes_detected": 1.0,
    }
