"""Build and validate a clean, open-boundary SSH chain.

Run from the project root with: python src/clean_system.py
The plots are written to results/figures/ so they can be inspected later.
"""

from pathlib import Path
import os

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def build_ssh_hamiltonian(n_cells: int, t1: float, t2: float) -> np.ndarray:
    """Return the real-space Hamiltonian for a 2*n_cells-site SSH chain.

    Site order is A1, B1, A2, B2, ... . The alternating nearest-neighbor
    bonds are t1 (within a cell) and t2 (between cells).
    """
    hamiltonian = np.zeros((2 * n_cells, 2 * n_cells), dtype=float)

    for site in range(2 * n_cells - 1):
        hopping = t1 if site % 2 == 0 else t2
        hamiltonian[site, site + 1] = hopping
        hamiltonian[site + 1, site] = hopping  # Hermitian conjugate

    return hamiltonian


def main() -> None:
    # The hopping amplitudes set the energy scale; all energies are in units
    # of t2 here.
    n_cells = 40
    t1 = 0.5
    t2 = 1.0
    hamiltonian = build_ssh_hamiltonian(n_cells, t1, t2)

    # eigh is designed for real symmetric / complex Hermitian matrices. It
    # returns eigenvalues in ascending order and normalized eigenvectors in
    # columns of `states`.
    energies, states = np.linalg.eigh(hamiltonian)

    # For this site ordering, chiral symmetry is +1 on A and -1 on B.
    chiral_operator = np.diag([1.0, -1.0] * n_cells)
    chiral_symmetry_error = np.max(
        np.abs(chiral_operator @ hamiltonian + hamiltonian @ chiral_operator)
    )
    pair_error = np.max(np.abs(energies + energies[::-1]))

    # Treat the first/last four unit cells (eight sites at each end) as edge
    # regions. This captures most of the clean edge-state tails for t1/t2=0.5.
    edge_site_count = 2 * 4
    edge_indices = np.r_[0:edge_site_count, 2 * n_cells - edge_site_count : 2 * n_cells]
    closest_indices = np.argsort(np.abs(energies))[:2]
    edge_weights = [
        float(np.sum(np.abs(states[edge_indices, index]) ** 2))
        for index in closest_indices
    ]
    iprs = [float(np.sum(np.abs(states[:, index]) ** 4)) for index in closest_indices]
    normalization_error = float(
        np.max(np.abs(np.sum(np.abs(states) ** 2, axis=0) - 1.0))
    )

    # Infinite-chain bulk gap is the gap between the two bands. For t2 > t1,
    # the band edges are at +/- (t2-t1), so the full gap is 2*(t2-t1).
    bulk_gap = 2.0 * abs(t2 - t1)
    localization_length = 1.0 / np.log(t2 / t1)

    print(f"Hamiltonian shape: {hamiltonian.shape}")
    print(f"Hamiltonian is Hermitian: {np.allclose(hamiltonian, hamiltonian.T)}")
    print(f"Maximum chiral-symmetry error: {chiral_symmetry_error:.3e}")
    print(f"Maximum +/- eigenvalue-pair error: {pair_error:.3e}")
    print(f"Maximum eigenvector normalization error: {normalization_error:.3e}")
    print(f"Expected clean bulk gap: {bulk_gap:.6f}")
    print(f"Expected edge localization length: {localization_length:.4f} unit cells")
    print("Lowest five energies:", energies[:5])
    print("Highest five energies:", energies[-5:])
    for index, edge_weight, ipr in zip(closest_indices, edge_weights, iprs):
        print(
            f"Near-zero state: E={energies[index]:+.6e}, "
            f"edge weight={edge_weight:.6f}, IPR={ipr:.6f}"
        )

    checks = {
        "80-site Hamiltonian": hamiltonian.shape == (80, 80),
        "Hermitian Hamiltonian": np.allclose(hamiltonian, hamiltonian.T),
        "Chiral symmetry": chiral_symmetry_error < 1e-12,
        "Paired spectrum": pair_error < 1e-10,
        "Two near-zero states": np.max(np.abs(energies[closest_indices])) < 1e-8,
        "Normalized eigenvectors": normalization_error < 1e-12,
        "Edge-localized near-zero states": min(edge_weights) > 0.90,
    }
    print("\nClean-model validation:")
    for name, passed in checks.items():
        print(f"  {'PASS' if passed else 'FAIL'}  {name}")
    if not all(checks.values()):
        raise SystemExit("Clean-model validation failed; inspect diagnostics before continuing.")
    print("  PASS  All clean-model checks")

    # Save figures instead of requiring an interactive plotting window.
    figure_dir = Path("results") / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7, 4))
    plt.scatter(np.arange(len(energies)), energies, s=15)
    plt.axhline(0, color="black", linestyle="--", linewidth=1)
    plt.xlabel("Eigenstate index")
    plt.ylabel("Energy (in units of the hopping scale)")
    plt.title("Clean SSH energy spectrum")
    plt.tight_layout()
    plt.savefig(figure_dir / "clean_spectrum.png", dpi=160)
    plt.close()

    sites = np.arange(2 * n_cells)
    plt.figure(figsize=(8, 4))
    for index in closest_indices:
        probability = np.abs(states[:, index]) ** 2
        plt.plot(sites, probability, "o-", markersize=3, label=f"E = {energies[index]:+.2e}")
    plt.xlabel("Site index (A1, B1, A2, B2, ...)")
    plt.ylabel("Probability")
    plt.title("Two eigenstates closest to zero energy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure_dir / "clean_edge_states.png", dpi=160)
    plt.close()
    print(f"Plots saved in: {figure_dir.resolve()}")


if __name__ == "__main__":
    main()
