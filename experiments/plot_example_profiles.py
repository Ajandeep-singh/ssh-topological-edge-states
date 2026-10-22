"""Save reproducible example wavefunction profiles for selected edge candidates."""

from pathlib import Path
import os
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.disorder import sample_hopping_disorder, sample_onsite_disorder
from src.ssh_model import build_ssh_hamiltonian, edge_weights


def main() -> None:
    n_cells, t1, t2 = 40, 0.5, 1.0
    strengths = (0.0, 0.25, 0.5)
    base_seed = 20261022
    sites = np.arange(2 * n_cells)
    profile_rows = []
    figure, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)

    for row, model in enumerate(("hopping", "onsite")):
        for column, strength in enumerate(strengths):
            sample_seed = base_seed + 10 * row + column
            rng = np.random.default_rng(sample_seed)
            if model == "hopping":
                t1_bonds, t2_bonds = sample_hopping_disorder(
                    n_cells, t1, t2, strength, rng
                )
                onsite = None
            else:
                t1_bonds = t2_bonds = None
                onsite = sample_onsite_disorder(n_cells, strength, rng)
            hamiltonian = build_ssh_hamiltonian(
                n_cells,
                t1,
                t2,
                intracell_hoppings=t1_bonds,
                intercell_hoppings=t2_bonds,
                onsite_energies=onsite,
            )
            energies, states = np.linalg.eigh(hamiltonian)
            weights = edge_weights(states, n_cells, edge_cells=4)
            candidates = np.flatnonzero(
                (np.abs(energies) < abs(t2 - t1)) & (weights >= 0.5)
            )
            if len(candidates) < 2:
                raise RuntimeError(
                    f"Could not identify two edge candidates for {model}, W={strength}"
                )
            selected = candidates[np.argsort(weights[candidates])[-2:]]
            selected = selected[np.argsort(energies[selected])]
            ax = axes[row, column]
            for rank, state_index in enumerate(selected, start=1):
                probability = np.abs(states[:, state_index]) ** 2
                ax.plot(
                    sites,
                    probability,
                    marker="o",
                    markersize=2,
                    linewidth=1,
                    label=f"E={energies[state_index]:+.3g}",
                )
                for site, value in enumerate(probability):
                    profile_rows.append(
                        {
                            "n_cells": n_cells,
                            "disorder_model": model,
                            "W": strength,
                            "sample_seed": sample_seed,
                            "candidate_rank": rank,
                            "energy": energies[state_index],
                            "site": site,
                            "probability": value,
                        }
                    )
            ax.set_title(f"{model.capitalize()} disorder, W={strength:g}")
            ax.set_xlabel("Site index")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
    axes[0, 0].set_ylabel("Probability")
    axes[1, 0].set_ylabel("Probability")
    figure.suptitle("Example candidate edge-state profiles (N=40, fixed sample seeds)")
    figure.tight_layout()

    figure_dir = PROJECT_ROOT / "results" / "figures"
    data_dir = PROJECT_ROOT / "results" / "data"
    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(figure_dir / "edge_state_profiles.png", dpi=160)
    plt.close(figure)
    pd.DataFrame(profile_rows).to_csv(
        data_dir / "example_edge_state_profiles.csv", index=False
    )
    print(f"Saved example profile plot to {figure_dir / 'edge_state_profiles.png'}")
    print(f"Saved profile values to {data_dir / 'example_edge_state_profiles.csv'}")
    print(f"Base sample seed: {base_seed}")


if __name__ == "__main__":
    main()
