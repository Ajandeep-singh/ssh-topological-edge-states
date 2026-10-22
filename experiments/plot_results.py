"""Plot disorder-sweep summary data created by disorder_sweep.py."""

from pathlib import Path
import argparse
import os

import pandas as pd
import matplotlib

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".matplotlib"))
matplotlib.use("Agg")
import matplotlib.pyplot as plt


METRICS = [
    ("edge_weight_mean", "Average edge-region probability", "Edge weight"),
    ("edge_energy_pair_midpoint", "Midpoint of candidate edge-state energies", "Energy"),
    ("edge_energy_1", "Lower-energy candidate edge state", "Energy"),
    ("edge_energy_2", "Higher-energy candidate edge state", "Energy"),
    ("edge_energy_mean_abs", "Mean absolute edge-state energy", "Energy"),
    ("edge_energy_splitting", "Edge-state energy splitting", "Energy"),
    ("edge_ipr_mean", "Mean IPR of selected edge states", "IPR"),
    ("bulk_gap_estimate", "Estimated gap between non-edge states", "Energy"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=Path("results/data/disorder_summary.csv"))
    parser.add_argument("--figure-dir", type=Path, default=Path("results/figures"))
    args = parser.parse_args()
    data = pd.read_csv(args.summary)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    for metric, title, ylabel in METRICS:
        fig, axes = plt.subplots(1, len(sorted(data["n_cells"].unique())), figsize=(5 * len(sorted(data["n_cells"].unique())), 4), squeeze=False)
        for ax, n_cells in zip(axes[0], sorted(data["n_cells"].unique())):
            subset = data[data["n_cells"] == n_cells]
            for model, label in (("hopping", "Hopping disorder"), ("onsite", "On-site disorder")):
                curve = subset[subset["disorder_model"] == model].sort_values("W")
                ax.errorbar(
                    curve["W"],
                    curve[f"{metric}_mean"],
                    yerr=curve[f"{metric}_sem"],
                    marker="o",
                    markersize=3,
                    capsize=2,
                    label=label,
                )
            ax.set_title(f"N = {n_cells} unit cells")
            ax.set_xlabel("Disorder strength W")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.25)
        axes[0, 0].legend()
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(args.figure_dir / f"{metric}.png", dpi=160)
        plt.close(fig)

    metric = "edge_energy_splitting"
    fig, axes = plt.subplots(
        1,
        len(sorted(data["n_cells"].unique())),
        figsize=(5 * len(sorted(data["n_cells"].unique())), 4),
        squeeze=False,
    )
    for ax, n_cells in zip(axes[0], sorted(data["n_cells"].unique())):
        subset = data[data["n_cells"] == n_cells]
        for model, label in (("hopping", "Hopping disorder"), ("onsite", "On-site disorder")):
            curve = subset[subset["disorder_model"] == model].sort_values("W")
            ax.plot(curve["W"], curve[f"{metric}_median"], marker="o", markersize=3, label=label)
            ax.fill_between(
                curve["W"].to_numpy(dtype=float),
                curve[f"{metric}_q25"].to_numpy(dtype=float),
                curve[f"{metric}_q75"].to_numpy(dtype=float),
                alpha=0.18,
            )
        ax.set_title(f"N = {n_cells} unit cells")
        ax.set_xlabel("Disorder strength W")
        ax.set_ylabel("Median energy splitting")
        ax.grid(alpha=0.25)
    axes[0, 0].legend()
    fig.suptitle("Typical edge-state splitting (shaded area: interquartile range)")
    fig.tight_layout()
    fig.savefig(args.figure_dir / "edge_energy_splitting_median_iqr.png", dpi=160)
    plt.close(fig)
    print(f"Saved metric plots in {args.figure_dir}")


if __name__ == "__main__":
    main()
