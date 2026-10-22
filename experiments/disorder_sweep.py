"""Run seeded SSH disorder ensembles and save realization/summary CSV files.

From the project root: python -m experiments.disorder_sweep
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.disorder import sample_hopping_disorder, sample_onsite_disorder
from src.ssh_model import analyze_spectrum, build_ssh_hamiltonian


def run_sweep(
    n_cells_values: list[int],
    strengths: np.ndarray,
    n_realizations: int,
    seed: int,
    t1: float = 0.5,
    t2: float = 1.0,
    edge_cells: int = 4,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run both disorder models; return per-realization and summary tables."""
    rng = np.random.default_rng(seed)
    rows = []
    for n_cells in n_cells_values:
        if 2 * edge_cells > n_cells:
            raise ValueError(f"N={n_cells} is too small for edge_cells={edge_cells}")
        for strength in strengths:
            for disorder_model in ("hopping", "onsite"):
                for realization in range(n_realizations):
                    if disorder_model == "hopping":
                        t1_bonds, t2_bonds = sample_hopping_disorder(
                            n_cells, t1, t2, float(strength), rng
                        )
                        onsite = None
                    else:
                        t1_bonds = t2_bonds = None
                        onsite = sample_onsite_disorder(n_cells, float(strength), rng)
                    hamiltonian = build_ssh_hamiltonian(
                        n_cells,
                        t1,
                        t2,
                        intracell_hoppings=t1_bonds,
                        intercell_hoppings=t2_bonds,
                        onsite_energies=onsite,
                    )
                    energies, states = np.linalg.eigh(hamiltonian)
                    metrics = analyze_spectrum(
                        energies,
                        states,
                        n_cells,
                        edge_cells=edge_cells,
                        energy_window=abs(t2 - t1),
                    )
                    rows.append(
                        {
                            "n_cells": n_cells,
                            "n_sites": 2 * n_cells,
                            "t1": t1,
                            "t2": t2,
                            "edge_cells": edge_cells,
                            "disorder_model": disorder_model,
                            "W": float(strength),
                            "realization": realization,
                            "seed": seed,
                            **metrics,
                        }
                    )

    raw = pd.DataFrame(rows)
    metric_columns = [
        "edge_energy_1",
        "edge_energy_2",
        "edge_energy_pair_midpoint",
        "edge_energy_mean_abs",
        "edge_energy_splitting",
        "edge_weight_mean",
        "edge_weight_min",
        "edge_ipr_mean",
        "bulk_gap_estimate",
        "edge_modes_detected",
    ]
    grouped = raw.groupby(["n_cells", "disorder_model", "W"], as_index=False)
    summary = grouped[metric_columns].agg(["mean", "std", "median", "count"]).reset_index()
    # Flatten pandas' multi-level aggregation column names for convenient CSV use.
    summary.columns = [
        "_".join(str(part) for part in column if part != "").rstrip("_")
        if isinstance(column, tuple)
        else str(column)
        for column in summary.columns
    ]
    for metric in metric_columns:
        summary[f"{metric}_sem"] = summary[f"{metric}_std"] / np.sqrt(summary[f"{metric}_count"])
    quartiles = raw.groupby(["n_cells", "disorder_model", "W"])[metric_columns].quantile(
        [0.25, 0.75]
    ).unstack(level=-1)
    quartiles.columns = [
        f"{metric}_q{int(quantile * 100)}" for metric, quantile in quartiles.columns
    ]
    summary = summary.merge(
        quartiles.reset_index(), on=["n_cells", "disorder_model", "W"]
    )
    return raw, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--realizations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--n-cells", type=int, nargs="+", default=[20, 40, 80])
    parser.add_argument("--w-max", type=float, default=0.50)
    parser.add_argument("--w-step", type=float, default=0.05)
    parser.add_argument("--edge-cells", type=int, default=4)
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "results" / "data")
    args = parser.parse_args()
    if args.realizations < 1 or args.w_step <= 0 or args.w_max < 0:
        parser.error("realizations and w-step must be positive; w-max cannot be negative")
    strengths = np.round(np.arange(0, args.w_max + args.w_step / 2, args.w_step), 10)
    raw, summary = run_sweep(
        args.n_cells, strengths, args.realizations, args.seed, edge_cells=args.edge_cells
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_dir / "disorder_realizations.csv"
    summary_path = args.output_dir / "disorder_summary.csv"
    metadata_path = args.output_dir / "run_metadata.json"
    raw.to_csv(raw_path, index=False)
    summary.to_csv(summary_path, index=False)
    import matplotlib

    metadata = {
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "matplotlib_version": matplotlib.__version__,
        "n_cells": args.n_cells,
        "n_realizations": args.realizations,
        "t1": 0.5,
        "t2": 1.0,
        "disorder_strengths": strengths.tolist(),
        "disorder_models": ["hopping", "onsite"],
        "seed": args.seed,
        "edge_cells": args.edge_cells,
        "edge_selection": "two highest edge weights among states with |E| < |t2-t1| and edge weight >= 0.5",
        "edge_selection_rate": float(raw["edge_modes_detected"].mean()),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Saved {len(raw)} realizations to {raw_path}")
    print(f"Saved {len(summary)} ensemble rows to {summary_path}")
    print(f"Saved run metadata to {metadata_path}")
    print(f"Base random seed: {args.seed}")


if __name__ == "__main__":
    main()
