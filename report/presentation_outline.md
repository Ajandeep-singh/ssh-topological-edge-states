# Short presentation: SSH edge states under disorder

Use this as an 8-slide, 5–7 minute presentation. All plots and values below were generated from the saved-data pipeline. Keep the stated classifier caveats with any interpretation.

1. **Title and question** — Topological Edge States in the SSH Model Under Disorder. Ask how hopping and on-site disorder affect edge localization and finite-chain energy splitting.
2. **The SSH chain** — Draw the `A_n`, `B_n` sites and alternating `t1`, `t2` bonds. State the open-boundary convention and the clean topological regime `t2 > t1`.
3. **Symmetry and edge profile** — Show the chiral operator `Gamma = diag(+1,-1,...)`, the anticommutation relation for hopping-only disorder, and the clean decay `|a_n| ∝ |t1/t2|^(n-1)`.
4. **Two disorder models** — Compare random bond changes with random site potentials. Explain why the first preserves chiral symmetry in this model and the second generally breaks it.
5. **Numerical method** — Show the chain lengths, disorder grid, fixed seed, ensemble count, edge probability window, and IPR definition. Explain that the two highest edge-weight states are candidates requiring profile inspection.
6. **Clean validation** — Show `results/figures/clean_spectrum.png` and `clean_edge_states.png`. Report the 80 by 80 matrix, `8.88e-16` maximum spectral-pair error, edge weight `0.996094`, and all seven checks passed. Compare the theoretical localization length `xi=1/ln(2)=1.4427` cells.
7. **Disorder results** — Show `results/figures/edge_energy_splitting_median_iqr.png`, `edge_weight_mean.png`, and `edge_state_profiles.png`. At `N=40, W=0.50`, hopping disorder had median splitting `1.52e-16` and mean `0.0212`; on-site disorder had median `0.2403` and mean `0.2575`. Both models' candidate pairs had median edge weight above `0.99` at this setting. Explain the hopping-disorder skew and the finite-chain candidate rule.
8. **Conclusion and limitations** — The sampled data show much greater typical splitting under on-site disorder than hopping disorder, consistent with the symmetry difference. State finite-size, sampling, edge-tracking, and gap-estimation limitations. Say this is a computational study of an established topic, not a novelty claim.
