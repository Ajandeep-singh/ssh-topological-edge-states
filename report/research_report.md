# Topological Edge States in the SSH Model Under Disorder

## Abstract

We numerically diagonalized finite open SSH chains with chiral-symmetry-preserving random hopping and chiral-symmetry-breaking random on-site potentials. The clean `N=40` chain passed matrix, spectral-pairing, normalization, and edge-localization checks. A seeded ensemble of 6,600 realizations covered `N=20, 40, 80`, 11 disorder amplitudes from 0 to 0.50, and 100 realizations per setting. Under the stated central-energy and edge-weight selection rule, typical on-site-disorder edge splitting increased with disorder, while hopping-disorder splitting was usually near numerical zero with a skewed tail of larger values. These results describe finite chains and heuristic candidate-state tracking; they do not locate a thermodynamic transition or establish novelty.

## 1. Introduction and research question

How do chiral-symmetry-preserving hopping disorder and chiral-symmetry-breaking on-site disorder affect edge-state localization and energy splitting in finite SSH chains?

The SSH chain is a canonical one-dimensional model with alternating nearest-neighbor hoppings. The original Su, Schrieffer, and Heeger work developed the dimerized-chain description in the context of solitons in polyacetylene [1]. Random hopping in one-dimensional nearest-neighbor chains has also been studied as off-diagonal disorder with sublattice symmetry [2]. Later work has treated disordered SSH chains with additional hopping ranges and transport observables [3], and direct work on the disordered SSH model with both on-site and hopping disorder has appeared [4].

This project is therefore framed as a reproducible finite-chain numerical study and learning project, not as a novel discovery. Its scope is the comparison of two disorder ensembles in the real, nearest-neighbor, single-particle open chain, using finite-size spectra and wavefunction diagnostics. It does not calculate a real-space winding invariant or Lyapunov exponent, and its finite-chain gap estimate should not be treated as a thermodynamic phase boundary.

**Research question:** How do chiral-symmetry-preserving hopping disorder and chiral-symmetry-breaking on-site disorder affect edge-state localization and energy splitting in finite SSH chains?

## 2. Model and symmetry

With site order `A1, B1, A2, B2, ...`, the open-chain Hamiltonian is

`H = sum_n t1,n (|A_n><B_n| + h.c.) + sum_n t2,n (|B_n><A_(n+1)| + h.c.) + sum_i epsilon_i |i><i|`,

where the intercell sum ends at `N-1`. The clean chain has `epsilon_i=0`, uniform `t1,n=t1`, and uniform `t2,n=t2`. Its Bloch bands are `E±(k) = ±sqrt(t1^2 + t2^2 + 2 t1 t2 cos(k))`. For this termination and an open gap, `t2>t1>0` is the topological phase. A semi-infinite left edge mode lives on one sublattice with amplitude proportional to `(-t1/t2)^(n-1)` and localization length `xi=1/ln(t2/t1)` in unit cells.

The chiral operator is `Gamma=diag(+1,-1,+1,-1,...)`. With only nearest-neighbor inter-sublattice hopping, `Gamma H Gamma = -H`, including for spatially varying real hopping amplitudes. On-site terms commute with `Gamma`, so generic nonzero `epsilon_i` break this symmetry. Chiral symmetry explains the positive/negative spectral pairing in the hopping-disorder model; it does not alone guarantee an identifiable isolated edge mode at every disorder strength.

## 3. Numerical methods

The starting parameters were `t1=0.5`, `t2=1.0`, open boundaries, `N=20, 40, 80`, and uniform disorder amplitudes `W=0, 0.05, ..., 0.50`. The ensemble used 100 realizations per setting and base seed `20261021`. It records each realization and summarizes mean, sample standard deviation, median, quartiles, count, and standard error of the mean.

Candidate states must satisfy `|E| < |t2-t1| = 0.5` and have probability at least 0.5 in the first and last four unit cells combined. The two qualifying states with the largest edge weights are selected. If fewer than two qualify, pair metrics are missing. This is a finite-chain heuristic, not an invariant; inspect spatial profiles, particularly at stronger disorder. The reported bulk-gap estimate removes those two candidates and measures the separation between the nearest remaining negative and positive eigenvalues. This estimate can be contaminated by bulk in-gap states and is not the edge-state splitting.

**Execution record:** Python 3.12.10, NumPy 2.5.4, Matplotlib 3.11.2, pandas 3.0.6. The clean validation passed all seven checks. The full sweep completed with 6,600 realization rows and 66 summary rows; no summary metrics were missing. The candidate-pair rule selected a pair in all 6,600 realizations. Exact parameters and software versions are recorded in `results/data/run_metadata.json`.

## 4. Results

### 4.1 Clean-chain validation

For `N=40`, `t1=0.5`, and `t2=1.0`, the Hamiltonian was 80 by 80 and symmetric. The maximum chiral anticommutator error was 0, the maximum positive/negative eigenvalue-pair error was `8.88e-16`, and the maximum eigenvector normalization error was `1.33e-15`. The near-zero energies were approximately `-6.8205e-13` and `+6.8219e-13`. Each had edge-region probability `0.996094` in the first and last four unit cells and IPR `0.300000`. The analytical bulk gap is 1, and the edge localization length is `1/ln(2) = 1.4427` unit cells.

### 4.2 Disorder ensembles

At `N=40`, the median edge-state splitting for hopping disorder stayed at numerical zero over most of the sweep. At `W=0.50`, its median was `1.52e-16`, while its mean was `0.0212` and sample standard deviation was `0.1099`, indicating a strongly right-skewed distribution with rare larger splittings. For on-site disorder, the `N=40` median splitting increased from `0.0356` at `W=0.10` to `0.1191` at `W=0.25` and `0.2403` at `W=0.50`. At `W=0.50`, the mean on-site splitting was `0.2575` with standard deviation `0.1742` across 100 realizations. The pair-midpoint mean was `-0.00693`; individual candidate energies were no longer constrained to opposite signs.

At `N=40`, `W=0.50`, the median selected-pair edge weight was `0.9961` for hopping disorder and `0.9936` for on-site disorder. Median IPRs were `0.5372` and `0.5345`, respectively. Median non-edge gap estimates were `0.4937` and `0.7398`. These gap estimates are finite-spectrum proxies after removing the candidate pair, not thermodynamic bulk gaps.

For the clean `N=80` chain, the splitting is below useful double-precision resolution; the computed value should be read as zero at machine precision, not as an accurate physical splitting. The hopping-disorder splitting distribution at high W also has near-zero medians but some larger outliers, so the mean alone is misleading. The report therefore includes mean/SEM and median/interquartile-range plots.

Plots include the clean spectrum and wavefunctions, example disorder profiles, average edge weight, signed candidate energies and their midpoint, mean absolute edge energy, edge splitting, IPR, and finite-size comparisons. Raw observations and summaries are in CSV; parameters and software versions are in `results/data/run_metadata.json`.

## 5. Discussion

The results are consistent with the symmetry-based expectation that hopping disorder retains chiral spectral pairing while on-site potentials break that constraint and move edge energies away from symmetric values. In this finite ensemble, on-site disorder produced a clear increase in typical edge-pair splitting as W rose, whereas typical hopping-disorder splitting stayed much smaller. The hopping-disorder mean at high W was affected by rare realizations, which is why medians and interquartile ranges are included.

High edge probability and IPR do not by themselves establish topological protection. The example profiles show that the selected states remain concentrated near the ends for the displayed parameter points, but selection uses a fixed clean-gap window and a finite edge region. The finite-size gap proxy also responds to disorder-localized in-gap states. These data do not identify a phase boundary or establish a disorder-driven topological transition.

## 6. Limitations

- The calculation is a finite, single-particle, nearest-neighbor model with independent uniform disorder; it omits interactions and experimental imperfections.
- Candidate edge modes are selected by a fixed clean half-gap window and a finite edge-region weight threshold. This is not a bulk invariant. All realizations met the rule here, but a classifier can still mistake a disorder-localized in-gap state for a topological boundary state.
- The sample profiles illustrate selected candidates for fixed seeds; they do not show every realization. A more complete robustness claim should test sensitivity to the edge-region width and energy window and should inspect distributions of profiles or overlaps.
- IPR measures concentration but cannot distinguish edge localization from localization in the bulk.
- The finite-spectrum gap proxy can be affected by localized in-gap states and is not a thermodynamic bulk gap or a phase-transition criterion.
- At `N=80`, clean edge splitting is below useful double-precision resolution. Treat its computed value as zero at machine precision.
- The ensemble size is 100 per setting. Rare events make mean splitting particularly sensitive to outliers; medians and interquartile ranges are shown alongside mean and SEM.
- The project is a computational study of an established model. It makes no novelty claim.

## 7. Conclusion

The clean finite SSH chain reproduces the expected midgap, end-localized pair. In the sampled disorder range, selected edge candidates remained strongly concentrated in the predefined edge regions. On-site disorder shifted and split the candidate pair more strongly as `W` increased, while hopping disorder generally retained near-zero paired energies; the hopping-disorder mean splitting at high W was sensitive to rare realizations. These observations are consistent with the contrasting symmetry constraints of the two Hamiltonians, but the finite-chain diagnostics do not establish a phase boundary or prove topological protection under all sampled disorder. The result is best presented as a reproducible numerical study with explicit state-tracking limits.

## References

[1] W. P. Su, J. R. Schrieffer, and A. J. Heeger, “Solitons in Polyacetylene,” *Physical Review Letters* **42**, 1698–1701 (1979), https://doi.org/10.1103/PhysRevLett.42.1698.

[2] A. Krishna and R. N. Bhatt, “Beyond universal behavior in the one-dimensional chain with random nearest-neighbor hopping,” *Physical Review B* **101**, 224203 (2020), https://doi.org/10.1103/PhysRevB.101.224203. This is a broader random-hopping chain study; it provides context for chiral/off-diagonal disorder, not a direct validation of the exact observables used here.

[3] B. Pérez-González, M. Bello, Á. Gómez-León, and G. Platero, “Interplay between long-range hopping and disorder in topological systems,” *Physical Review B* **99**, 035146 (2019), https://doi.org/10.1103/PhysRevB.99.035146. This extends the SSH setting with long-range hopping, so its model differs from the present nearest-neighbor chain.

[4] M. Hilke, “Disordered Su-Schrieffer-Heeger model,” *Physical Review B* **114**, 245111 (2026), https://doi.org/10.1103/by8p-dn9w. This recent work treats both on-site and hopping disorder and studies Lyapunov exponents and a real-space winding number; those are outside the present project’s numerical scope.
