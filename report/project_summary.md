# Project summary

## Topological Edge States in the SSH Model Under Disorder

Built a reproducible Python exact-diagonalization study of the finite, open-boundary Su-Schrieffer-Heeger chain. The project compares independent uniform hopping disorder, which preserves chiral symmetry in the nearest-neighbor model, with independent uniform on-site disorder, which generally breaks it.

The clean `N=40`, `t1=0.5`, `t2=1.0` chain passed all seven implemented checks: Hamiltonian size and Hermiticity, chiral symmetry, positive/negative spectral pairing, two near-zero states, eigenvector normalization, and edge localization. A seeded ensemble then covered `N=20, 40, 80`, 11 disorder amplitudes from 0 to 0.50, and 100 realizations for each setting and disorder model (6,600 realizations total). The repository includes per-realization data, ensemble summaries with mean/SEM and median/quartiles, figures, metadata, and example wavefunction profiles.

At `N=40` and `W=0.50`, the candidate-pair median splitting was `1.52e-16` for hopping disorder and `0.2403` for on-site disorder. The hopping-disorder mean was `0.0212`, revealing rare larger values in a strongly skewed distribution. The median selected-pair edge weights exceeded `0.99` for both models at this point.

These findings are finite-chain numerical observations under a documented candidate-state rule; they do not establish a thermodynamic transition, a topological invariant, or novelty. The broad comparison has prior literature. The project is best described as a reproducible computational study demonstrating how symmetry, localization diagnostics, state tracking, finite-size effects, and ensemble statistics work together.

## Short CV description

Implemented a seeded NumPy exact-diagonalization study of 6,600 disordered finite SSH chains, comparing chiral-symmetric hopping disorder with chiral-symmetry-breaking on-site disorder. Validated the clean spectrum and edge-localized wavefunctions; analyzed edge weight, IPR, finite-size energy splitting, and gap proxies with ensemble uncertainty and reproducible plots.
