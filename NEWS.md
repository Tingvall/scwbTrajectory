# scwbTrajectory 0.1.5

- Optimizes trajectory preranked GSEA by reusing deterministic permutation
  plans, caching gene ranks, evaluating global and local nulls with batched
  matrix multiplication, and streaming the authenticated feature universe in
  bounded chunks. The optimized kernels preserve the declared permutation
  assignments, caller RNG state, exclusions, tie policy, and output contract.

- Adds `scwbtrajectory::trajectory_signal_intervals` as the canonical name for
  permutation-qualified elevated and reduced signal intervals. The operation
  separates the global inferential gate, local tests, and merged signal
  intervals in its terminal names.
- Retains `scwbtrajectory::trajectory_signal_regions` as a deprecated
  compatibility identity with its 0.1.4 terminal contract.
- Adds `scwbtrajectory::trajectory_preranked_gsea` for global and fixed-window
  preranked GSEA over authenticated matrix and feature-set registry inputs,
  including true leading edges and merged enrichment intervals.
- Pins the GSEA engine to fgsea 1.34.2 and records its complete numerical,
  feature-universe, support, permutation, and tie configuration.
