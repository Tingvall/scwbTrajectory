# scwbTrajectory 0.1.7

- Adds `trajectory_geneset_scores` as the canonical reusable gene-set scoring
  operation. It performs the same per-gene partial depth adjustment and
  standardization as gene-set dynamics without running a dynamics analysis.
- The resulting `adjusted_scores` can be consumed directly by generic signal
  dynamics or signal intervals. Adjustment identities remain attached so
  intervals do not remove depth twice.
- Refactors both established gene-set dynamics operations to use the same
  internal scoring kernel while retaining their public contracts.
- Exposes `adjustment_spline_df` across gene-set scoring, gene-set dynamics,
  signal intervals, and optional GSEA depth adjustment. The reproducible
  default remains 4; accepted values are integers from 2 through 10.

# scwbTrajectory 0.1.6

- Adds partial depth adjustment to both gene-set dynamics interfaces, enabled
  by default with `depth_field = "nCount_RNA"`. Each gene retains its fitted
  pseudotime spline while only the centered log-depth component is removed;
  adjusted genes are standardized before set averaging.
- Publishes cell-level adjusted gene-set scores and gene-level adjustment
  diagnostics, and uses the adjusted values for curves, dynamics, drivers, and
  permutation analyses.
- Enables partial depth adjustment by default for signal intervals and the
  deprecated signal-regions alias. Already-adjusted signals are recognized and
  not adjusted twice.
- Adds optional GSEA depth adjustment, disabled by default for the recommended
  authenticated `SCT_data` input. Raw-count UCell remains an independent
  comparison rather than the primary depth-controlled signal.

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
