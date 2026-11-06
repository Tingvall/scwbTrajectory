# Trajectory signal and enrichment operations

## Method distinction

`scwbtrajectory::trajectory_signal_dynamics` is descriptive. It summarizes a
precomputed cell-level signal over a whole selected trajectory and publishes
curve-oriented outputs. It does not qualify local intervals.

`scwbtrajectory::trajectory_signal_intervals` is inferential. Its `global_test`
is the permutation-qualified gate for its `local_tests`. Supported adjacent
local tests are merged into `signal_intervals`; the global test is not a second
global signal-analysis method. `trajectory_signal_regions` is the deprecated
0.1.4 compatibility identity for this same scientific method.

Both interval identities enable partial depth adjustment by default. A raw
signal must supply an authenticated `depth` cell vector; `depth_field` defaults
to `nCount_RNA`. `adjustment_spline_df` defaults to 4 and accepts integers from
2 to 10. The operation fits `signal ~ ns(pseudotime, df=adjustment_spline_df) +
centered_log1p(depth)` and removes only the depth term. Signals already marked
as depth-adjusted with the same field are not adjusted again.

## Reusable gene-set scores

`scwbtrajectory::trajectory_geneset_scores` performs gene-level partial depth
adjustment and score construction without requiring a whole-trajectory
dynamics analysis. It publishes `adjusted_scores`, `adjustment_diagnostics`,
and `fit_summary`. The score family retains its matrix, pseudotime, membership,
and depth identities and can be supplied independently to
`trajectory_signal_dynamics_explicit`, `trajectory_signal_intervals`, or both.
Intervals recognize the adjustment metadata and do not adjust these scores a
second time. `feature_sets` is required, ordered, unique, and nonempty; the
published family uses `pseudotime_key::feature_set` keys in parameter order and
effective cells in authenticated input order.
The combined gene-set dynamics operations use the same authenticated key rule;
their registered model-backed form additionally requires one explicit lineage.

`scwbtrajectory::trajectory_preranked_gsea` ranks genes from the authenticated
matrix feature universe. Its global scope produces one full-trajectory ranking.
Its local scope produces one complete ranking for each supported window on a
grid fixed before feature sets are examined. It applies preranked GSEA to each
ranking and merges adjacent supported window results into
`enrichment_intervals` by feature set and enrichment direction.

## Signal intervals

The canonical terminal closure is:

- `signal_curve`: windowed signal values and support.
- `global_test`: whole-trajectory permutation test.
- `local_tests`: unmerged elevated or reduced tests.
- `signal_intervals`: supported adjacent local tests merged into intervals.
- `adjustment_summary`: continuous-covariate adjustment diagnostics.
- `fit_summary`: lineage, support, permutation, parameter, and provenance
  summary.

Automatic permutation mode is blocked when `sample_id` is supplied and
unblocked otherwise. Covariate adjustment never includes pseudotime in the
technical-covariate design.

## Preranked GSEA

Use the authenticated `SCT_data` matrix for GSEA. `adjust_depth` defaults to
`FALSE`; when enabled with an authenticated `depth` cell vector, the same joint
pseudotime/depth model removes only the depth term before ranking.

`analysis_scope` is `global`, `local`, or `both`; `both` is the default.
`monotonic_association` uses a signed Spearman association, while
`early_late_contrast` uses the adjusted late mean minus early mean over explicit
trajectory fractions. Positive ranking statistics always indicate increasing
expression in the selected orientation.

Local raw effects are adjusted window means minus the adjusted trajectory
background outside the window. Ranking statistics are standardized against the
declared blocked or unblocked pseudotime-permutation null. The same cell axis,
sample blocks, covariate design, orientation, and support thresholds apply to
all genes.

`q_within_ranking` is adjusted within a ranking. `q_global_family` is adjusted
over the complete tested local feature-set-by-window family. These are declared
`p.adjust` corrections and are not described as maximum-statistic correction.

The terminal closure is:

- `global_enrichment`: full-trajectory GSEA results.
- `window_enrichment`: feature-set results in supported fixed windows.
- `enrichment_intervals`: supported adjacent windows merged by set and
  direction.
- `leading_edge`: true fgsea leading-edge feature membership.
- `ranked_features`: included and explicitly excluded matrix features for every
  ranking.
- `fit_summary`: engine, universe, set, support, permutation, seed, and tie
  diagnostics.

The operation requires fgsea 1.34.2, executes it serially with explicit seeds,
and preserves caller RNG state. It consumes only the authenticated matrix-bound
feature-set registry snapshot and performs no runtime collection downloads.
Original ranking values are never jittered. Ties use authenticated matrix
feature order, and leading-edge boundaries crossing a tied group are flagged.

Ranking qualification uses one deterministic permutation-index plan per
declared seed. Global monotonic statistics reuse centered gene ranks, while
endpoint and local-window effects use contrast-weight matrices. Null moments
are accumulated online from batched matrix products, and registered execution
streams every authenticated matrix feature in bounded chunks without reducing
the feature universe.
