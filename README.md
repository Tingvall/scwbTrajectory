# scwbTrajectory

This repository is the explicitly activated scworkbench
Trajectory operation pack. Its integrated ABI-v3 manifest retains the accepted
S1/S2/S4 and accepted Monocle3 operations and adds the Phase 5 S5 Palantir
fit/project methods alongside the S3 non-adoption Slingshot methods:
`slingshot`, `slingshot_fit`, `slingshot_lineage_weights`, and
`slingshot_project`. `slingshot_fit` additionally exposes Slingshot's optional
`omega` and `omega_scale` controls; leaving them absent preserves Slingshot's
existing defaults. Call all methods under the `scwbtrajectory::` namespace after
explicit activation. `slingshot_adopt` remains a core operation; this pack
does not claim external-state adoption.

`trajectory_pseudotime_from_cell_vector` adapts one authenticated, promoted,
continuous scalar `cell_vector` into a one-key `cell_vector_family`. It copies
only the finite trajectory-cell rows without orienting, negating, normalizing,
rescaling, shifting, or imputing values, publishes exact input/finite/missing
counts, and retains the source artifact's authenticated ancestry. Cells outside
the trajectory are omitted from the scalar artifact and therefore remain absent
from the adapted family; plots may left-join the selected lineage to the full
cell universe to render those cells separately. Fitted lineage-family
pseudotime is the default trajectory source. This adapter is an explicit
fallback when fitted lineage output is unavailable or scalar pseudotime is
specifically requested; downstream consumers must still select its lineage key
explicitly.

`trajectory_occupancy` describes where support-qualified cells fall along one
explicit pseudotime key. It applies global equal-width bins, reports counts and
within-group fractions, optionally adds group-by-sample strata, and preserves
complete accepted, withheld, and accepted-but-missing accounting. The raw
pseudotime input must retain the complete analysis-cell axis; only cells marked
accepted with finite pseudotime enter bins. The operation is descriptive-only
and does not perform cell-level inferential comparisons.

The signal operations answer distinct questions. `trajectory_signal_dynamics`
describes whole-trajectory curves for precomputed cell-level signals.
`trajectory_geneset_scores` is the canonical reusable score-construction
operation. It adjusts each gene for `nCount_RNA` by default using an
authenticated `depth` cell vector, jointly fitting a natural pseudotime spline
and centered `log1p(depth)`, subtracting only the fitted depth component,
standardizing each adjusted gene, and averaging genes within each set. Cell-level scores and
gene-level diagnostics are retained as `adjusted_scores` and
`adjustment_diagnostics`. Callers must provide an explicit ordered
`feature_sets` vector; score keys are authenticated as
`pseudotime_key::feature_set` in that order. Those scores can be sent independently to
`trajectory_signal_dynamics_explicit`, `trajectory_signal_intervals`, or both.
The established `trajectory_geneset_dynamics` operations retain their combined
scoring-and-dynamics contracts and use the same scoring kernel. They also
require explicit ordered `feature_sets`; model-backed gene-set dynamics
requires one explicit `lineage` so its adjusted-score keys are fully granted.
`adjustment_spline_df` controls the pseudotime spline flexibility for all
depth-adjustable operations; it defaults to 4 and accepts integers from 2 to
10. The selected value is retained in output provenance.
`trajectory_signal_intervals` tests whether a signal has globally qualified
trajectory dynamics, performs local elevated and reduced tests, and merges
supported adjacent local tests into signal intervals. Its terminals are
`signal_curve`, `global_test`, `local_tests`, `signal_intervals`,
`adjustment_summary`, and `fit_summary`. Optional continuous technical
covariates are fitted without pseudotime. Depth adjustment is enabled by
default and uses an authenticated `depth` cell vector representing
`nCount_RNA`; a signal already marked with the same depth adjustment is reused
without a second adjustment. With sample identities,
`permutation_mode = "auto"` shuffles within samples; otherwise it records an
unblocked analysis. `trajectory_signal_regions` is the deprecated 0.1.4
compatibility identity and retains its historical terminal contract; it is not
a separate scientific method.

`trajectory_preranked_gsea` performs gene-set enrichment from an authenticated
expression matrix and the matrix-bound, checksummed feature-set registry
snapshot. `SCT_data` is recommended. Optional partial depth adjustment is
disabled by default because SCT expression is already depth-normalized; when
enabled, it requires an authenticated `depth` cell vector and defaults its
recorded field identity to `nCount_RNA`.
`analysis_scope = "global"` builds one full-trajectory ranking;
`"local"` builds rankings over a fixed grid of overlapping pseudotime windows;
and `"both"` runs both under the same lineage, orientation, cell selection,
covariate, sample-block, support, and feature-universe policies. It publishes
`global_enrichment`, `window_enrichment`, `enrichment_intervals`,
`leading_edge`, `ranked_features`, and `fit_summary`. Adjacent supported window
results merge only within the same feature set and enrichment direction. The
leading-edge terminal uses fgsea's actual leading-edge membership and reports
peak-window and interval-consensus membership. Ranking ties are never jittered;
the authenticated feature order resolves them deterministically. The operation
does not download feature-set collections at runtime.

All three methods accept generic primitives. No Bdev lineage, signature
collection, or study-specific threshold is embedded.

`trajectory_single_changepoint_regions` detects one change point per feature along
an ordered pseudotime axis. Its numeric core runs in a private Python process using
`ruptures`; the interpreter comes from `scwbTrajectory.changepoint_python`, else
`SCWORKBENCH_CHANGEPPOINT_PYTHON`, else `python3`, and the probe publishes the
interpreter path and the library versions it actually used in `run_summary`.
Published ranks and pseudotime boundaries are half-open: boundary rank `b` sits
before the cell at 0-based rank `b`, and `trend_density` publishes the count of
boundary `b` at cell rank `b`. A region is gated by the density of its own
direction: an up region by the up density, a down region by the down density. The
combined density is published as a diagnostic and gates nothing, because mixing the
directions lets a strong up population carry a lone down breakpoint over
`change_threshold`. A feature with missing matrix values follows `non_finite_policy`:
`drop` publishes it unanalysed with its `missing_value_count`, `error` rejects the
run and names the feature and its count of missing values.

The accepted pack that the standard profile qualified has 25 operations at version
0.1.1. This branch carries 32 operations at 0.1.9 and is experimental; it is not
what that profile qualifies. Model-backed
`trajectory_geneset_dynamics` and `trajectory_signal_dynamics` require `model`.
Their separate `_explicit` operations require typed pseudotime and
`pseudotime_key`; these are distinct identities, not a shortcut through the
model interface. Inspect the qualified live interface before execution.

Installation is inert. Activate the installed manifest explicitly with
`wb_operation_pack_activate()`. The pack consumes only authenticated typed
inputs and publishes complete declared output closures through core-owned
atomic staging. It does not receive a bundle, database, registry, SQL, or
committed-path handle.

`palantir_project` requires `pickle_execution_trust =
"operator_verified_source_v1"`. This is an explicit operator trust decision,
not a signature, sandbox, or proof that a hash-valid pickle is safe to execute.
The private Python process receives bounded typed inputs and a checked
stage-local model, never a bundle path or host model path.

The 2D path projection declares its source space explicitly as a typed input.
The old core surface inferred that source from `path$space_id`; callers should
pass the same space as `source_space`. This is the only S1 call-shape
difference and keeps state authority inside the ABI reader boundary.

All outputs use `exact_canonical_value` replay equivalence. Branch comparison
publishes `branch_de_table` with `branch_summary`; route smoothing publishes
`trend_table` with `fit_summary`; path scoring publishes `path_scores` with
`fit_summary`; projection publishes `curves`. The engine stages and audits each
declared output closure as one transition before any session output is visible.

The S1 ABI adds two generic validation declarations without adding a reader or
writer: `declared_keys_with_optional_prefix_v1` authenticates fixed family keys
after an optional string prefix, and `data_dependent_family_v1` admits a
runtime-sized path-family row axis while preserving kind and semantic checks.
Projection also uses the engine's generic persisted-2D-space normalization,
which converts the public `cell_id`/`x`/`y` table view to the authenticated
row-named coordinate matrix transported to the pack.
