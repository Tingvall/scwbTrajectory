# scwbTrajectory

This repository is the explicitly activated scworkbench
Trajectory operation pack. Its integrated ABI-v3 manifest retains the accepted
S1/S2/S4 and accepted Monocle3 operations and adds the Phase 5 S5 Palantir
fit/project methods alongside the S3 non-adoption Slingshot methods:
`slingshot`, `slingshot_fit`, `slingshot_lineage_weights`, and
`slingshot_project`. Call all methods under the `scwbtrajectory::` namespace after
explicit activation. `slingshot_adopt` remains a core operation; this pack
does not claim external-state adoption.

`trajectory_single_changepoint_regions` detects one change point per feature along
an ordered pseudotime axis. Its numeric core runs in a private Python process using
`ruptures`; the interpreter comes from `scwbTrajectory.changepoint_python`, else
`SCWORKBENCH_CHANGEPPOINT_PYTHON`, else `python3`, and the probe publishes the
interpreter path and the library versions it actually used in `run_summary`.
Published ranks and pseudotime boundaries are half-open: boundary rank `b` sits
before the cell at 0-based rank `b`, and `trend_density` publishes the count of
boundary `b` at cell rank `b`.

The accepted pack that the standard profile qualified has 25 operations at version
0.1.1. This branch carries 26 operations at 0.1.2 and is experimental; it is not
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
