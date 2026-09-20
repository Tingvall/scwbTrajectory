# scwbTrajectory

This repository is the provisional, explicitly activated scworkbench
Trajectory operation pack. Its Phase 5 S1 pilot owns four ABI-v3 operations:
`trajectory_branch_compare`, `trajectory_route_smooth`, `space_path_score`,
and `space_path_project_2d`.

Installation is inert. Activate the installed manifest explicitly with
`wb_operation_pack_activate()`. The pack consumes only authenticated typed
inputs and publishes complete declared output closures through core-owned
atomic staging. It does not receive a bundle, database, registry, SQL, or
committed-path handle.

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
