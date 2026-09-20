# Phase 5 S1 operation contract

`scwbTrajectory` is explicitly activated and owns these ABI-v3 names:

| operation | typed inputs | parameters | terminal closure |
|---|---|---|---|
| `scwbtrajectory::trajectory_branch_compare` | `matrix`, `cell_vector_family` pseudotime, `relation_family` membership | branch keys, optional feature-set/features, `min_pct`, adjustment method, chunk size | `branch_de_table`, `branch_summary` |
| `scwbtrajectory::trajectory_route_smooth` | `matrix`, `cell_vector` pseudotime, `relation_family` membership | route key, optional feature-set/features, bins, summary, minimum bin size, weight policy, probability floor, chunk size | `trend_table`, `fit_summary` |
| `scwbtrajectory::space_path_score` | `space`, `space_path_family` | orthogonal projection, Euclidean distance, optional score-key prefix | `path_scores`, `fit_summary` |
| `scwbtrajectory::space_path_project_2d` | `space_path_family`, explicit source `space`, target `space_2d` | local interpolation, optional point count | `curves` |

Every terminal declares `exact_canonical_value`. Matrix operations consume an
authenticated ordered feature partition and the full effective cell axis;
feature selection occurs only after all granted chunks have been read. Optional
feature-set selection consumes an authenticated registry snapshot. Space and
path operations reject source-space mismatches, insufficient shared cells,
non-finite coordinates, malformed family rows, and non-2D targets before the
transition is published.

The old core projection inferred the source space from `path$space_id` and
returned source/path identifiers in its ephemeral parameter list. The pack
requires `source_space` explicitly. ABI-v3 currently retains `project_mode` and
`approx_points`; source, plot, and path identities remain durably available as
authenticated typed input descriptors and provenance bindings rather than
duplicated output metadata. This is the bounded S1 call/parameter difference;
scientific output identities are byte-for-byte equivalent to the frozen core
fixture.
