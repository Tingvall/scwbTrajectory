# Phase 5 S5 Monocle3 candidate

This candidate moves `monocle3_fit` and `monocle3_project` to the
`scwbtrajectory::` operation pack. It does not include Palantir. The migrated
count remains 19/24 pending independent acceptance and Controller integration.

`tests/fixtures/trajectory_s5_monocle3_old_core_oracle.json` combines the
executed old-core fit oracle with the independently captured old-core project
curves. The fit comparison requires the exact terminal closure, 24 ordered
pseudotime cells and canonical content identity, 25 ordered relation rows,
their exact keys, cells, weights, and canonical content identity. The model
RDS bytes varied across old-core runs, so the model check uses its contained
tree identity, engine and type, source linkage, promotion, relocation, G1
reopen, and project execution. Project curves require exact lineage key and
point order and elementwise bound `abs(observed-reference) <= 1e-8 +
1e-8*abs(reference)` for x/y coordinates. The
oracle records that policy before the pack comparison.

The fit runner reads an authenticated `space`, declares the fitted axis,
computes Monocle3 pseudotime and principal-graph membership, declares sealed
relation rows, and writes the complete output group through core terminal and
G2 writers. The project runner reads a promoted, authenticated G1 model and
an authenticated 2D space before writing its curve terminal. Neither runner
receives a bundle handle or accepts a caller-provided model tree path.

Method-owned tests live in `tests/test_s5_monocle3_contract.R`. They cover
registered execution, exact oracle comparison, promotion, finalization,
relocation, reopened project, read-only recipient validation, exact name
ownership, invalid inputs, zero-row rejection, model tree and selected-axis
tampering, source-space tampering, missing-pack predispatch rejection, and
publication failure with retry. The core composite test exercises the real
installed pack boundary; generic G1/G2 and relation authority tests remain
in core.
