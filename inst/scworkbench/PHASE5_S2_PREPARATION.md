# Phase 5 S2 preparation

Base: scwbTrajectory `9f3cf7c6e463cc332eb371ca370e635bcbc0131d`.
Old-core numerical oracle: scworkbench
`ba24635611756240f21f6817044883b832f12701`.

The eight cases in `tests/fixtures/trajectory_s2_core_reference.json` use a
24-cell matrix, two lineage keys, one feature set, and two signal keys. The
model cases use a synthetic `model_ref$pseudotime_table` to exercise method
equations without optional fit packages. They do **not** establish G1 model
reopen, model-tree integrity, or derived fitted-axis acceptance. Each terminal
has a frozen core canonical content identity, structural fingerprint, kind,
table type, and row count. `trajectory_s2_compare_outputs()` checks the entire
output closure, including empty driver terminals. The dynamics cases pin
bootstrap and permutation seeds and check the source kind separately.
The pack also owns an explicit-family selector and shared bin/aggregate
equations, with focused comparisons against the old-core method helpers.
One focused test executes the old-core operations to verify that the frozen
fixture can be reproduced. A separate test executes the two direct pack
runners against the same frozen comparator. At S2 integration, replace the
old-core caller with the namespaced pack runner; the test must not continue to
depend on removed core method functions.

| Operation | Prepared now | Remaining for S2 integration |
|---|---|
| `trajectory_de` | Model method oracle for `de_table` and `fit_summary`, feature chunks of two. | Pack runner using bounded matrix and feature-set grants; actual G1 reopened S3/S5 model; derived fitted-axis acceptance for subset fits. |
| `trajectory_smooth` | Model method oracle for `trend_table` and `fit_summary`, five bins and bounded feature chunks. | Pack runner and actual G1 reopened model; derived fitted-axis acceptance for subset fits. |
| `trajectory_signature_smooth` | Model method oracle for `signature_table` and `fit_summary` against an authenticated feature-set name. | Pack runner using registry snapshot and bounded matrix; actual G1 reopened model and subset-axis support. |
| `trajectory_geneset_dynamics` | Separate model and explicit-pseudotime six-terminal oracles. An unregistered ABI-v3-shaped runner executes explicit mode through typed pseudotime, optional membership, feature-set snapshot, and bounded matrix chunks. All six output identities match the frozen oracle. | Integration-owned registration and manifest; atomic publication audit; actual G1 model fixture and subset-axis support for model mode. |
| `trajectory_vector_smooth` | Model method oracle for `signal_table` and `fit_summary`, selected signal key. | Pack runner and actual G1 reopened model; derived fitted-axis acceptance for subset fits. |
| `trajectory_signal_dynamics` | Separate model and explicit-pseudotime four-terminal oracles. An unregistered ABI-v3-shaped runner executes explicit mode through typed signal/pseudotime inputs and optional membership. All four output identities match the frozen oracle. | Integration-owned registration and manifest; atomic publication audit; actual G1 model fixture and subset-axis support for model mode. |

The explicit-pseudotime inputs are complete-axis `cell_vector_family`
primitives. The direct runners consume only the accepted typed readers and
return declared primitives to terminal writers; they do not receive a bundle,
database, registry handle, or committed path. They remain deliberately
unregistered because shared registration and manifest files belong to the
integration owner. The focused writer-fault test checks that a failed staged
write propagates and leaves the durable bundle unchanged. A full ABI atomic
rollback test requires integration-owned registration. This preparation does
not remove any core operation or assert migration completion.

The explicit runner tests compare all ten terminal identities against the
frozen core oracle, reject absent or unknown pseudotime and method keys and
insufficient or wrong-kind membership, and repeat seeded bootstrap and
permutation calculations. Seed changes alter the expected stochastic fields;
the runners restore the caller's RNG state, including an initially absent
seed.

## Smallest generic dependency exposed by model-backed S2

G1 exists at the oracle core pin, but its self-linked reader recomputes
`fitted_cell_order_identity` from the **full** durable parent-space axis.
Reproduction: bind a model to a parent space with 24 cells, fit only 12 cells,
then reopen that model for any S2 consumer. The selected 12-cell identity and
the reader's recomputed 24-cell identity differ before method dispatch. This
is the same derived-fitted-axis issue identified for S3, not a new S2-specific
capability. The smallest remedy is the already authorized generic G1/G2
ordered-subset sidecar and sealed fitted-axis linkage, with exact reopen
verification. S2 must not pad the excluded cells, alter the durable parent
space, or bypass G1 checks. A positive fixture must reopen a sealed 12-of-24
model and run all model-backed consumers; an adversarial fixture must reject a
reordered or tampered fitted-axis sidecar before dispatch.

## Integration checks still owed

Once the generic derived-axis capability and fitted-model producer are
accepted, replace the synthetic model oracle with a retained S3 or S5 model,
reopen it, and compare the six pack terminal closures. Add pack promotion,
rollback, replay-plan, pack-free inspection, and missing-pack dispatch checks
at the grouped S2 gate. Preserve each terminal's `exact_canonical_value`
policy; the frozen method oracle alone does not prove operation-level or
biological equivalence.
