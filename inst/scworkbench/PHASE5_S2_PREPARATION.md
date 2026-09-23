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
| `trajectory_de` | Direct Slingshot and Monocle3 model runners write `de_table` and `fit_summary` from G1 pseudotime and a typed matrix. Frozen Slingshot and both Monocle3 branch oracles match. | Register G1 and matrix grants; grouped publication audit. |
| `trajectory_smooth` | Direct Slingshot and Monocle3 model runners write `trend_table` and `fit_summary`. Frozen Slingshot and both Monocle3 branch oracles match. This runner also consumes a real relocated G1 private read. | Register G1 and matrix grants; grouped publication audit. |
| `trajectory_signature_smooth` | Direct Slingshot and Monocle3 model runners use a typed feature-set snapshot and bounded matrix. Frozen Slingshot and both Monocle3 branch oracles match. | Register G1, registry, and matrix grants; grouped publication audit. |
| `trajectory_geneset_dynamics` | Direct explicit, Slingshot, and Monocle3 model runners write six terminals. Frozen explicit and Slingshot identities and both Monocle3 branch oracles match. | Register source modes and grants; grouped publication audit. |
| `trajectory_vector_smooth` | Direct Slingshot and Monocle3 model runners use a typed signal family. Frozen Slingshot and both Monocle3 branch oracles match. | Register signal and G1 grants; grouped publication audit. |
| `trajectory_signal_dynamics` | Direct explicit, Slingshot, and Monocle3 model runners write four terminals. Frozen explicit and Slingshot identities and both Monocle3 branch oracles match. | Register source modes and grants; grouped publication audit. |

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

### Authenticated typed-reader boundary at the core pin

At scworkbench `ba24635611756240f21f6817044883b832f12701`, the
`wb_session_run_op()` path resolves required and present optional inputs and
checks their declared kinds in `R/session_context.R:434-456` and `:308`.
Before invoking any operation runner, it calls `.wb_typed_effective_cells()`
at `R/session_context.R:458-479`. That function validates each complete
resolved input **before** intersecting it with the requested session cells
(`R/getters_coverage.R:231-259`). For `cell_vector_family` pseudotime and
signal, it rejects invalid or foreign bundle cell IDs and repeated key/cell
pairs (`R/getters_coverage.R:131-139,204-208`). For `cell_vector` membership,
it rejects foreign and duplicate cell IDs (`:201-203`). The parent axis for
these primitive kinds is the canonical bundle cell universe from
`wb_cells(bundle, "all")` at `:172-173`; it is not the effective session
subset. A row inside that bundle axis but outside the requested subset is
valid source data and is filtered during binding.

The ABI-v3 runner binds these primitives through the v2 binder at
`R/operation_pack_v3_transition.R:3009-3012`. The binder filters membership
while retaining its source row order, and orders family rows by effective
cell order with original row order as a stable tie break
(`R/operation_pack_v2_transition.R:121-125,155-162`). The bound payload's
kind, content, and cell-order identities enter the input descriptor at
`R/operation_pack_v2_transition.R:224-245`; the manifest closure and kinds
are checked again at `:1356-1364`. The native primitive reader returns a
deep copy of that bound input at `:808-816`. These identity checks bind the
delivered payload and its order; a direct primitive supplied by a caller is
not a durable artifact with an independently verified source identity.

`tests/test_s2_explicit_execution_contract.R` exercises the core coverage
function against foreign pseudotime, signal, and membership cells and against
duplicate family key/cell and membership cell rows. All fail at core coverage
before dispatch. The same test accepts a complete-axis source with a 12-cell
session subset and checks the bound row order. The existing core contract in
`tests/test_typed_effective_cell_coverage_contract.R:77-108` checks subset
ordering and duplicate vector rejection. Direct-runner mocks bypass core
coverage, so only method-owned malformed values are asserted there: mixed
known/unknown feature-set and signal-key selectors fail before the first
writer. Repeated pseudotime key/cell rows belong to the core coverage check.
The staged-writer fault and seed/oracle checks remain direct-runner tests;
they do not claim a full ABI transaction guarantee.

### Model-derived inputs and direct-runner evidence

The six Slingshot model runners call one `model_payload` reader with the
accepted G1 shape. Core's self-linked G1 adapter verifies the private model
copy and supplies `runtime_private_path(primary_path)` and a reconstructed
`fitted_source_space` matrix (`scworkbench/R/operation_pack_v3_models.R`,
accepted derived-axis branch, lines 665-713). The pack requires a valid
fitted-space matrix and compares its exact ordered rownames with those of
`slingshot::slingPseudotime()` before selecting effective session cells.
The resulting pseudotime matrix provides finite values and `L1`...`Ln`
lineage keys.
Slingshot curves and lineage weights influence the values returned by
`slingPseudotime()`, but none of these six old-core methods uses their raw
coordinates or weights after that call. The pack never reads a caller's
`model_ref$pseudotime_table`, `model_path`, or bundle field.

| Operation | Model-derived information used by the old-core method |
|---|---|
| `trajectory_de` | Per-lineage pseudotime and ordered fitted cells drive Spearman tests and slopes. The model engine is recorded in both terminal metadata objects. |
| `trajectory_smooth` | Per-lineage pseudotime and cells determine bins and matrix slices. No other model metadata enters its outputs. |
| `trajectory_signature_smooth` | The same pseudotime and cells determine bins and signature scores; the feature-set members come from the typed registry snapshot. |
| `trajectory_vector_smooth` | The same pseudotime and cells align the typed signal family before aggregation. |
| `trajectory_geneset_dynamics` | Pseudotime and lineage selection drive scores, curves, bootstrap intervals and permutation nulls. The source kind, selected lineage key, and post-selection cell count enter output metadata. |
| `trajectory_signal_dynamics` | Pseudotime and lineage selection drive signal curves and bootstrap intervals. The source kind, selected lineage key, and post-selection cell count enter output metadata. |

`tests/test_s2_model_execution_contract.R` serializes a minimal valid
`SlingshotDataSet` with two exact pseudotime curves and passes it through a
mock G1 reader with the same callable interface. This synthetic model is not
a fitted or reopened G1 artifact. Its six direct runners match every frozen
old-core terminal kind, table type, row count, and canonical identity. Tests
also compare alternate lineage and parameter choices to the old core, check
bounded matrix reads, deterministic seeds and RNG restoration, and reject
unknown keys, invalid model payloads, and disjoint fitted axes before any
terminal write. The mock G1 reader supplies a fitted-space matrix separately
from the model RDS. A reordered axis, a different subset of the same parent
universe, an absent or malformed fitted-space field all fail with zero
terminal writes. A full fitted axis with a smaller valid session subset
passes and matches the old-core method output. These tests establish pack
numerical equivalence and axis binding under a
synthetic Slingshot projection; they do not establish real G1 reopen or ABI
publication.

The pack derives Monocle3 root-to-leaf keys and closest-vertex membership
directly from the CDS loaded through the verified private G1 reader. It
compares the global pseudotime's exact ordered cells with G1's authenticated
fitted axis before deriving the relation or writing terminals. It requires a
named connected tree, valid roots, complete fitted-cell membership, and
valid vertices. Caller pseudotime, lineage mapping, model path, and model
metadata cannot replace the CDS. The old-core relation in
`scworkbench/R/op_trajectory_de.R:43-70` and `R/op_monocle3.R:38-174`
provides the independent branch oracle. The 24-cell synthetic branching CDS
tests all six model runners on both branch keys and a smaller session subset.
The 18 old-core cases are frozen in
`tests/fixtures/trajectory_s2_monocle3_reference.json` with an
`exact_canonical_value` comparison policy. The test compares every
terminal's table, metadata, canonical identity, kind, and output order with
both the frozen record and the old core, and checks the old-core effective
parameters against the record. Reordered or foreign fitted cells
and invalid graph or membership fail with no terminal writes.

## G1 integration dependency

The old oracle pin predates the accepted generic G1/G2 derived fitted-axis
capability. Its self-linked reader compared the selected fitted-cell identity
with the full parent-space axis. The accepted generic work now reconstructs
the authenticated fitted subset and exposes it as `fitted_source_space` to
the private model reader. This pack checks the model's exact ordered cells
against that matrix. The direct numerical fixtures still use a mock G1
reader. A separate lifecycle test uses the accepted generic G1 fixture pack
to seal a Slingshot RDS on a 12-of-24 fitted axis, promote and finalize its
bundle, relocate the bundle, and obtain a verified private G1 reader with
the reconstructed fitted axis. That read executes the unregistered
`trajectory_smooth` runner, and both terminal identities match the old core.
Altering either the retained model tree or
selected-axis evidence rejects during G1 preparation, before S2 computation
or terminal writes. This is a real retained-state reopen with a synthetic
fitted model; it is not registered S2 dispatch or proof for all six methods.
The pack does not synthesize excluded cells or bypass G1 checks.

## Integration checks still owed

At grouped S2 integration, execute all six consumers through registered
dispatch on retained S3 or S5 fitted models. Add pack promotion,
rollback, replay-plan, pack-free inspection, and missing-pack dispatch checks
at the grouped S2 gate. Preserve each terminal's `exact_canonical_value`
policy; the frozen method oracle alone does not prove operation-level or
biological equivalence.
