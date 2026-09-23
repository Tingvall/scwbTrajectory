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
| `trajectory_de` | Unregistered Slingshot model runner writes `de_table` and `fit_summary` from G1 lineage pseudotime, a bounded typed matrix, and a feature-set snapshot. Both frozen identities match. | Register the G1 and matrix grants; real fitted-axis reopen; Monocle3 branch projection. |
| `trajectory_smooth` | Unregistered Slingshot model runner writes `trend_table` and `fit_summary` using the shared bin and aggregate equations. Both frozen identities match. | Same G1, matrix, registration, and Monocle3 work. |
| `trajectory_signature_smooth` | Unregistered Slingshot model runner uses an authenticated feature-set snapshot and bounded matrix reads. Both frozen identities match. | Same G1, matrix, registration, and Monocle3 work. |
| `trajectory_geneset_dynamics` | Unregistered explicit and Slingshot model runners share the six-terminal numerical method. Both frozen model and explicit terminal closures match. | Register each source mode and its grants; real G1 reopen, publication audit, and Monocle3 branch projection. |
| `trajectory_vector_smooth` | Unregistered Slingshot model runner uses a typed signal family. Both frozen identities match. | Register signal and G1 grants; real reopen and Monocle3 branch projection. |
| `trajectory_signal_dynamics` | Unregistered explicit and Slingshot model runners share the four-terminal numerical method. Both frozen model and explicit terminal closures match. | Register each source mode and its grants; real G1 reopen, publication audit, and Monocle3 branch projection. |

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

The Monocle3 route remains at the G1 projection interface. The old core
gets global `monocle3::pseudotime(cds)` and separately derives a keyed
cell-to-lineage relation from the ordered principal graph and closest-vertex
mapping (`scworkbench/R/op_trajectory_de.R:43-70`,
`R/op_monocle3.R:38-174`). Reproducer: an authenticated `monocle3_cds`
with two branches has one global pseudotime vector but requires that keyed
relation before any of the six consumers can select lineage cells. The
Slingshot matrix projection alone cannot represent it. A pack-owned
projection from a verified CDS, or a typed relation bound to the same fit,
must be reviewed with the completed G1 contract. An unsupported model type
fails before terminal writes; no synthetic pseudotime field is accepted as
an alternate production source.

## G1 integration dependency

The old oracle pin predates the accepted generic G1/G2 derived fitted-axis
capability. Its self-linked reader compared the selected fitted-cell identity
with the full parent-space axis. The accepted generic work now reconstructs
the authenticated fitted subset and exposes it as `fitted_source_space` to
the private model reader. This pack checks the model's exact ordered cells
against that matrix, but its direct tests use a mock G1 reader. A real S2
fixture must still reopen a sealed 12-of-24 fitted model, execute the six
consumers, and reject reordered or tampered fitted-axis evidence before
dispatch. The pack does not synthesize excluded cells or bypass G1 checks.

## Integration checks still owed

Once the generic derived-axis capability and fitted-model producer are
accepted, replace the synthetic model oracle with a retained S3 or S5 model,
reopen it, and compare the six pack terminal closures. Add pack promotion,
rollback, replay-plan, pack-free inspection, and missing-pack dispatch checks
at the grouped S2 gate. Preserve each terminal's `exact_canonical_value`
policy; the frozen method oracle alone does not prove operation-level or
biological equivalence.
