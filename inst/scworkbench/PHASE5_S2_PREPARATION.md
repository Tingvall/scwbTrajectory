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
fixture can be reproduced. Direct pack runners use the same comparator.
The pack manifest now declares eight grouped S2 entries: six canonical model
names and two explicit-pseudotime dynamics names. The latter execute through
registered dispatch and pass terminal publication checks. Model entries
execute through registered staged and promoted G1 reads on the proposed core
delta below. All model terminals have frozen-oracle publication coverage.

| Operation | Prepared now | Remaining for S2 integration |
|---|---|
| `trajectory_de` | Registered staged and promoted Slingshot and Monocle3 reads match both terminals and frozen branch oracles; terminal group publishes. The runner processes one authenticated feature partition at a time. | Independent core-delta acceptance and grouped integration. |
| `trajectory_smooth` | Registered staged and promoted Slingshot and Monocle3 reads match both terminals; terminal group publishes. | Independent core-delta acceptance and grouped integration. |
| `trajectory_signature_smooth` | Registered staged and promoted Slingshot and Monocle3 reads match both terminals; terminal group publishes. | Independent core-delta acceptance and grouped integration. |
| `trajectory_geneset_dynamics` | Registered model reads match all six terminals and publish; explicit mode retains frozen publication and rollback checks. | Independent core-delta acceptance and grouped integration. |
| `trajectory_vector_smooth` | Registered model reads match both terminals and publish from an authenticated signal artifact. | Independent core-delta acceptance and grouped integration. |
| `trajectory_signal_dynamics` | Registered model reads match all four terminals and publish from an authenticated signal artifact; explicit mode retains its checks. | Independent core-delta acceptance and grouped integration. |

The DE runner scans the authenticated feature axis once in declared chunks.
It densifies at most the selected rows in the current chunk, discards that
block after accumulating per-feature rows, and adjusts p-values per lineage
after the scan. The frozen old-core fixture
`trajectory_s2_de_streaming_reference.json` covers all features and lineages
and a reordered explicit feature selection; chunk sizes 2 and 3 produce the
same terminal identities. The focused block observer records a largest dense
block no greater than `chunk_size`. Peak memory still includes the current
feature block, lineage pseudotime vectors, and the final per-feature output
rows. The accepted typed-reader discipline requests every authenticated
feature partition once, so a small feature selection still scans the full
feature axis, though it densifies only selected rows.

The explicit-pseudotime inputs are complete-axis `cell_vector_family`
primitives. The direct runners consume only the accepted typed readers and
return declared primitives to terminal writers; they do not receive a bundle,
database, registry handle, or committed path. The registered explicit modes
also pass promotion fault rollback and retry through the core ABI. This
preparation does not remove any core operation or assert migration completion.

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
reader. A separate lifecycle test uses a pack-owned copy of the accepted
generic G1 fixture pack at `tests/fixtures/s2-g1-derived-axis`, copied from
scworkbench `5433162180017f6cd1b58badad305cfc0188a7b0` fixture tree
`c30d6274a54bd964545d62140b7a98477b8894c6`. The test pins SHA-256 for
all four fixture files and fails on missing or changed inputs. It uses G1
to seal a Slingshot RDS on a 12-of-24 fitted axis, promote and finalize its
bundle, relocate the bundle, and obtain a verified private G1 reader with
the reconstructed fitted axis. That read executes the
`trajectory_smooth` runner directly, and both terminal identities match the old core.
Altering either the retained model tree or
selected-axis evidence rejects during G1 preparation, before S2 computation
or terminal writes. This is a real retained-state reopen with a synthetic
fitted model; it is not registered S2 dispatch or proof for all six methods.
The pack does not synthesize excluded cells or bypass G1 checks.

## Integration checks still owed

The proposed generic core delta is scworkbench
`ac7656c74af42160b6a1342ddfcd86f5c1b336a1`, directly descended from
accepted G1 core `434d55a280b0e03e90530bf272add63dfdb3064b`. PyBench
remains `9ad8a7a3bd2298fdc6f426032e54ccb07eee0402`. The delta excludes
transient deserialized `model_obj` from the contained `model_ref` request
identity. It has **not** received independent acceptance under this order.

With that exact pair, the pack test now executes all six registered model
operations through both `wb_ref_run()` and promoted `wb_ref_artifact()` G1
reads. It compares all declared terminal identities with the frozen old-core
Slingshot and both Monocle3 branch oracles. A replayable test signal comes
from the verified G1 sibling artifact. Every promoted terminal group passes
`exact_canonical_value` content verification and replay-input planning.
Slingshot publication faults roll back and retry for all six operations.
Model and terminal artifacts load after read-only relocation with the pack
registration removed. Missing-pack dispatch fails. Invalid lineage, model
type, contained primary bytes, and selected-axis evidence reject with no
durable output, provenance, transition-row, or payload-file residue; both
altered files are restored byte-for-byte and to their original permissions.

Remaining work is focused independent acceptance of the small core delta,
then grouped manifest integration with Forge-1's separate S3 branch and the
independent S2 migration audit. This pack branch does not remove old-core
owners or claim the combined catalog. Phase 5 remains **9/24**.

## Focused and package qualification

The following package check belongs to the predecessor pack commit
`d17df85bd53d09db2a21b2ee39717e87d5e338db`, before the DE streaming
correction. The isolated R library `/tmp/forge2-s2-ac7656c-lib` was installed
from scworkbench `ac7656c74af42160b6a1342ddfcd86f5c1b336a1` and that
pack commit. The exact PyBench source was
`/home/dirk/Projects/.pybench-worktrees/phase5-g1-model-reference-dispatch-quill-20260923`
at `9ad8a7a3bd2298fdc6f426032e54ccb07eee0402`.

From the pack worktree:

```sh
R_LIBS=/tmp/forge2-s2-ac7656c-lib Rscript -e 'test_dir <- "tests"; source("tests/helper_trajectory_fixture.R",local=.GlobalEnv); source("tests/test_s2_registered_publication_contract.R",local=.GlobalEnv); n <- sort(ls(.GlobalEnv,pattern="^test_trajectory_s2_registered_")); for(x in n){get(x,.GlobalEnv)();cat("PASS",x,"\n")};cat("RESULT",length(n),"passed\n")'
```

Result: **8/8 passed, exit 0**. The broader pack contract runner found
**34/34 passed** inside the final package check. From the core delta
worktree, with that R library, the exact PyBench source on `PYTHONPATH`, and
`SCWORKBENCH_PYBENCH_PYTHON` set to the project virtual environment:

```sh
Rscript inst/tests/run_contract_tests.R --test-file=tests/test_operation_pack_v3_g1_model_reference_dispatch_repair_contract.R
```

Result: **3/3 passed, exit 0**. From the pinned PyBench worktree:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /home/dirk/Projects/pybench/.venv/bin/python -m pytest -q -p no:cacheprovider tests/test_phase5_model_capabilities_v3.py
```

Result: **31/31 passed, exit 0**. The pack was built in
`/tmp/forge2-s2-ac7656c-final-check` with `R CMD build
--no-build-vignettes <pack-worktree>` and checked there with the exact R
library and PyBench source using `R CMD check --no-manual
scwbTrajectory_0.1.1.tar.gz`: both **exit 0**, check **Status: OK**.
Changed R files parse and `git diff --check` passes. The earlier exploratory
check, before test-file source portability and import declarations were fixed,
ended with 1 ERROR, 1 WARNING, and 1 NOTE; it is superseded by the final
clean check.

Sable and Meridian-9 should independently review the two-file core delta
`434d55a..ac7656c` and this pinned S2 pack candidate as separate decisions.
The core delta is not accepted by the pack qualification alone.

For the DE streaming correction, the temporary test library was rebuilt from
a clean archive of core `ac7656c` and the current pack source. The focused
`test_trajectory_s2_de_streaming_block_bound_contract` passed (**1/1**);
the registered S2 suite passed **9/9**. The new frozen DE oracle was generated
with old core `ba24635611756240f21f6817044883b832f12701`, and old-core
chunk sizes 2 and 3 produced identical identities for both cases. Changed
R and JSON files parsed and `git diff --check` passed. A new broader package
check awaits the final accepted core/pack pair, as directed by the streaming
work order. Sable should review the bounded-block correction before that
paired qualification is used for migration acceptance.
