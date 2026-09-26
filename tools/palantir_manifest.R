# Candidate S5 definitions. Core retains the historical operations until a
# separate ownership order; these names are qualified pack operations only.
palantir_policy_identity <- function(value) {
  value$policy_identity <- getFromNamespace(
    ".wb_operation_pack_v2_identity_scoped", "scworkbench")(
      value$schema, "policy excluding policy_identity", value)
  value
}
palantir_one_cell <- function(name) list(output = name,
  policy = "authenticated_output_axes_v1", shape = "one_per_cell",
  cell_axis = list(source = "effective_cells_v1", input = "fit_space"),
  row_order = "input_cell_order_v1", missing_values = "forbidden")
palantir_relation_policy <- palantir_policy_identity(list(
  schema = "scworkbench.operation_pack.sealed_relation_family_policy.v1",
  output = "branch_membership", policy = "sealed_relation_family_v1",
  axis_source = list(source = "g2_derived_fitted_axis_v1",
    model_output = "model"),
  key_order = "first_occurrence_v1", row_order = "declared_exact_v1",
  pair_multiplicity = "unique_key_cell_v1",
  numeric_fields = "declare_presence_finite_v1",
  declaration_discipline = "exactly_once_before_write",
  maximum_axis_cells = 1000000L, maximum_keys = 1024L,
  maximum_rows = 1000000L, maximum_id_bytes = 1024L,
  maximum_control_bytes = 67108864L, maximum_chunk_pairs = 8192L))
palantir_model <- list(name = "model", category = "terminal_artifact",
  kind = "model_ref", storage_profile = "model_ref_contained_tree_v1",
  allocated_id_key = "model_id", required = TRUE,
  promotion_group = "palantir_fit",
  contained_model_sibling_closure_policy = list(
    policy = "staged_contained_model_sibling_closure_v1",
    engine = "palantir", model_type = "palantir_fit",
    primary_path = "model.pkl", source_space_input = "fit_space",
    dimensions_parameter = "dims",
    required_sibling_terminals = list("pseudotime", "entropy"),
    file_order = "utf8_path_lexicographic_v1", maximum_file_count = 8L,
    maximum_total_bytes = 1073741824L),
  derived_fitted_axis_policy = palantir_policy_identity(list(
    schema = "scworkbench.operation_pack.derived_fitted_axis_policy.v1",
    mode = "source_relative_subset_v1", source_space_input = "fit_space",
    model_output = "model", required_siblings = list(
      list(name = "pseudotime", kind = "cell_vector",
        coverage_policy = "exact_fitted_axis_v1"),
      list(name = "entropy", kind = "cell_vector",
        coverage_policy = "exact_fitted_axis_v1")),
    declaration_discipline = "exactly_once_before_cell_indexed_writes",
    maximum_selected_cells = 1000000L,
    maximum_combined_id_bytes = 67108864L,
    maximum_id_bytes = 1024L, maximum_family_keys = 64L)))
palantir_vector <- function(name) {
  value <- output(name, "cell_vector", "palantir_fit")
  value$value_kind <- "continuous"
  value
}
palantir_fit_outputs <- c("pseudotime", "branch_membership",
  "branch_summary", "entropy", "model")
palantir_fit <- list(name = "palantir_fit",
  qualified_name = "scwbtrajectory::palantir_fit",
  title = "Palantir fitted trajectory",
  description = paste(
    "Fit Palantir from an authenticated typed space in private Python",
    "computation; publish one contained model, pseudotime, ordered branch",
    "probabilities, branch summary, and entropy atomically."),
  runner = ".run_palantir_fit", inputs = list(input("fit_space", "space")),
  params_schema = list(type = "object", properties = list(
    root_cell = list(type = "string", minLength = 1L),
    terminal_cells = list(type = "array", minItems = 1L,
      uniqueItems = TRUE, items = list(type = "string", minLength = 1L)),
    dims = list(type = "array", minItems = 2L, uniqueItems = TRUE,
      items = list(type = "integer", minimum = 1L)),
    knn = integer(30L, 2L, 1000L, FALSE),
    n_diffusion_components = integer(10L, 2L, 200L, FALSE),
    n_waypoints = integer(300L, 1L, 10000L, FALSE),
    use_early_cell_as_start = boolean(TRUE),
    seed = integer(20L, 0L, 999999L, FALSE)),
    required = list("root_cell", "terminal_cells", "dims"),
    additionalProperties = FALSE),
  dimension_policies = list(palantir_one_cell("pseudotime"),
    palantir_relation_policy, table_policy("branch_summary"),
    palantir_one_cell("entropy")),
  outputs = list(palantir_model, palantir_vector("pseudotime"),
    output("branch_membership", "relation_family", "palantir_fit"),
    output("branch_summary", "obs_table", "palantir_fit",
      "palantir_branch_summary"), palantir_vector("entropy")),
  dependencies = list(dependency("scworkbench", "wb_prim_model_ref")),
  promotion_closure = closure(palantir_fit_outputs),
  replay_equivalence = list(model = "loadable_reference",
    pseudotime = "exact_canonical_value",
    branch_membership = "exact_canonical_value",
    branch_summary = "exact_canonical_value",
    entropy = "exact_canonical_value"))
palantir_project <- list(name = "palantir_project",
  qualified_name = "scwbtrajectory::palantir_project",
  title = "Project Palantir fitted branches",
  description = paste(
    "Project a core-authenticated Palantir model into an authenticated 2D",
    "display space. Its contained pickle is executed only after an explicit",
    "operator-verified-source attestation; checksum validity alone never",
    "establishes safety for an untrusted bundle."),
  runner = ".run_palantir_project",
  inputs = list(input("model", "model_ref"),
    input("plot_space", "space_2d")),
  params_schema = list(type = "object", properties = list(
    prob_min = number(0.05, 0, 1),
    smoothness = number(1, 0.05, 20),
    min_cells_per_branch = integer(50L, 2L, 1000000L, FALSE),
    approx_points = integer(200L, 0L, 5000L, FALSE),
    pickle_execution_trust = list(type = "string",
      enum = list("operator_verified_source_v1"))),
    required = list("pickle_execution_trust"),
    additionalProperties = FALSE),
  dimension_policies = list(family_policy("curves")),
  outputs = list(output("curves", "path_2d_family", "projection")),
  dependencies = list(dependency("scworkbench", "wb_prim_model_ref")),
  promotion_closure = closure("curves"),
  replay_equivalence = list(curves = "exact_canonical_value"),
  self_linked_contained_model_reads = list(list(name = "fitted_model",
    input = "model", reader = "read_self_linked_contained_model_payload",
    allowed_models = list(list(engine = "palantir",
      model_type = "palantir_fit")),
    linkage_policy = "authenticated_existing_source_and_cell_linkage_v1",
    maximum_file_count = 8L, maximum_total_bytes = 1073741824L,
    call_discipline = "exactly_once")))
palantir_operations <- list(palantir_fit, palantir_project)
