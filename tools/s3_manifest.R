s3_identity <- function(value, scheme, scope, field) {
  value[[field]] <- getFromNamespace(
    ".wb_operation_pack_v2_identity_scoped", "scworkbench")(
      scheme, scope, value)
  value
}
s3_model_output <- function(group, source, siblings,
    derived = FALSE) {
  result <- list(name = "model", category = "terminal_artifact",
    kind = "model_ref", storage_profile = "model_ref_contained_tree_v1",
    allocated_id_key = "model_id", required = TRUE,
    promotion_group = group,
    contained_model_sibling_closure_policy = list(
      policy = "staged_contained_model_sibling_closure_v1",
      engine = "slingshot", model_type = "slingshot",
      primary_path = "model.rds", source_space_input = source,
      dimensions_parameter = NULL,
      required_sibling_terminals = as.list(siblings),
      file_order = "utf8_path_lexicographic_v1",
      maximum_file_count = 8L, maximum_total_bytes = 1073741824))
  if (derived) {
    policy <- list(
      schema = "scworkbench.operation_pack.derived_fitted_axis_policy.v1",
      mode = "source_relative_subset_v1", source_space_input = source,
      model_output = "model",
      required_siblings = lapply(intersect(siblings,
        c("pseudotime", "lineage_weights")), function(name) list(
          name = name, kind = "cell_vector_family",
          coverage_policy = "per_key_source_relative_subset_v1")),
      declaration_discipline = "exactly_once_before_cell_indexed_writes",
      maximum_selected_cells = 1000000L,
      maximum_combined_id_bytes = 67108864L,
      maximum_id_bytes = 1024L, maximum_family_keys = 64L)
    result$derived_fitted_axis_policy <- s3_identity(policy,
      policy$schema, "policy excluding policy_identity", "policy_identity")
  }
  result
}
s3_sparse_policy <- function(name, source) {
  scheme <- "scworkbench.operation_pack.standalone_sparse_cell_family_policy.v1"
  policy <- list(schema = scheme, output = name,
    policy = "standalone_sparse_cell_family_v1",
    axis_source = source,
    key_order = "declared_verbatim_v1",
    row_order = "key_major_then_source_relative_cell_v1",
    per_key_coverage = "per_key_source_relative_subset_v1",
    missing_values = "forbidden",
    declaration_discipline = "exactly_once_before_write",
    maximum_axis_cells = 1000000L, maximum_keys = 64L,
    maximum_rows = 1000000L, maximum_id_bytes = 1024L,
    maximum_control_bytes = 67108864L)
  s3_identity(policy, scheme, "policy excluding policy_identity",
    "policy_identity")
}
s3_read_model <- list(name = "fitted_model", input = "model",
  reader = "read_self_linked_contained_model_payload",
  allowed_models = list(list(engine = "slingshot", model_type = "slingshot")),
  linkage_policy = "authenticated_existing_source_and_cell_linkage_v1",
  maximum_file_count = 8L, maximum_total_bytes = 1073741824,
  call_discipline = "exactly_once")
s3_family_output <- function(name, group) {
  value <- output(name, "cell_vector_family", group)
  value$value_kind <- "continuous"
  value
}
s3_params <- function(include_fit = FALSE) {
  properties <- list(
    start.clus = string(), end.clus = strings(),
    exclude_clusters = strings(), stretch = number(2, 0, 10),
    reweight = boolean(TRUE))
  if (include_fit) properties <- c(properties, list(
    include_clusters = strings(), exclude_cells = strings(),
    shrink = number(1, 0, 1), maxit = integer(10L, 1L, 200L, FALSE)))
  list(type = "object", properties = properties,
    additionalProperties = FALSE)
}
s3_operation <- function(name, title, description, inputs, params, outputs,
    group, policies = list(), model_read = FALSE, adoption = NULL) {
  names <- vapply(outputs, `[[`, character(1), "name")
  result <- list(name = name, qualified_name = paste0("scwbtrajectory::", name),
    title = title, description = description, runner = paste0(".run_", name),
    inputs = inputs, params_schema = params,
    dimension_policies = policies, outputs = outputs,
    dependencies = list(dependency("scworkbench", "wb_prim_model_ref"),
      dependency("slingshot", "slingshot", ">=2.8.0,<3.0.0"),
      dependency("SingleCellExperiment", "SingleCellExperiment",
        ">=1.0.0,<3.0.0")),
    promotion_closure = closure(names),
    replay_equivalence = stats::setNames(lapply(names, function(x)
      if (identical(x, "model")) "loadable_reference" else
        "exact_canonical_value"), names))
  # The optional public-interface metadata restricts names to underscores.
  # Leave it absent so the catalog projects the established start.clus and
  # end.clus parameter names directly from the admitted parameter schema.
  if (model_read) result$self_linked_contained_model_reads <-
    list(s3_read_model)
  if (!is.null(adoption)) result$supported_external_state_adoptions <-
    list(adoption)
  result
}

s3_legacy <- s3_operation("slingshot", "Slingshot trajectory",
  "Infer lineage pseudotime and display curves in an authenticated 2D space.",
  list(input("space_2d", "space_2d"), input("clusters", "cell_vector")),
  s3_params(), list(s3_family_output("pseudotime", "trajectory"),
    output("curves", "path_2d_family", "trajectory")), "trajectory",
  list(s3_sparse_policy("pseudotime", list(
    source = "typed_space_input_v1", input = "space_2d")),
    family_policy("curves")))

s3_fit <- s3_operation("slingshot_fit", "Slingshot fitted model",
  "Fit Slingshot to an authenticated space and publish a contained model with its lineage families.",
  list(input("fit_space", "space"), input("clusters", "cell_vector")),
  s3_params(TRUE), list(s3_model_output("fit", "fit_space",
    c("pseudotime", "lineage_weights"), TRUE),
    s3_family_output("pseudotime", "fit"),
    s3_family_output("lineage_weights", "fit")), "fit")

s3_weights <- s3_operation("slingshot_lineage_weights",
  "Slingshot lineage weights",
  "Extract sparse lineage weights from an authenticated fitted model without refitting.",
  list(input("model", "model_ref")),
  list(type = "object", properties = list(), additionalProperties = FALSE),
  list(s3_family_output("lineage_weights", "weights")), "weights",
  list(s3_sparse_policy("lineage_weights", list(
    source = "verified_g1_model_fitted_axis_v1",
    model_read = "fitted_model"))), model_read = TRUE)

s3_project <- s3_operation("slingshot_project",
  "Project Slingshot curves",
  "Project verified fitted curves into an authenticated two-dimensional display space.",
  list(input("model", "model_ref"), input("plot_space", "space_2d")),
  list(type = "object", properties = list(stretch = nullable_number(0, 10),
    shrink = nullable_number(0, 1),
    approx_points = integer(200L, 0L, 5000L, FALSE)),
    additionalProperties = FALSE),
  list(output("curves", "path_2d_family", "projection")),
  "projection", list(family_policy("curves")), model_read = TRUE)

s3_operations <- list(s3_legacy, s3_fit, s3_weights, s3_project)
