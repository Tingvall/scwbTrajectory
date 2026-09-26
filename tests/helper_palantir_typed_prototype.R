# A transient contract fixture invokes the unregistered scwbTrajectory runners
# through the real core ABI. It is never a production pack registration.
palantir_proto_array <- function(x) unname(as.list(x))

palantir_proto_output <- function(name, kind, group = "fit",
    profile = "terminal_primitive_v1", value_kind = NULL,
    table_type = NULL) {
  value <- list(name = name, category = "terminal_artifact", kind = kind,
    storage_profile = profile, allocated_id_key = paste0(name, "_id"),
    required = TRUE, promotion_group = group)
  if (!is.null(value_kind)) value$value_kind <- value_kind
  if (!is.null(table_type)) value$table_type <- table_type
  value
}

palantir_proto_manifest <- function(namespace = "scwbpalantirproto",
    package_name = "scwbPalantirProto") {
  closure <- sort(c("pseudotime", "branch_membership", "branch_summary",
    "entropy", "model"))
  relation_policy <- list(
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
    maximum_control_bytes = 67108864L, maximum_chunk_pairs = 8192L)
  relation_policy$policy_identity <- getFromNamespace(
    ".wb_operation_pack_v2_identity_scoped", "scworkbench")(
      relation_policy$schema, "policy excluding policy_identity",
      relation_policy)
  one_cell <- function(output) list(output = output,
    policy = "authenticated_output_axes_v1", shape = "one_per_cell",
    cell_axis = list(source = "effective_cells_v1", input = "fit_space"),
    row_order = "input_cell_order_v1", missing_values = "forbidden")
  params_fit <- list(type = "object", properties = list(
    root_cell = list(type = "string", minLength = 1L),
    terminal_cells = list(type = "array", items = list(type = "string",
      minLength = 1L), minItems = 1L, uniqueItems = TRUE),
    dims = list(type = "array", items = list(type = "integer",
      minimum = 1L), minItems = 2L, uniqueItems = TRUE),
    knn = list(type = "integer", minimum = 2L),
    n_diffusion_components = list(type = "integer", minimum = 2L),
    n_waypoints = list(type = "integer", minimum = 1L),
    use_early_cell_as_start = list(type = "boolean"),
    seed = list(type = "integer", minimum = 0L)),
    required = palantir_proto_array(c("root_cell", "terminal_cells",
      "dims", "knn", "n_diffusion_components", "n_waypoints",
      "use_early_cell_as_start", "seed")), additionalProperties = FALSE)
  model <- palantir_proto_output("model", "model_ref",
    profile = "model_ref_contained_tree_v1")
  model$contained_model_sibling_closure_policy <- list(
    policy = "staged_contained_model_sibling_closure_v1",
    engine = "palantir", model_type = "palantir_fit",
    primary_path = "model.pkl", source_space_input = "fit_space",
    dimensions_parameter = "dims",
    # G2's sibling-axis check is exact one-per-cell; the variable-row
    # relation is instead sealed independently and group-closed atomically.
    required_sibling_terminals = palantir_proto_array(c("pseudotime",
      "entropy")),
    file_order = "utf8_path_lexicographic_v1", maximum_file_count = 8L,
    maximum_total_bytes = 1073741824L)
  model$derived_fitted_axis_policy <- list(
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
    maximum_id_bytes = 1024L, maximum_family_keys = 64L)
  model$derived_fitted_axis_policy$policy_identity <- getFromNamespace(
    ".wb_operation_pack_v2_identity_scoped", "scworkbench")(
      model$derived_fitted_axis_policy$schema,
      "policy excluding policy_identity",
      model$derived_fitted_axis_policy)
  fit <- list(name = "fit_typed", qualified_name =
      paste0(namespace, "::fit_typed"),
    title = "Private Palantir fit fixture",
    description = "Test-only typed Palantir fit computation.",
    runner = ".fit_typed", inputs = list(list(name = "fit_space",
      kind = "space", required = TRUE, alignment = "cell_ordered")),
    params_schema = params_fit,
    dimension_policies = list(one_cell("pseudotime"), relation_policy,
      list(output = "branch_summary", policy = "data_dependent_table_v1"),
      one_cell("entropy")),
    outputs = list(model,
      palantir_proto_output("pseudotime", "cell_vector",
        value_kind = "continuous"),
      palantir_proto_output("branch_membership", "relation_family"),
      palantir_proto_output("branch_summary", "obs_table",
        table_type = "palantir_branch_summary"),
      palantir_proto_output("entropy", "cell_vector",
        value_kind = "continuous")), dependencies = list(),
    promotion_closure = stats::setNames(lapply(closure, function(x)
      palantir_proto_array(closure)), closure),
    replay_equivalence = list(model = "loadable_reference",
      pseudotime = "exact_canonical_value",
      branch_membership = "exact_canonical_value",
      branch_summary = "exact_canonical_value",
      entropy = "exact_canonical_value"))
  params_project <- list(type = "object", properties = list(
    prob_min = list(type = "number", minimum = 0, maximum = 1),
    smoothness = list(type = "number", minimum = 0.05, maximum = 20),
    min_cells_per_branch = list(type = "integer", minimum = 2L),
    approx_points = list(type = "integer", minimum = 0L,
      maximum = 5000L),
    pickle_execution_trust = list(type = "string",
      enum = list("operator_verified_source_v1"))),
    required = palantir_proto_array(c("prob_min", "smoothness",
      "min_cells_per_branch", "approx_points",
      "pickle_execution_trust")),
    additionalProperties = FALSE)
  project <- list(name = "project_typed", qualified_name =
      paste0(namespace, "::project_typed"),
    title = "Private Palantir project fixture",
    description = "Test-only typed Palantir project computation.",
    runner = ".project_typed", inputs = list(
      list(name = "model", kind = "model_ref", required = TRUE,
        alignment = "model_linked"),
      list(name = "plot_space", kind = "space_2d", required = TRUE,
        alignment = "cell_ordered")), params_schema = params_project,
    dimension_policies = list(list(output = "curves",
      policy = "data_dependent_family_v1")),
    outputs = list(palantir_proto_output("curves", "path_2d_family",
      group = "projection", profile = "path_2d_family_v1")),
    dependencies = list(),
    promotion_closure = list(curves = list("curves")),
    replay_equivalence = list(curves = "exact_canonical_value"),
    self_linked_contained_model_reads = list(list(name = "fitted_model",
      input = "model", reader =
        "read_self_linked_contained_model_payload",
      allowed_models = list(list(engine = "palantir",
        model_type = "palantir_fit")),
      linkage_policy =
        "authenticated_existing_source_and_cell_linkage_v1",
      maximum_file_count = 8L, maximum_total_bytes = 1073741824L,
      call_discipline = "exactly_once")))
  list(schema_version = "scworkbench.operation_pack.v3",
    pack_abi_version = "3",
    identity_scheme = "scworkbench.operation_pack.canonical_json.v1",
    namespace = namespace,
    package = list(name = package_name, version = "0.0.1"),
    backend = "R", compatibility = list(
      scworkbench = list(minimum = "0.1.0",
        maximum_exclusive = "0.2.0"),
      pack_abi = list(minimum = "3", maximum = "3")),
    entry_point = paste0(package_name, "::register_operation_pack"),
    operations = list(fit, project))
}

.palantir_proto_activation_state <- new.env(parent = emptyenv())
.palantir_proto_activation_state$next_id <- 0L

palantir_proto_activate <- function() {
  .palantir_proto_activation_state$next_id <-
    .palantir_proto_activation_state$next_id + 1L
  suffix <- as.character(.palantir_proto_activation_state$next_id)
  package_name <- paste0("scwbPalantirProto", suffix)
  namespace <- paste0("scwbpalantirproto", suffix)
  root <- tempfile("scwb-palantir-fixture-")
  source <- file.path(root, "source")
  library <- file.path(root, "library")
  dir.create(file.path(source, "R"), recursive = TRUE)
  dir.create(file.path(source, "inst", "scworkbench"), recursive = TRUE)
  dir.create(library)
  writeLines(c(paste0("Package: ", package_name), "Type: Package",
    "Title: Transient Palantir Typed Prototype Fixture",
    "Version: 0.0.1", "Authors@R: person('Contract', 'Fixture', email = 'noreply@example.invalid', role = c('aut', 'cre'))",
    "Description: Temporary test-only Palantir typed ABI fixture.",
    "License: MIT", "Encoding: UTF-8",
    "Imports: scworkbench, scwbTrajectory"),
    file.path(source, "DESCRIPTION"))
  writeLines("export(register_operation_pack)",
    file.path(source, "NAMESPACE"))
  writeLines(c(
    ".fit_typed <- function(readers, writers, params, context)",
    "  getFromNamespace('.trajectory_palantir_typed_fit', 'scwbTrajectory')(readers, writers, params, context)",
    ".project_typed <- function(readers, writers, params, context)",
    "  getFromNamespace('.trajectory_palantir_typed_project', 'scwbTrajectory')(readers, writers, params, context)",
    "register_operation_pack <- function() scworkbench::wb_operation_pack_register(",
    paste0("  system.file('scworkbench', 'operation-pack.json', package = '",
      package_name, "'),"),
    paste0("  list('", namespace, "::fit_typed' = .fit_typed,"),
    paste0("       '", namespace,
      "::project_typed' = .project_typed))")),
    file.path(source, "R", "register.R"))
  manifest <- file.path(source, "inst", "scworkbench",
    "operation-pack.json")
  jsonlite::write_json(palantir_proto_manifest(namespace, package_name),
    manifest,
    auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  conformance <- scworkbench::wb_operation_pack_manifest_validate(manifest)
  stopifnot(isTRUE(conformance$ok))
  log <- file.path(root, "install.log")
  status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL",
    paste0("--library=", shQuote(library)), shQuote(source)),
    stdout = log, stderr = log)
  if (!identical(status, 0L)) stop(paste(readLines(log,
    warn = FALSE), collapse = "\n"), call. = FALSE)
  previous <- .libPaths()
  .libPaths(c(library, previous))
  scworkbench::wb_operation_pack_activate(
    system.file("scworkbench", "operation-pack.json",
      package = package_name, lib.loc = library),
    lib.loc = library)
  list(root = root, previous = previous,
    fit = paste0(namespace, "::fit_typed"),
    project = paste0(namespace, "::project_typed"),
    cleanup = function() {
      try(scworkbench::wb_operation_pack_unload(namespace),
        silent = TRUE)
      if (package_name %in% loadedNamespaces())
        try(unloadNamespace(package_name), silent = TRUE)
      .libPaths(previous)
      unlink(root, recursive = TRUE, force = TRUE)
    })
}
