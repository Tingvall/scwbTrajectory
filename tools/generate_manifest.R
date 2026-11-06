args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)

string <- function(default = NULL, choices = NULL, nullable = is.null(default)) {
  out <- list(type = if (nullable) c("string", "null") else "string")
  if (!is.null(default)) out$default <- default
  if (!is.null(choices)) out$enum <- as.list(choices)
  out
}
number <- function(default, minimum = NULL, maximum = NULL) {
  out <- list(type = "number", default = default)
  if (!is.null(minimum)) out$minimum <- minimum
  if (!is.null(maximum)) out$maximum <- maximum
  out
}
integer <- function(default, minimum = NULL, maximum = NULL,
                    nullable = is.null(default)) {
  out <- list(type = if (nullable) c("integer", "null") else "integer")
  if (!is.null(default)) out$default <- default
  if (!is.null(minimum)) out$minimum <- minimum
  if (!is.null(maximum)) out$maximum <- maximum
  out
}
strings <- function() list(type = c("array", "null"), uniqueItems = TRUE,
  items = list(type = "string", minLength = 1L))
required_strings <- function() list(type = "array", minItems = 1L,
  uniqueItems = TRUE, items = list(type = "string", minLength = 1L))
boolean <- function(default = FALSE) list(type = "boolean", default = default)
endpoint <- function() list(type = c("object", "string", "null"))
artifact_ref <- function() list(type = c("object", "string", "null"))
input <- function(name, kind, required = TRUE) list(name = name, kind = kind,
  required = required, alignment = switch(kind,
    matrix = "cell_and_feature_ordered", space = "cell_ordered",
    space_2d = "cell_ordered", cell_vector = "cell_ordered",
    cell_vector_family = "cell_ordered", relation_family = "cell_ordered",
    space_path_family = "associated_space",
    path_2d_family = "canonical_row_ordered",
    model_ref = "model_linked", obs_table = "none"))
output <- function(name, kind, group, table_type = NULL) {
  value <- list(name = name, category = "terminal_artifact", kind = kind,
    storage_profile = if (kind == "path_2d_family") "path_2d_family_v1" else
      if (kind == "space_path_family") "space_path_family_v1" else
        "terminal_primitive_v1", allocated_id_key = paste0(name, "_id"),
    required = TRUE, promotion_group = group)
  if (!is.null(table_type)) value$table_type <- table_type
  if (kind == "cell_vector_family") value$value_kind <- "continuous"
  value
}
dependency <- function(name, probe, spec = ">=0.1.0,<0.2.0") list(
  requirement_id = paste0(tolower(name), "_runtime"), kind = "r_package",
  name = name, version = list(scheme = "r_package_version_constraint_v1",
    specifier = spec), disposition = "required", import_probe = probe)
retained_parameter <- function(name, source = "parameter_value",
                               parameter = name, required = FALSE) {
  out <- list(name = name, source = source, required = required)
  if (!is.null(parameter)) out$parameter <- parameter
  out
}
retained_row <- function(name, output, field) list(name = name,
  source = "output_first_row", output = output, field = field)
retained_meta <- function(name, output, field) list(name = name,
  source = "output_meta", output = output, field = field)
table_policy <- function(name) list(output = name,
  policy = "data_dependent_table_v1")
matrix_read <- function(chunk_parameter) list(name = "matrix_chunks",
  input = "matrix", reader = "read_typed_slices", kind = "matrix",
  supported_axes = list("features"), chunk_axis = "features",
  chunk_size_parameter = chunk_parameter,
  call_discipline = "ordered_partition_exactly_once",
  empty_selection_policy = "reject", duplicate_policy = "reject",
  selection = list(cells = list(policy = "effective_input_axis_v1"),
    features = list(policy = "effective_input_axis_v1")))
feature_registry <- list(name = "feature_sets",
  category = "authenticated_core_registry_snapshot",
  kind = "feature_set_registry", reader = "read_feature_set_registry",
  matrix_input = "matrix", no_selector_policy = "all_namespaces",
  empty_match_policy = "allow")
interface <- function(intent, inputs, params, outputs) list(
  intent = as.list(intent),
  inputs = lapply(inputs, function(name) list(name = name, label = name,
    help = paste("Typed", name, "input."))),
  parameters = lapply(params, function(item) list(name = item[[1L]],
    type = item[[2L]], help = item[[3L]])),
  outputs = lapply(outputs, function(name) list(name = name,
    description = paste("Declared", name, "terminal."))))

branch_params <- list(branch_a = string(), branch_b = string(),
  feature_set = string(), features = strings(), min_pct = number(0, 0, 1),
  p_adjust_method = string("BH", nullable = FALSE),
  chunk_size = integer(128L, 1L, 100000L, FALSE))
branch_outputs <- c("branch_de_table", "branch_summary")
branch <- list(name = "trajectory_branch_compare",
  qualified_name = "scwbtrajectory::trajectory_branch_compare",
  title = "Trajectory Branch Compare",
  description = "Compare two trajectory branches feature-by-feature using authenticated matrix, pseudotime, and branch-membership inputs.",
  runner = ".run_trajectory_branch_compare",
  inputs = list(input("matrix", "matrix"),
    input("pseudotime", "cell_vector_family"),
    input("branch_membership", "relation_family")),
  params_schema = list(type = "object", properties = branch_params,
    additionalProperties = FALSE),
  registry_reads = list(feature_registry),
  dimension_policies = lapply(branch_outputs, table_policy),
  outputs = list(output("branch_de_table", "feature_table", "branch",
    "trajectory_branch_compare"),
    output("branch_summary", "obs_table", "branch",
      "trajectory_branch_compare_summary")),
  dependencies = list(dependency("scworkbench", "wb_prim_feature_table")),
  promotion_closure = stats::setNames(rep(list(as.list(sort(branch_outputs,
      method = "radix"))), 2L),
    branch_outputs), replay_equivalence = stats::setNames(
      as.list(rep("exact_canonical_value", 2L)), branch_outputs),
  typed_slice_reads = list(matrix_read("chunk_size")),
  retained_parameters = list(retained_meta("matrix_id", "branch_de_table",
      "matrix_id"), retained_parameter("branch_a", required = TRUE),
    retained_parameter("branch_b", required = TRUE),
    retained_row("n_cells_a", "branch_summary", "n_cells_a"),
    retained_row("n_cells_b", "branch_summary", "n_cells_b"),
    retained_row("n_overlap_dropped", "branch_summary", "n_overlap_dropped"),
    retained_row("n_features", "branch_summary", "n_features_tested"),
    retained_parameter("min_pct", "effective_parameter_value", required = TRUE),
    retained_parameter("p_adjust_method", "effective_parameter_value",
      required = TRUE)),
  public_interface = interface(c("trajectory", "comparison", "statistics"),
    c("matrix", "pseudotime", "branch_membership"), list(
      c("branch_a", "string", "Branch/lineage key A."),
      c("branch_b", "string", "Branch/lineage key B."),
      c("feature_set", "string", "Optional matrix feature-set name."),
      c("features", "string_vec", "Optional explicit feature IDs."),
      c("min_pct", "dbl", "Minimum detection fraction."),
      c("p_adjust_method", "string", "P-value adjustment method."),
      c("chunk_size", "int", "Authenticated feature chunk size.")),
    branch_outputs))

route_params <- list(route_key = string(), feature_set = string(),
  features = strings(), n_bins = integer(50L, 5L, 1000L, FALSE),
  summary_fun = string("mean", c("mean", "median"), FALSE),
  min_cells_per_bin = integer(5L, 1L, 100000L, FALSE),
  weight_mode = string("probability", c("probability", "binary_threshold"),
    FALSE), prob_min = number(0, 0, 1),
  chunk_size = integer(128L, 1L, 100000L, FALSE))
route_outputs <- c("trend_table", "fit_summary")
route <- list(name = "trajectory_route_smooth",
  qualified_name = "scwbtrajectory::trajectory_route_smooth",
  title = "Trajectory Route Smoothing",
  description = "Summarize feature expression along an authenticated pseudotime route.",
  runner = ".run_trajectory_route_smooth",
  inputs = list(input("matrix", "matrix"), input("pseudotime", "cell_vector"),
    input("branch_membership", "relation_family")),
  params_schema = list(type = "object", properties = route_params,
    additionalProperties = FALSE), registry_reads = list(feature_registry),
  dimension_policies = lapply(route_outputs, table_policy),
  outputs = list(output("trend_table", "obs_table", "route",
    "trajectory_route_smooth"), output("fit_summary", "obs_table", "route",
      "trajectory_route_smooth_summary")),
  dependencies = list(dependency("scworkbench", "wb_prim_obs_table")),
  promotion_closure = stats::setNames(rep(list(as.list(sort(route_outputs,
      method = "radix"))), 2L),
    route_outputs), replay_equivalence = stats::setNames(
      as.list(rep("exact_canonical_value", 2L)), route_outputs),
  typed_slice_reads = list(matrix_read("chunk_size")),
  retained_parameters = list(retained_row("matrix_id", "fit_summary",
      "matrix_id"), retained_parameter("route_key", required = TRUE),
    retained_parameter("feature_set"),
    retained_row("n_features_input", "fit_summary", "n_features_input"),
    retained_row("n_cells_route", "fit_summary", "n_cells_route"),
    retained_row("n_bins", "fit_summary", "n_bins"),
    retained_parameter("summary_fun", "effective_parameter_value",
      required = TRUE), retained_parameter("min_cells_per_bin",
      "effective_parameter_value", required = TRUE),
    retained_parameter("weight_mode", "effective_parameter_value",
      required = TRUE), retained_parameter("prob_min",
      "effective_parameter_value", required = TRUE)),
  public_interface = interface(c("trajectory", "summary", "route"),
    c("matrix", "pseudotime", "branch_membership"), list(
      c("route_key", "string", "Route key in branch membership."),
      c("feature_set", "string", "Optional matrix feature-set name."),
      c("features", "string_vec", "Optional explicit feature IDs."),
      c("n_bins", "int", "Number of pseudotime bins."),
      c("summary_fun", "choice", "Mean or median summary."),
      c("min_cells_per_bin", "int", "Minimum cells per bin."),
      c("weight_mode", "choice", "Probability or binary weights."),
      c("prob_min", "dbl", "Minimum route probability."),
      c("chunk_size", "int", "Authenticated feature chunk size.")),
    route_outputs))

score_params <- list(score_mode = string("orthogonal_projection",
    "orthogonal_projection", FALSE), distance_metric = string("euclidean",
    "euclidean", FALSE), score_prefix = string())
score_outputs <- c("path_scores", "fit_summary")
score <- list(name = "space_path_score",
  qualified_name = "scwbtrajectory::space_path_score",
  title = "Score cells against analytical path",
  description = "Score authenticated space cells by projection onto an analytical path.",
  runner = ".run_space_path_score",
  inputs = list(input("space", "space"),
    input("path", "space_path_family")),
  params_schema = list(type = "object", properties = score_params,
    additionalProperties = FALSE),
  dimension_policies = list(list(output = "path_scores",
    policy = "authenticated_output_axes_v1", shape = "cell_key_cartesian",
    cell_axis = list(source = "effective_cells_v1", input = "space"),
    key_axis = list(source = "declared_keys_with_optional_prefix_v1",
      keys = as.list(c("progress", "distance", "nearest_segment",
        "nearest_anchor")), prefix_parameter = "score_prefix"),
    row_order = "key_major_then_cell_v1", missing_values = "forbidden"),
    table_policy("fit_summary")),
  outputs = list(output("path_scores", "cell_vector_family", "score"),
    output("fit_summary", "obs_table", "score",
      "space_path_score_summary")),
  dependencies = list(dependency("scworkbench", "wb_prim_cell_vector_family")),
  promotion_closure = stats::setNames(rep(list(as.list(sort(score_outputs,
      method = "radix"))), 2L),
    score_outputs), replay_equivalence = stats::setNames(
      as.list(rep("exact_canonical_value", 2L)), score_outputs),
  retained_parameters = list(retained_row("space_id", "fit_summary",
      "space_id"), retained_row("path_id", "fit_summary", "path_id"),
    retained_parameter("score_mode", "effective_parameter_value",
      required = TRUE), retained_parameter("distance_metric",
      "effective_parameter_value", required = TRUE),
    retained_row("n_cells", "fit_summary", "cells_scored")),
  public_interface = interface(c("trajectory", "geometry", "path"),
    c("space", "path"), list(
      c("score_mode", "choice", "Orthogonal-projection scoring mode."),
      c("distance_metric", "choice", "Euclidean distance metric."),
      c("score_prefix", "string", "Optional prefix for score keys.")),
    score_outputs))

project_params <- list(project_mode = string("local_interpolation",
    "local_interpolation", FALSE), approx_points = integer(NULL, 2L, 50000L))
project <- list(name = "space_path_project_2d",
  qualified_name = "scwbtrajectory::space_path_project_2d",
  title = "Project analytical path into 2D",
  description = "Project an analytical path from its authenticated source space into an authenticated 2D display space.",
  runner = ".run_space_path_project_2d",
  inputs = list(input("path", "space_path_family"),
    input("source_space", "space"), input("plot_space", "space_2d")),
  params_schema = list(type = "object", properties = project_params,
    additionalProperties = FALSE),
  dimension_policies = list(list(output = "curves",
    policy = "data_dependent_family_v1")),
  outputs = list(output("curves", "path_2d_family", "projection")),
  dependencies = list(dependency("scworkbench", "wb_prim_path_2d_family"),
    dependency("RANN", "nn2", ">=2.6.0,<3.0.0")),
  promotion_closure = list(curves = list("curves")),
  replay_equivalence = list(curves = "exact_canonical_value"),
  retained_parameters = list(retained_parameter("project_mode",
    "effective_parameter_value", required = TRUE),
    retained_parameter("approx_points", "effective_parameter_value")),
  public_interface = interface(c("trajectory", "geometry", "path", "display"),
    c("path", "source_space", "plot_space"), list(
      c("project_mode", "choice", "Local-interpolation projection mode."),
      c("approx_points", "int", "Optional display point count.")), "curves"))

source(file.path(root, "tools", "s4_manifest.R"), local = TRUE)
source(file.path(root, "tools", "s3_manifest.R"), local = TRUE)
source(file.path(root, "tools", "palantir_manifest.R"), local = TRUE)
source(file.path(root, "tools", "changepoint_manifest.R"), local = TRUE)
source(file.path(root, "tools", "preranked_gsea_manifest.R"), local = TRUE)

# Preserve the independently accepted S1/S2/S4 and S3 definitions. The
# accepted Monocle3 definitions are frozen from integrated commit 7e8c3fc;
# never regenerate or reinterpret them while appending Palantir.
manifest <- jsonlite::read_json(file.path(root, "tools",
  "integrated_s1_s2_s4_manifest.json"), simplifyVector = FALSE)
manifest$compatibility$scworkbench$minimum <- "0.1.1"
stopifnot(identical(length(manifest$operations), 17L),
  identical(vapply(manifest$operations, `[[`, character(1),
    "qualified_name"), c(vapply(list(branch, route, score, project),
      `[[`, character(1), "qualified_name"),
      vapply(s4_operations, `[[`, character(1), "qualified_name"),
      paste0("scwbtrajectory::", c("trajectory_de", "trajectory_smooth",
        "trajectory_signature_smooth", "trajectory_geneset_dynamics",
        "trajectory_vector_smooth", "trajectory_signal_dynamics",
        "trajectory_geneset_dynamics_explicit",
        "trajectory_signal_dynamics_explicit")))))

# The 0.1.6 depth-control update changes the two accepted S2 gene-set
# consumers in lockstep while retaining their distinct pseudotime sources.
upgrade_geneset_depth <- function(operation) {
  explicit <- grepl("_explicit$", operation$name)
  prefix_parameter <- if (explicit) "pseudotime_key" else "lineage"
  operation$params_schema$properties$feature_sets <- required_strings()
  if (!explicit) operation$params_schema$properties$lineage <- list(
    type = "string", minLength = 1L)
  operation$params_schema$required <- unique(c(
    operation$params_schema$required %||% list(), "feature_sets",
    if (explicit) list() else "lineage"))
  operation$params_schema$properties$adjust_depth <- boolean(TRUE)
  operation$params_schema$properties$depth_field <- list(
    type = c("string", "null"), default = "nCount_RNA", minLength = 1L)
  operation$params_schema$properties$adjustment_spline_df <-
    integer(4L, 2L, 10L, FALSE)
  operation$inputs <- c(operation$inputs,
    list(input("depth", "cell_vector", FALSE)))
  operation$provenance_reads <- c(operation$provenance_reads, list(list(
    name = "depth_ancestry",
    category = "authenticated_core_provenance_snapshot",
    kind = "artifact_ancestry", reader = "read_artifact_ancestry",
    input = "depth", activation = "when_resolved_from_artifact",
    lifecycle_policy = "require_current")))
  extra_outputs <- list(
    output("adjusted_scores", "cell_vector_family",
      "trajectory_geneset_dynamics"),
    output("adjustment_diagnostics", "obs_table",
      "trajectory_geneset_dynamics",
      "trajectory_geneset_depth_adjustment"))
  fit_index <- match("fit_summary", vapply(operation$outputs, `[[`,
    character(1), "name"))
  operation$outputs <- append(operation$outputs, extra_outputs,
    after = fit_index - 1L)
  operation$dimension_policies <- c(operation$dimension_policies,
    list(list(output = "adjusted_scores",
      policy = "authenticated_output_axes_v1",
      shape = "cell_key_cartesian",
      cell_axis = list(source = "effective_cells_v1", input = "matrix"),
      key_axis = list(source = "parameter_keys_with_optional_prefix_v1",
        keys_parameter = "feature_sets",
        prefix_parameter = prefix_parameter, separator = "::"),
      row_order = "key_major_then_cell_v1", missing_values = "forbidden"),
      table_policy("adjustment_diagnostics")))
  operation$retained_parameters <- lapply(
    operation$retained_parameters %||% list(), function(item) {
      if (as.character(item$name %||% "")[1L] %in%
          c("feature_sets", prefix_parameter)) item$required <- TRUE
      item
    })
  operation$dependencies <- lapply(operation$dependencies, function(item) {
    if (identical(as.character(item$kind)[1L], "r_package") &&
        identical(as.character(item$name)[1L], "scworkbench"))
      item$version$specifier <- ">=0.1.1,<0.2.0"
    item
  })
  names <- vapply(operation$outputs, `[[`, character(1), "name")
  operation$promotion_closure <- closure(names)
  operation$replay_equivalence <- equivalence(names)
  operation
}
for (qualified in c("scwbtrajectory::trajectory_geneset_dynamics",
    "scwbtrajectory::trajectory_geneset_dynamics_explicit")) {
  index <- match(qualified, vapply(manifest$operations, `[[`, character(1),
    "qualified_name"))
  stopifnot(!is.na(index))
  manifest$operations[[index]] <- upgrade_geneset_depth(
    manifest$operations[[index]])
}

geneset_score_params <- list(
  namespace = list(type = "string", minLength = 1L),
  pseudotime_key = list(type = "string", minLength = 1L),
  feature_sets = required_strings(),
  score_method = string("mean", c("mean", "sum", "rank_mean"), FALSE),
  adjust_depth = boolean(TRUE),
  depth_field = string("nCount_RNA", nullable = FALSE),
  adjustment_spline_df = integer(4L, 2L, 10L, FALSE),
  min_set_size = integer(5L, 1L, 10000000L, FALSE),
  max_set_size = integer(5000L, 1L, 10000000L, FALSE),
  pseudotime_direction = string("as_is", c("as_is", "invert"), FALSE),
  chunk_size = integer(128L, 1L, 100000L, FALSE))
geneset_score_outputs <- c("adjusted_scores", "adjustment_diagnostics",
  "fit_summary")
geneset_score_provenance <- lapply(c("pseudotime", "membership", "depth"),
  function(name) list(name = paste0(name, "_ancestry"),
    category = "authenticated_core_provenance_snapshot",
    kind = "artifact_ancestry", reader = "read_artifact_ancestry",
    input = name, activation = "when_resolved_from_artifact",
    lifecycle_policy = "require_current"))
geneset_scores <- list(
  name = "trajectory_geneset_scores",
  qualified_name = "scwbtrajectory::trajectory_geneset_scores",
  title = "Trajectory Gene-set Scores",
  description = paste(
    "Construct reusable cell-level gene-set scores after optional per-gene",
    "partial depth adjustment that preserves fitted pseudotime structure."),
  runner = ".run_trajectory_geneset_scores",
  inputs = list(input("matrix", "matrix"),
    input("pseudotime", "cell_vector_family"),
    input("membership", "cell_vector", FALSE),
    input("depth", "cell_vector", FALSE)),
  params_schema = list(type = "object", properties = geneset_score_params,
    required = list("namespace", "pseudotime_key", "feature_sets"),
    additionalProperties = FALSE),
  registry_reads = list(feature_registry),
  provenance_reads = geneset_score_provenance,
  typed_slice_reads = list(matrix_read("chunk_size")),
  dimension_policies = list(list(output = "adjusted_scores",
    policy = "authenticated_output_axes_v1", shape = "cell_key_cartesian",
    cell_axis = list(source = "effective_cells_v1", input = "matrix"),
    key_axis = list(source = "parameter_keys_with_optional_prefix_v1",
      keys_parameter = "feature_sets", prefix_parameter = "pseudotime_key",
      separator = "::"), row_order = "key_major_then_cell_v1",
    missing_values = "forbidden"),
    table_policy("adjustment_diagnostics"),
    table_policy("fit_summary")),
  outputs = list(
    output("adjusted_scores", "cell_vector_family",
      "trajectory_geneset_scores"),
    output("adjustment_diagnostics", "obs_table",
      "trajectory_geneset_scores", "trajectory_geneset_depth_adjustment"),
    output("fit_summary", "obs_table", "trajectory_geneset_scores",
      "trajectory_geneset_scores_summary")),
  dependencies = list(dependency("scworkbench",
    "wb_prim_cell_vector_family", ">=0.1.1,<0.2.0")),
  promotion_closure = closure(geneset_score_outputs),
  replay_equivalence = equivalence(geneset_score_outputs),
  retained_parameters = c(lapply(names(geneset_score_params), function(name)
    retained_parameter(name, "effective_parameter_value",
      required = name %in% c("namespace", "pseudotime_key", "feature_sets") ||
        !is.null(geneset_score_params[[name]]$default))), list(
    retained_meta("matrix_id", "adjusted_scores", "matrix_id"),
    retained_meta("pseudotime_artifact_id", "adjusted_scores",
      "pseudotime_artifact_id"),
    retained_meta("depth_artifact_id", "adjusted_scores",
      "depth_artifact_id"),
    retained_row("n_cells_used", "fit_summary", "n_cells_used"),
    retained_row("n_feature_sets_scored", "fit_summary",
      "n_feature_sets_scored"))),
  public_interface = interface(c("trajectory", "gene_set", "scoring"),
    c("matrix", "pseudotime", "membership", "depth"), list(
      c("namespace", "string", "Feature-set registry namespace."),
      c("pseudotime_key", "string", "Selected pseudotime lineage key."),
      c("feature_sets", "string_vec", "Ordered feature-set names."),
      c("score_method", "choice", "Mean, sum, or rank-mean scoring."),
      c("adjust_depth", "bool", "Apply per-gene partial depth adjustment."),
      c("depth_field", "string", "Depth identity; defaults to nCount_RNA."),
      c("adjustment_spline_df", "int",
        "Pseudotime spline flexibility; defaults to 4."),
      c("min_set_size", "int", "Minimum retained feature-set size."),
      c("max_set_size", "int", "Maximum retained feature-set size."),
      c("pseudotime_direction", "choice", "As-is or inverted pseudotime."),
      c("chunk_size", "int", "Authenticated feature chunk size.")),
    geneset_score_outputs))
monocle3_operations <- jsonlite::read_json(file.path(root, "tools",
  "accepted_s5_monocle3_operations.json"), simplifyVector = FALSE)
definition_identity <- getFromNamespace(".wb_operation_pack_v3_id",
  "scworkbench")
stopifnot(identical(vapply(monocle3_operations, `[[`, character(1),
  "qualified_name"), c("scwbtrajectory::monocle3_fit",
    "scwbtrajectory::monocle3_project")),
  identical(vapply(monocle3_operations, function(operation)
    definition_identity("scworkbench.operation_pack.definition.v3",
      operation)$sha256, character(1)), c(
        "5fb8c122ed1bee4edb1d0918b4cc89b288a61dd1b59218bdf021d76a92ca08e0",
        "c6fd76e839d7d473b933d9edc7ada746998c73b5ad9750d1d82ce34916d8576e")))
manifest$operations <- c(manifest$operations, s3_operations,
  monocle3_operations, palantir_operations, changepoint_operations,
  list(geneset_scores), trajectory_analysis_operations)
# Package and manifest versions move together. Version 0.1.8 authenticates
# dynamic multi-feature-set score-family axes.
manifest$package$version <- "0.1.8"

path <- file.path(root, "inst", "scworkbench", "operation-pack.json")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE,
  null = "null", digits = NA), path, useBytes = TRUE)
