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
input <- function(name, kind, required = TRUE) list(name = name, kind = kind,
  required = required, alignment = switch(kind,
    matrix = "cell_and_feature_ordered", space = "cell_ordered",
    space_2d = "cell_ordered", cell_vector = "cell_ordered",
    cell_vector_family = "cell_ordered", relation_family = "cell_ordered",
    space_path_family = "associated_space"))
output <- function(name, kind, group, table_type = NULL) {
  value <- list(name = name, category = "terminal_artifact", kind = kind,
    storage_profile = if (kind == "path_2d_family") "path_2d_family_v1" else
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

manifest <- list(schema_version = "scworkbench.operation_pack.v3",
  pack_abi_version = "3",
  identity_scheme = "scworkbench.operation_pack.canonical_json.v1",
  namespace = "scwbtrajectory",
  package = list(name = "scwbTrajectory", version = "0.1.0"), backend = "R",
  compatibility = list(scworkbench = list(minimum = "0.1.0",
    maximum_exclusive = "0.2.0"), pack_abi = list(minimum = "3", maximum = "3")),
  entry_point = "scwbTrajectory::register_operation_pack",
  operations = list(branch, route, score, project))

path <- file.path(root, "inst", "scworkbench", "operation-pack.json")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE,
  null = "null", digits = NA), path, useBytes = TRUE)
