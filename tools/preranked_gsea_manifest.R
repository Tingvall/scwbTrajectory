# Canonical signal-interval identity and trajectory-resolved preranked GSEA.
interval_output_names <- c("signal_curve", "global_test", "local_tests",
  "signal_intervals", "adjustment_summary", "fit_summary")
interval_params <- s3_regions_params
interval_params$global_q_cutoff <- interval_params$dynamic_q_cutoff
interval_params$local_q_cutoff <- interval_params$region_q_cutoff
interval_params$dynamic_q_cutoff <- NULL
interval_params$region_q_cutoff <- NULL
interval_provenance <- s3_regions_provenance
interval_operation <- list(
  name = "trajectory_signal_intervals",
  qualified_name = "scwbtrajectory::trajectory_signal_intervals",
  title = "Trajectory Signal Intervals",
  description = paste(
    "Test whole-trajectory dynamics of authenticated cell-level signals and",
    "report permutation-supported elevated or reduced pseudotime intervals."),
  runner = ".run_trajectory_signal_intervals",
  inputs = list(input("signal", "cell_vector_family"),
    input("pseudotime", "cell_vector_family"),
    input("covariates", "cell_vector_family", FALSE),
    input("depth", "cell_vector", FALSE),
    input("sample_id", "cell_vector", FALSE),
    input("membership", "cell_vector", FALSE)),
  params_schema = list(type = "object", properties = interval_params,
    required = list("pseudotime_key"), additionalProperties = FALSE),
  provenance_reads = interval_provenance,
  dimension_policies = lapply(interval_output_names, table_policy),
  outputs = list(
    output("signal_curve", "obs_table", "trajectory_signal_intervals",
      "trajectory_signal_intervals_curve"),
    output("global_test", "obs_table", "trajectory_signal_intervals",
      "trajectory_signal_intervals_global_test"),
    output("local_tests", "obs_table", "trajectory_signal_intervals",
      "trajectory_signal_intervals_local_tests"),
    output("signal_intervals", "obs_table", "trajectory_signal_intervals",
      "trajectory_signal_intervals"),
    output("adjustment_summary", "obs_table",
      "trajectory_signal_intervals",
      "trajectory_signal_intervals_adjustment"),
    output("fit_summary", "obs_table", "trajectory_signal_intervals",
      "trajectory_signal_intervals_summary")),
  dependencies = list(dependency("scworkbench", "wb_prim_obs_table")),
  promotion_closure = closure(interval_output_names),
  replay_equivalence = equivalence(interval_output_names),
  retained_parameters = c(lapply(names(interval_params), function(name)
    retained_parameter(name, "effective_parameter_value",
      required = identical(name, "pseudotime_key") ||
        !is.null(interval_params[[name]]$default))), list(
    retained_row("n_pseudotime_finite", "fit_summary",
      "n_pseudotime_finite"),
    retained_row("n_signal_keys_analyzed", "fit_summary",
      "n_signal_keys_analyzed"),
    retained_row("n_signal_intervals", "fit_summary",
      "n_signal_intervals"))),
  public_interface = interface(c("trajectory", "signal", "intervals",
    "permutation"), c("signal", "pseudotime", "covariates", "depth",
      "sample_id", "membership"), list(
      c("pseudotime_key", "string", "Required explicit lineage key."),
      c("signal_keys", "string_vec", "Optional signal keys."),
      c("covariate_keys", "string_vec", "Optional covariate keys."),
      c("pseudotime_direction", "choice", "As-is or inverted orientation."),
      c("window_width_fraction", "dbl", "Signal window width fraction."),
      c("window_step_fraction", "dbl", "Signal window step fraction."),
      c("min_cells_per_window", "int", "Minimum cells per window."),
      c("min_samples_per_window", "int", "Minimum samples per window."),
      c("summary_fun", "choice", "Window mean or median."),
      c("n_permutations", "int", "Permutation count."),
      c("seed", "int", "Deterministic base seed."),
      c("covariate_adjustment", "choice", "None, linear, or spline."),
      c("covariate_spline_df", "int", "Spline degrees of freedom."),
      c("adjust_depth", "bool", "Remove only the fitted depth component."),
      c("depth_field", "string", "Depth identity; defaults to nCount_RNA."),
      c("adjustment_spline_df", "int",
        "Pseudotime spline flexibility; defaults to 4."),
      c("permutation_mode", "choice", "Automatic, blocked, or unblocked."),
      c("include_endpoint_extrema", "bool", "Allow endpoint extrema."),
      c("global_q_cutoff", "dbl", "Global signal-test q cutoff."),
      c("local_q_cutoff", "dbl", "Local signal-test q cutoff."),
      c("merge_gap_fraction", "dbl", "Maximum interval merge gap.")),
    interval_output_names))

gsea_output_names <- c("global_enrichment", "window_enrichment",
  "enrichment_intervals", "leading_edge", "ranked_features",
  "fit_summary")
gsea_params <- list(
  namespace = string(),
  pseudotime_key = string(),
  analysis_scope = string("both", c("global", "local", "both"), FALSE),
  pseudotime_direction = string("as_is", c("as_is", "invert"), FALSE),
  global_rank_method = string("monotonic_association",
    c("monotonic_association", "early_late_contrast"), FALSE),
  early_fraction = number(0.2, 1e-06, 0.499999),
  late_fraction = number(0.2, 1e-06, 0.499999),
  window_width_fraction = number(0.15, 1e-06, 0.999999),
  window_step_fraction = number(0.05, 1e-06, 1),
  min_cells_per_window = integer(20L, 2L, 10000000L, FALSE),
  min_samples_per_window = integer(1L, 1L, 1000000L, FALSE),
  permutation_mode = string("auto", c("auto", "blocked", "unblocked"),
    FALSE),
  n_permutations = integer(1000L, 2L, 1000000L, FALSE),
  seed = integer(1L, 1L, 2147483646L, FALSE),
  covariate_adjustment = string("none",
    c("none", "linear", "natural_spline"), FALSE),
  covariate_spline_df = integer(4L, 2L, 20L, FALSE),
  adjust_depth = boolean(FALSE),
  depth_field = list(type = c("string", "null"),
    default = "nCount_RNA", minLength = 1L),
  adjustment_spline_df = integer(4L, 2L, 10L, FALSE),
  min_set_size = integer(5L, 1L, 10000000L, FALSE),
  max_set_size = integer(5000L, 1L, 10000000L, FALSE),
  score_type = string("std", c("std", "pos", "neg"), FALSE),
  gsea_param = number(1, 1e-12),
  sample_size = integer(101L, 3L, 1000000L, FALSE),
  n_perm_simple = integer(1000L, 1L, 10000000L, FALSE),
  eps = number(1e-10, 0, 1),
  p_adjust_method = string("BH", stats::p.adjust.methods, FALSE),
  global_q_cutoff = number(0.05, 1e-12, 1),
  local_q_cutoff = number(0.05, 1e-12, 1),
  merge_gap_fraction = number(0.05, 0, 1),
  consensus_leading_edge_fraction = number(0.5, 1e-12, 1),
  chunk_size = integer(128L, 1L, 1000000L, FALSE))
gsea_provenance <- lapply(c("pseudotime", "membership", "sample_id",
  "covariates", "depth"), function(name) list(
    name = paste0(name, "_ancestry"),
    category = "authenticated_core_provenance_snapshot",
    kind = "artifact_ancestry", reader = "read_artifact_ancestry",
    input = name, activation = "when_resolved_from_artifact",
    lifecycle_policy = "require_current"))
gsea_public_type <- function(schema) {
  type <- as.character(unlist(schema$type, use.names = FALSE))
  if ("integer" %in% type) return("int")
  if ("number" %in% type) return("dbl")
  if ("boolean" %in% type) return("bool")
  if ("string" %in% type) return(if (is.null(schema$enum)) "string" else
    "choice")
  "list"
}
gsea_operation <- list(
  name = "trajectory_preranked_gsea",
  qualified_name = "scwbtrajectory::trajectory_preranked_gsea",
  title = "Trajectory Preranked GSEA",
  description = paste(
    "Run full-trajectory and fixed-window preranked GSEA using an",
    "authenticated matrix feature universe and feature-set registry snapshot."),
  runner = ".run_trajectory_preranked_gsea",
  inputs = list(input("matrix", "matrix"),
    input("pseudotime", "cell_vector_family"),
    input("membership", "cell_vector", FALSE),
    input("sample_id", "cell_vector", FALSE),
    input("covariates", "cell_vector_family", FALSE),
    input("depth", "cell_vector", FALSE)),
  params_schema = list(type = "object", properties = gsea_params,
    required = list("namespace", "pseudotime_key"),
    additionalProperties = FALSE),
  registry_reads = list(feature_registry),
  provenance_reads = gsea_provenance,
  typed_slice_reads = list(matrix_read("chunk_size")),
  dimension_policies = lapply(gsea_output_names, table_policy),
  outputs = list(
    output("global_enrichment", "obs_table", "trajectory_preranked_gsea",
      "trajectory_preranked_gsea_global_enrichment"),
    output("window_enrichment", "obs_table", "trajectory_preranked_gsea",
      "trajectory_preranked_gsea_window_enrichment"),
    output("enrichment_intervals", "obs_table", "trajectory_preranked_gsea",
      "trajectory_preranked_gsea_enrichment_intervals"),
    output("leading_edge", "obs_table", "trajectory_preranked_gsea",
      "trajectory_preranked_gsea_leading_edge"),
    output("ranked_features", "obs_table", "trajectory_preranked_gsea",
      "trajectory_preranked_gsea_ranked_features"),
    output("fit_summary", "obs_table", "trajectory_preranked_gsea",
      "trajectory_preranked_gsea_summary")),
  dependencies = list(dependency("scworkbench", "wb_prim_obs_table"),
    dependency("fgsea", "fgseaMultilevel", "==1.34.2"),
    dependency("BiocParallel", "SerialParam", ">=1.42.0")),
  promotion_closure = closure(gsea_output_names),
  replay_equivalence = equivalence(gsea_output_names),
  retained_parameters = c(lapply(names(gsea_params), function(name)
    retained_parameter(name, "effective_parameter_value",
      required = name %in% c("namespace", "pseudotime_key") ||
        !is.null(gsea_params[[name]]$default))), list(
    retained_row("matrix_id", "fit_summary", "matrix_id"),
    retained_row("fgsea_version", "fit_summary", "fgsea_version"),
    retained_row("feature_universe_size", "fit_summary",
      "feature_universe_size"))),
  public_interface = interface(c("trajectory", "enrichment", "gsea",
    "permutation"), c("matrix", "pseudotime", "membership", "sample_id",
    "covariates", "depth"), lapply(names(gsea_params), function(name)
      c(name, gsea_public_type(gsea_params[[name]]),
        paste("Trajectory GSEA parameter", name))),
    gsea_output_names))

trajectory_analysis_operations <- list(interval_operation, gsea_operation)
