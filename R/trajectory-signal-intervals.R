# Canonical terminology for permutation-qualified signal intervals.  The
# numerical implementation remains shared with the 0.1.4
# trajectory_signal_regions compatibility operation.
.trajectory_intervals_rename_columns <- function(table, mapping) {
  for (old in names(mapping)) {
    hit <- match(old, names(table))
    if (!is.na(hit)) names(table)[hit] <- mapping[[old]]
  }
  table
}

.trajectory_signal_intervals_outputs <- function(result) {
  curve <- result$outputs$signal_table$table
  global <- .trajectory_intervals_rename_columns(
    result$outputs$dynamic_tests$table,
    c(p_dynamic = "p_global", q_dynamic = "q_global",
      dynamic_significant = "global_significant",
      observed_dynamic_range = "observed_signal_range",
      dynamic_exceedances = "global_exceedances"))
  local <- .trajectory_intervals_rename_columns(
    result$outputs$candidate_regions$table,
    c(candidate_id = "local_test_id", merged_region_id = "interval_id",
      region_exceedances = "local_exceedances",
      region_supported = "local_supported",
      p_dynamic = "p_global", q_dynamic = "q_global",
      dynamic_significant = "global_significant"))
  intervals <- .trajectory_intervals_rename_columns(
    result$outputs$regions$table,
    c(region_id = "interval_id",
      constituent_candidate_ids = "constituent_local_test_ids",
      p_dynamic = "p_global", q_dynamic = "q_global"))
  adjustment <- result$outputs$adjustment_summary$table
  summary <- .trajectory_intervals_rename_columns(
    result$outputs$fit_summary$table,
    c(n_candidate_regions = "n_local_tests",
      n_supported_regions = "n_signal_intervals",
      dynamic_q_cutoff = "global_q_cutoff",
      region_q_cutoff = "local_q_cutoff",
      region_null = "local_null",
      dynamic_null = "global_null"))
  if ("region_supported" %in% names(curve))
    names(curve)[names(curve) == "region_supported"] <- "local_supported"
  if ("region_id" %in% names(curve))
    names(curve)[names(curve) == "region_id"] <- "interval_id"

  list(
    signal_curve = scworkbench::wb_prim_obs_table(curve,
      table_type = "trajectory_signal_intervals_curve",
      entity_kind = "trajectory_window", meta = list(
        operation = "trajectory_signal_intervals",
        lineage_id = as.character(summary$lineage_id[[1L]]),
        pseudotime_artifact_id = as.character(
          summary$pseudotime_artifact_id[[1L]]), overlapping_windows = TRUE)),
    global_test = scworkbench::wb_prim_obs_table(global,
      table_type = "trajectory_signal_intervals_global_test",
      entity_kind = "signal"),
    local_tests = scworkbench::wb_prim_obs_table(local,
      table_type = "trajectory_signal_intervals_local_tests",
      entity_kind = "signal_local_test"),
    signal_intervals = scworkbench::wb_prim_obs_table(intervals,
      table_type = "trajectory_signal_intervals",
      entity_kind = "signal_interval"),
    adjustment_summary = scworkbench::wb_prim_obs_table(adjustment,
      table_type = "trajectory_signal_intervals_adjustment",
      entity_kind = "signal"),
    fit_summary = scworkbench::wb_prim_obs_table(summary,
      table_type = "trajectory_signal_intervals_summary",
      entity_kind = "lineage"))
}

.run_trajectory_signal_intervals <- function(readers, writers, params,
                                             context) {
  legacy_params <- params
  if (!is.null(legacy_params$global_q_cutoff)) {
    legacy_params$dynamic_q_cutoff <- legacy_params$global_q_cutoff
    legacy_params$global_q_cutoff <- NULL
  }
  if (!is.null(legacy_params$local_q_cutoff)) {
    legacy_params$region_q_cutoff <- legacy_params$local_q_cutoff
    legacy_params$local_q_cutoff <- NULL
  }
  result <- .trajectory_s2_call(.trajectory_signal_regions,
    list(effective_cells = as.character(context$effective_cells),
      signal = readers$signal(), pseudotime = readers$pseudotime(),
      covariates = .trajectory_regions_optional_reader(readers,
        "covariates"),
      depth = .trajectory_regions_optional_reader(readers, "depth"),
      sample_id = .trajectory_regions_optional_reader(readers,
        "sample_id"),
      membership = .trajectory_regions_optional_reader(readers,
        "membership"),
      signal_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "signal_ancestry"),
      pseudotime_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "pseudotime_ancestry"),
      covariates_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "covariates_ancestry"),
      depth_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "depth_ancestry"),
      sample_id_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "sample_id_ancestry"),
      membership_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "membership_ancestry")), legacy_params)
  result$outputs <- .trajectory_signal_intervals_outputs(result)
  result$params$global_q_cutoff <- params$global_q_cutoff %||% 0.05
  result$params$local_q_cutoff <- params$local_q_cutoff %||% 0.05
  .trajectory_s2_publish(result, writers, c("signal_curve", "global_test",
    "local_tests", "signal_intervals", "adjustment_summary",
    "fit_summary"))
}
