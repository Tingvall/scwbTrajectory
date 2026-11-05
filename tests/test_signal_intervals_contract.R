.trajectory_intervals_fixture <- function(n = 120L) {
  cells <- sprintf("interval_cell_%03d", seq_len(n))
  pseudotime <- seq(0, 1, length.out = n)
  signals <- cbind(
    peak = 8 * exp(-((pseudotime - 0.3) / 0.05)^2),
    trough = -8 * exp(-((pseudotime - 0.7) / 0.05)^2))
  list(cells = cells,
    signal = scworkbench::wb_prim_cell_vector_family(
      rep(colnames(signals), each = n), rep(cells, ncol(signals)),
      as.numeric(signals), value_type = "continuous"),
    pseudotime = scworkbench::wb_prim_cell_vector_family(
      rep("lineage", n), cells, pseudotime, value_type = "continuous"))
}

test_trajectory_signal_intervals_manifest_contract <- function() {
  trajectory_activate_pack()
  manifest <- scworkbench::wb_operation_pack_inspect(
    "scwbtrajectory")$manifest
  operations <- stats::setNames(manifest$operations,
    vapply(manifest$operations, `[[`, character(1), "name"))
  canonical <- operations$trajectory_signal_intervals
  legacy <- operations$trajectory_signal_regions
  stopifnot(!is.null(canonical), !is.null(legacy),
    identical(canonical$qualified_name,
      "scwbtrajectory::trajectory_signal_intervals"),
    identical(canonical$title, "Trajectory Signal Intervals"),
    identical(vapply(canonical$outputs, `[[`, character(1), "name"),
      c("signal_curve", "global_test", "local_tests", "signal_intervals",
        "adjustment_summary", "fit_summary")),
    identical(vapply(legacy$outputs, `[[`, character(1), "name"),
      c("signal_table", "dynamic_tests", "candidate_regions", "regions",
        "adjustment_summary", "fit_summary")),
    canonical$params_schema$properties$global_q_cutoff$default == 0.05,
    canonical$params_schema$properties$local_q_cutoff$default == 0.05,
    identical(canonical$params_schema$properties$adjust_depth$default, TRUE),
    identical(canonical$params_schema$properties$depth_field$default,
      "nCount_RNA"),
    all(vapply(canonical$replay_equivalence, identical, logical(1),
      "exact_canonical_value")))
  TRUE
}

test_trajectory_signal_intervals_numerical_alias_contract <- function() {
  x <- .trajectory_intervals_fixture()
  method <- getFromNamespace(".trajectory_signal_regions",
    "scwbTrajectory")
  translate <- getFromNamespace(".trajectory_signal_intervals_outputs",
    "scwbTrajectory")
  set.seed(919L)
  before <- .Random.seed
  legacy <- method(effective_cells = x$cells, signal = x$signal,
    pseudotime = x$pseudotime, pseudotime_key = "lineage",
    window_width_fraction = 0.12, window_step_fraction = 0.02,
    min_cells_per_window = 6L, n_permutations = 99L, seed = 31L,
    adjust_depth = FALSE,
    dynamic_q_cutoff = 0.1, region_q_cutoff = 0.2)
  canonical <- translate(legacy)
  stopifnot(identical(before, .Random.seed),
    identical(names(canonical), c("signal_curve", "global_test",
      "local_tests", "signal_intervals", "adjustment_summary",
      "fit_summary")),
    identical(canonical$signal_curve$table$value,
      legacy$outputs$signal_table$table$value),
    identical(canonical$global_test$table$p_global,
      legacy$outputs$dynamic_tests$table$p_dynamic),
    identical(canonical$local_tests$table$p_empirical,
      legacy$outputs$candidate_regions$table$p_empirical),
    identical(canonical$signal_intervals$table$effect_size,
      legacy$outputs$regions$table$effect_size),
    all(c("local_test_id", "local_supported", "interval_id",
      "p_global", "q_global") %in% names(canonical$local_tests$table)),
    all(c("interval_id", "constituent_local_test_ids") %in%
      names(canonical$signal_intervals$table)))
  TRUE
}

test_trajectory_signal_intervals_registered_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  time <- seq(0, 1, length.out = length(x$cells))
  signal <- scworkbench::wb_prim_cell_vector_family(
    rep("peak", length(x$cells)), x$cells,
    8 * exp(-((time - 0.4) / 0.12)^2), value_type = "continuous")
  pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("L1", length(x$cells)), x$cells, time,
    value_type = "continuous")
  before <- list(artifacts = scworkbench::wb_art_list(x$bundle),
    fingerprint = getFromNamespace(".wb_session_source_identity",
      "scworkbench")(x$bundle)$scientific_state_fingerprint)
  result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbtrajectory::trajectory_signal_intervals", inputs = list(
      signal = signal, pseudotime = pseudotime, pseudotime_key = "L1",
      adjust_depth = FALSE,
      window_width_fraction = 0.3, window_step_fraction = 0.05,
      min_cells_per_window = 3L, n_permutations = 39L, seed = 37L,
      global_q_cutoff = 0.2, local_q_cutoff = 0.2, min_cells = 1L))
  after <- list(artifacts = scworkbench::wb_art_list(x$bundle),
    fingerprint = getFromNamespace(".wb_session_source_identity",
      "scworkbench")(x$bundle)$scientific_state_fingerprint)
  stopifnot(identical(names(result$run$outputs), c("signal_curve",
      "global_test", "local_tests", "signal_intervals",
      "adjustment_summary", "fit_summary")), identical(before, after))
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    result$session, outputs = names(result$run$outputs),
    created_by = "signal-intervals-contract")
  stopifnot(all(promoted$promoted$status == "ok"),
    setequal(promoted$promoted$output, names(result$run$outputs)))
  TRUE
}
