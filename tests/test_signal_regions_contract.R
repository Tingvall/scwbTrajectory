.trajectory_regions_fixture <- function(n = 120L, missing_pt = 0L) {
  cells <- sprintf("cell_%03d", seq_len(n))
  pseudotime <- seq(0, 1, length.out = n)
  pt_value <- pseudotime
  if (missing_pt > 0L) pt_value[seq.int(n - missing_pt + 1L, n)] <- NA_real_
  signals <- cbind(
    known_peak = 8 * exp(-((pseudotime - 0.30) / 0.05)^2),
    known_trough = -8 * exp(-((pseudotime - 0.70) / 0.05)^2),
    null_constant = rep(0, n))
  list(cells = cells, pseudotime_value = pt_value,
    signal = scworkbench::wb_prim_cell_vector_family(
      rep(colnames(signals), each = n), rep(cells, ncol(signals)),
      as.numeric(signals), value_type = "continuous"),
    pseudotime = scworkbench::wb_prim_cell_vector_family(
      rep("lineage_1", n), cells, pt_value, value_type = "continuous"))
}

.trajectory_regions_call <- function(x, ...) {
  method <- getFromNamespace(".trajectory_signal_regions",
    "scwbTrajectory")
  method(effective_cells = x$cells, signal = x$signal,
    pseudotime = x$pseudotime, pseudotime_key = "lineage_1",
    window_width_fraction = 0.12, window_step_fraction = 0.02,
    min_cells_per_window = 6L, n_permutations = 99L, seed = 17L,
    adjust_depth = FALSE,
    dynamic_q_cutoff = 0.1, region_q_cutoff = 0.2, ...)
}

.trajectory_regions_promote <- function(bundle, primitive, output) {
  run <- scworkbench::wb_new_run("metadata_cell_vector",
    params = list(vector_name = paste("signal-regions", output)),
    outputs = list(vector = primitive),
    notes = "Owned generic trajectory-signal-regions fixture.")
  promoted <- scworkbench::wb_promote_analysis(bundle,
    stats::setNames(list(run), run$run_id), head_run_id = run$run_id,
    outputs = "vector", created_by = "signal-regions-contract")
  as.character(promoted$promoted$artifact_id[[1L]])
}

.trajectory_regions_state <- function(bundle) list(
  artifacts = scworkbench::wb_art_list(bundle),
  fingerprint = getFromNamespace(".wb_session_source_identity",
    "scworkbench")(bundle)$scientific_state_fingerprint)

test_trajectory_signal_regions_manifest_contract <- function() {
  trajectory_activate_pack()
  manifest <- scworkbench::wb_operation_pack_inspect(
    "scwbtrajectory")$manifest
  operations <- stats::setNames(manifest$operations,
    vapply(manifest$operations, `[[`, character(1), "name"))
  operation <- operations$trajectory_signal_regions
  stopifnot(!is.null(operation),
    identical(operation$qualified_name,
      "scwbtrajectory::trajectory_signal_regions"),
    identical(operation$params_schema$required, list("pseudotime_key")),
    identical(vapply(operation$inputs, `[[`, character(1), "name"),
      c("signal", "pseudotime", "covariates", "depth", "sample_id",
        "membership")),
    identical(vapply(operation$inputs, `[[`, logical(1), "required"),
      c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)),
    identical(vapply(operation$outputs, `[[`, character(1), "name"),
      c("signal_table", "dynamic_tests", "candidate_regions", "regions",
        "adjustment_summary", "fit_summary")),
    identical(operation$outputs[[1L]]$table_type,
      "trajectory_signal_dynamics_curve"),
    identical(operation$params_schema$properties$summary_fun$default,
      "median"),
    identical(as.integer(operation$params_schema$properties$
      n_permutations$default), 1000L),
    identical(operation$params_schema$properties$permutation_mode$default,
      "auto"),
    identical(operation$params_schema$properties$adjust_depth$default, TRUE),
    identical(operation$params_schema$properties$depth_field$default,
      "nCount_RNA"),
    length(operation$provenance_reads) == 6L,
    all(vapply(operation$replay_equivalence, identical, logical(1),
      "exact_canonical_value")))
  TRUE
}

test_trajectory_signal_regions_synthetic_determinism_contract <- function() {
  x <- .trajectory_regions_fixture()
  set.seed(991L)
  before_seed <- .Random.seed
  first <- .trajectory_regions_call(x)
  after_seed <- .Random.seed
  second <- .trajectory_regions_call(x)
  stopifnot(identical(before_seed, after_seed),
    identical(serialize(first$outputs, NULL, version = 3L),
      serialize(second$outputs, NULL, version = 3L)))

  dynamic <- first$outputs$dynamic_tests$table
  regions <- first$outputs$regions$table
  curve <- first$outputs$signal_table$table
  peak <- regions[regions$signal_key == "known_peak" &
    regions$direction == "elevated", , drop = FALSE]
  trough <- regions[regions$signal_key == "known_trough" &
    regions$direction == "reduced", , drop = FALSE]
  null <- dynamic[dynamic$signal_key == "null_constant", , drop = FALSE]
  stopifnot(nrow(peak) >= 1L, nrow(trough) >= 1L,
    min(abs(peak$peak_pseudotime - 0.30)) <= 0.04,
    min(abs(trough$peak_pseudotime - 0.70)) <= 0.04,
    all(peak$effect_size > 0), all(trough$effect_size < 0),
    identical(as.numeric(null$p_dynamic), 1),
    !isTRUE(null$dynamic_significant),
    !any(regions$signal_key == "null_constant"),
    all(c("feature_set", "pseudotime", "value", "value_raw",
      "n_cells", "local_direction", "p_empirical", "q_adjusted",
      "region_id") %in% names(curve)),
    identical(first$outputs$signal_table$table_type,
      "trajectory_signal_dynamics_curve"))
  TRUE
}

test_trajectory_signal_regions_covariate_block_and_missing_contract <-
    function() {
  n <- 120L
  x <- .trajectory_regions_fixture(n, missing_pt = 10L)
  time <- seq(0, 1, length.out = n)
  technical <- rep(seq(-2, 2, length.out = 20L), 6L)
  adjusted_peak <- 5 * technical +
    8 * exp(-((time - 0.42) / 0.05)^2)
  x$signal <- scworkbench::wb_prim_cell_vector_family(
    rep("adjusted_peak", n), x$cells, adjusted_peak,
    value_type = "continuous")
  covariates <- scworkbench::wb_prim_cell_vector_family(
    rep("technical_depth", n), x$cells, technical,
    value_type = "continuous")
  # Every block spans pseudotime so the blocked null has genuine support.
  sample_values <- rep(paste0("sample_", seq_len(6L)), length.out = n)
  samples <- scworkbench::wb_prim_cell_vector(
    stats::setNames(sample_values, x$cells), value_type = "categorical")
  result <- .trajectory_regions_call(x, covariates = covariates,
    sample_id = samples, covariate_adjustment = "linear",
    permutation_mode = "auto", min_samples_per_window = 2L)
  adjustment <- result$outputs$adjustment_summary$table
  summary <- result$outputs$fit_summary$table
  regions <- result$outputs$regions$table
  stopifnot(identical(summary$permutation_mode_used, "blocked"),
    isTRUE(summary$sample_blocks_available),
    identical(as.integer(summary$n_pseudotime_finite), 110L),
    identical(as.integer(summary$n_pseudotime_missing), 10L),
    identical(adjustment$adjustment_method, "linear"),
    adjustment$adjusted_max_abs_spearman <
      adjustment$raw_max_abs_spearman,
    !grepl("pseudotime", adjustment$adjustment_model, fixed = TRUE),
    any(regions$direction == "elevated" &
      abs(regions$peak_pseudotime - 0.42) <= 0.06),
    all(regions$n_samples >= 2L))

  permute <- getFromNamespace(".trajectory_regions_permute",
    "scwbTrajectory")
  set.seed(44L)
  values <- seq_len(12L)
  blocks <- rep(letters[1:3], each = 4L)
  permuted <- permute(values, blocks)
  stopifnot(all(vapply(split(seq_along(values), blocks), function(index)
    setequal(values[index], permuted[index]), logical(1))))
  TRUE
}

test_trajectory_signal_regions_registered_lifecycle_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  time <- seq(0, 1, length.out = length(x$cells))
  signals <- cbind(peak = 7 * exp(-((time - 0.30) / 0.10)^2),
    trough = -7 * exp(-((time - 0.70) / 0.10)^2))
  signal <- scworkbench::wb_prim_cell_vector_family(
    rep(colnames(signals), each = length(x$cells)),
    rep(x$cells, ncol(signals)), as.numeric(signals),
    value_type = "continuous")
  pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("L1", length(x$cells)), x$cells, time,
    value_type = "continuous")
  signal_id <- .trajectory_regions_promote(x$bundle, signal, "signal")
  pseudotime_id <- .trajectory_regions_promote(x$bundle, pseudotime,
    "pseudotime")
  samples <- scworkbench::wb_prim_cell_vector(stats::setNames(
    rep(paste0("sample_", 1:4), length.out = length(x$cells)), x$cells),
    value_type = "categorical")
  before <- .trajectory_regions_state(x$bundle)
  result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbtrajectory::trajectory_signal_regions", inputs = list(
      signal = scworkbench::wb_ref_artifact(signal_id, "vector"),
      pseudotime = scworkbench::wb_ref_artifact(pseudotime_id,
        "vector"), sample_id = samples, pseudotime_key = "L1",
      window_width_fraction = 0.30, window_step_fraction = 0.05,
      min_cells_per_window = 3L, min_samples_per_window = 2L,
      n_permutations = 39L, seed = 29L,
      adjust_depth = FALSE,
      dynamic_q_cutoff = 0.2, region_q_cutoff = 0.2,
      min_cells = 1L))
  outputs <- result$run$outputs
  stopifnot(identical(names(outputs), c("signal_table", "dynamic_tests",
      "candidate_regions", "regions", "adjustment_summary",
      "fit_summary")),
    identical(as.character(outputs$fit_summary$table$signal_artifact_id),
      signal_id),
    identical(as.character(
      outputs$fit_summary$table$pseudotime_artifact_id), pseudotime_id),
    identical(before, .trajectory_regions_state(x$bundle)))

  plot_data <- scworkbench::wb_wrangle_plot_dataset(
    "pseudotime_signature_line", inputs = list(table =
      outputs$signal_table), bundle = x$bundle)
  stopifnot(is.data.frame(plot_data$table), nrow(plot_data$table) > 0L,
    identical(plot_data$dataset_type,
      "trajectory_signal_dynamics_curve"))

  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    result$session, outputs = names(outputs),
    created_by = "signal-regions-contract")
  finalized <- scworkbench::wb_finalize_recipient_bundle(x$bundle)
  scworkbench::wb_close_bundle(x$bundle)
  x$bundle <- NULL
  recipient <- scworkbench::wb_bundle_validation_report(x$root,
    "recipient")
  x$bundle <- scworkbench::wb_open_bundle(x$root, read_only = TRUE)
  stopifnot(setequal(promoted$promoted$output, names(outputs)),
    all(promoted$promoted$status == "ok"),
    isTRUE(finalized$ready),
    isTRUE(recipient$ok))
  TRUE
}

test_trajectory_signal_regions_rejection_no_state_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  time <- seq(0, 1, length.out = length(x$cells))
  signal <- scworkbench::wb_prim_cell_vector_family(
    rep("signal", length(x$cells)), x$cells, sin(2 * pi * time))
  pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("L1", length(x$cells)), x$cells, time)
  before <- .trajectory_regions_state(x$bundle)
  run <- function(extra = list()) scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_signal_regions", inputs = c(list(
      signal = signal, pseudotime = pseudotime,
      pseudotime_key = "L1", window_width_fraction = 0.3,
      adjust_depth = FALSE,
      window_step_fraction = 0.1, min_cells_per_window = 3L,
      n_permutations = 9L, min_cells = 1L), extra))
  missing_key <- try(run(list(pseudotime_key = "missing")), silent = TRUE)
  invalid_width <- try(run(list(window_width_fraction = -0.1)),
    silent = TRUE)
  blocked_without_samples <- try(run(list(permutation_mode = "blocked")),
    silent = TRUE)
  duplicate <- signal
  duplicate$cell_id[[2L]] <- duplicate$cell_id[[1L]]
  duplicate_signal <- try(run(list(signal = duplicate)), silent = TRUE)
  stopifnot(inherits(missing_key, "try-error"),
    inherits(invalid_width, "try-error"),
    inherits(blocked_without_samples, "try-error"),
    inherits(duplicate_signal, "try-error"),
    identical(before, .trajectory_regions_state(x$bundle)))
  TRUE
}
