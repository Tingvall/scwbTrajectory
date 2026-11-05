.depth_adjustment_synthetic <- function(n = 120L) {
  cells <- sprintf("depth_cell_%03d", seq_len(n))
  pseudotime <- seq(0, 1, length.out = n)
  depth <- exp(6 + 0.8 * pseudotime + 0.15 * sin(7 * pseudotime))
  centered <- log1p(depth) - mean(log1p(depth))
  basis <- splines::ns(pseudotime, df = 4L)
  biological <- as.numeric(1 + basis %*% c(0.7, -0.4, 1.1, 0.3))
  values <- biological + 2.5 * centered
  list(cells = cells, pseudotime = pseudotime, depth = depth,
    biological = biological, values = values)
}

test_trajectory_depth_partial_component_contract <- function() {
  x <- .depth_adjustment_synthetic()
  matrix <- rbind(gene_depth = x$values,
    gene_constant = rep(3, length(x$cells)))
  colnames(matrix) <- x$cells
  method <- getFromNamespace(".trajectory_depth_adjust_matrix",
    "scwbTrajectory")
  adjusted <- method(matrix, x$pseudotime, x$depth, spline_df = 4L,
    standardize = FALSE)
  stopifnot(abs(adjusted$diagnostics$depth_coefficient[[1L]] - 2.5) < 1e-8,
    max(abs(adjusted$values[1L, ] - x$biological)) < 1e-8,
    identical(adjusted$diagnostics$status[[2L]],
      "zero_adjusted_variance"))
  standardized <- method(matrix, x$pseudotime, x$depth,
    spline_df = 4L, standardize = TRUE)
  stopifnot(abs(mean(standardized$values[1L, ])) < 1e-12,
    abs(stats::sd(standardized$values[1L, ]) - 1) < 1e-12,
    identical(as.numeric(standardized$values[2L, ]),
      rep(0, length(x$cells))))
  TRUE
}

test_trajectory_geneset_depth_scoring_contract <- function() {
  x <- .depth_adjustment_synthetic(60L)
  matrix <- rbind(gene_a = x$values,
    gene_b = x$biological - 1.25 *
      (log1p(x$depth) - mean(log1p(x$depth))))
  colnames(matrix) <- x$cells
  pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("lineage", length(x$cells)), x$cells, x$pseudotime,
    value_type = "continuous")
  registry <- list(sets = list(list(namespace = "synthetic",
    feature_set = "depth_test", memberships = list(
      list(feature_id = "gene_a"), list(feature_id = "gene_b")))))
  method <- getFromNamespace(".trajectory_s2_geneset_dynamics",
    "scwbTrajectory")
  common <- list(matrix = matrix, registry = registry,
    effective_cells = x$cells, matrix_id = "SCT_data",
    pseudotime = pseudotime, pseudotime_key = "lineage",
    namespace = "synthetic", min_set_size = 1L, n_bins = 10L,
    min_cells_per_bin = 1L, top_driver_genes = 1L)
  adjusted <- do.call(method, c(common, list(depth = stats::setNames(
    x$depth, x$cells), adjust_depth = TRUE)))
  stopifnot(identical(names(adjusted$outputs), c("signature_table",
      "dynamics_table", "top_sets", "driver_table", "driver_curve_table",
      "adjusted_scores", "adjustment_diagnostics", "fit_summary")),
    isTRUE(adjusted$outputs$adjusted_scores$meta$depth_adjusted),
    identical(adjusted$outputs$adjusted_scores$meta$
      effective_score_method, "adjusted_z_mean"),
    nrow(adjusted$outputs$adjustment_diagnostics$table) == 2L,
    max(abs(adjusted$outputs$adjustment_diagnostics$table$
      depth_coefficient - c(2.5, -1.25))) < 1e-8)
  interval_method <- getFromNamespace(".trajectory_signal_regions",
    "scwbTrajectory")
  intervals <- interval_method(effective_cells = x$cells,
    signal = adjusted$outputs$adjusted_scores, pseudotime = pseudotime,
    pseudotime_key = "lineage", signal_keys = "lineage::depth_test",
    adjust_depth = TRUE, window_width_fraction = 0.2,
    window_step_fraction = 0.05, min_cells_per_window = 6L,
    n_permutations = 19L, seed = 53L)
  stopifnot(identical(intervals$outputs$adjustment_summary$table$
    depth_adjustment_status, "already_adjusted"))
  raw <- do.call(method, c(common, list(adjust_depth = FALSE)))
  stopifnot(max(abs(as.numeric(raw$outputs$adjusted_scores$value) -
      colMeans(matrix))) < 1e-12,
    !isTRUE(raw$outputs$adjusted_scores$meta$depth_adjusted),
    nrow(raw$outputs$adjustment_diagnostics$table) == 0L)
  TRUE
}

test_trajectory_signal_interval_depth_and_reuse_contract <- function() {
  x <- .depth_adjustment_synthetic()
  signal <- scworkbench::wb_prim_cell_vector_family(
    rep("signal", length(x$cells)), x$cells, x$values,
    value_type = "continuous")
  pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("lineage", length(x$cells)), x$cells, x$pseudotime,
    value_type = "continuous")
  depth <- scworkbench::wb_prim_cell_vector(stats::setNames(
    x$depth, x$cells), value_type = "continuous")
  method <- getFromNamespace(".trajectory_signal_regions",
    "scwbTrajectory")
  applied <- method(effective_cells = x$cells, signal = signal,
    pseudotime = pseudotime, pseudotime_key = "lineage", depth = depth,
    adjust_depth = TRUE, depth_field = "nCount_RNA",
    window_width_fraction = 0.15, window_step_fraction = 0.05,
    min_cells_per_window = 8L, n_permutations = 19L, seed = 71L)
  summary <- applied$outputs$adjustment_summary$table
  stopifnot(identical(summary$depth_adjustment_status, "applied"),
    abs(summary$depth_coefficient - 2.5) < 1e-8,
    identical(applied$outputs$fit_summary$table$depth_field,
      "nCount_RNA"))

  adjusted_signal <- scworkbench::wb_prim_cell_vector_family(
    rep("signal", length(x$cells)), x$cells, x$biological,
    value_type = "continuous", meta = list(depth_adjusted = TRUE,
      depth_field = "nCount_RNA", adjustment_spline_df = 4L))
  reused <- method(effective_cells = x$cells, signal = adjusted_signal,
    pseudotime = pseudotime, pseudotime_key = "lineage",
    adjust_depth = TRUE, depth_field = "nCount_RNA",
    window_width_fraction = 0.15, window_step_fraction = 0.05,
    min_cells_per_window = 8L, n_permutations = 19L, seed = 71L)
  stopifnot(identical(reused$outputs$adjustment_summary$table$
      depth_adjustment_status, "already_adjusted"),
    identical(reused$outputs$fit_summary$table$depth_adjustment_status,
      "already_adjusted"))
  TRUE
}

test_trajectory_depth_manifest_defaults_contract <- function() {
  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbTrajectory"),
    simplifyVector = FALSE)
  operations <- stats::setNames(manifest$operations,
    vapply(manifest$operations, `[[`, character(1), "name"))
  for (name in c("trajectory_geneset_dynamics",
      "trajectory_geneset_dynamics_explicit")) {
    operation <- operations[[name]]
    stopifnot(identical(operation$params_schema$properties$
        adjust_depth$default, TRUE),
      identical(operation$params_schema$properties$depth_field$default,
        "nCount_RNA"),
      identical(vapply(utils::tail(operation$outputs, 3L), `[[`,
        character(1), "name"), c("adjusted_scores",
          "adjustment_diagnostics", "fit_summary")),
      "depth" %in% vapply(operation$inputs, `[[`, character(1), "name"),
      any(vapply(operation$provenance_reads, function(read)
        identical(read$name, "depth_ancestry"), logical(1))))
  }
  stopifnot(identical(operations$trajectory_signal_intervals$
      params_schema$properties$adjust_depth$default, TRUE),
    identical(operations$trajectory_signal_regions$params_schema$properties$
      adjust_depth$default, TRUE),
    identical(operations$trajectory_preranked_gsea$params_schema$properties$
      adjust_depth$default, FALSE))
  TRUE
}

test_trajectory_geneset_depth_registered_lifecycle_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("L1", length(x$cells)), x$cells,
    seq(0, 1, length.out = length(x$cells)), value_type = "continuous")
  depth <- scworkbench::wb_prim_cell_vector(stats::setNames(
    as.numeric(colSums(x$values)), x$cells), value_type = "continuous")
  result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbtrajectory::trajectory_geneset_dynamics_explicit",
    inputs = list(matrix = "RNA", pseudotime = pseudotime, depth = depth,
      namespace = "trajectory", pseudotime_key = "L1",
      min_set_size = 1L, n_bins = 6L, min_cells_per_bin = 1L,
      min_cells = 1L))
  outputs <- result$run$outputs
  stopifnot(isTRUE(outputs$adjusted_scores$meta$depth_adjusted),
    identical(outputs$adjusted_scores$meta$depth_field, "nCount_RNA"),
    nrow(outputs$adjustment_diagnostics$table) == 2L)
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    result$session, outputs = names(outputs), created_by = "depth-contract")
  stopifnot(nrow(promoted$promoted) == length(outputs),
    all(promoted$promoted$status == "ok"))
  TRUE
}

test_trajectory_signal_interval_depth_registered_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  pseudotime_value <- seq(0, 1, length.out = length(x$cells))
  depth_value <- as.numeric(colSums(x$values))
  centered <- log1p(depth_value) - mean(log1p(depth_value))
  signal <- scworkbench::wb_prim_cell_vector_family(
    rep("peak", length(x$cells)), x$cells,
    6 * exp(-((pseudotime_value - 0.45) / 0.15)^2) + centered,
    value_type = "continuous")
  pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("L1", length(x$cells)), x$cells, pseudotime_value,
    value_type = "continuous")
  depth <- scworkbench::wb_prim_cell_vector(stats::setNames(
    depth_value, x$cells), value_type = "continuous")
  result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbtrajectory::trajectory_signal_intervals", inputs = list(
      signal = signal, pseudotime = pseudotime, depth = depth,
      pseudotime_key = "L1", window_width_fraction = 0.3,
      window_step_fraction = 0.05, min_cells_per_window = 3L,
      n_permutations = 19L, seed = 83L, global_q_cutoff = 0.2,
      local_q_cutoff = 0.2, min_cells = 1L))
  stopifnot(identical(result$run$outputs$adjustment_summary$table$
      depth_adjustment_status, "applied"),
    isTRUE(result$run$outputs$fit_summary$table$adjust_depth),
    identical(result$run$outputs$fit_summary$table$depth_field,
      "nCount_RNA"))
  TRUE
}

test_trajectory_depth_rejection_contract <- function() {
  x <- .depth_adjustment_synthetic(40L)
  method <- getFromNamespace(".trajectory_depth_adjust_matrix",
    "scwbTrajectory")
  matrix <- matrix(x$values, nrow = 1L,
    dimnames = list("gene", x$cells))
  negative <- x$depth; negative[[3L]] <- -1
  missing <- x$depth; missing[[4L]] <- NA_real_
  stopifnot(inherits(try(method(matrix, x$pseudotime, negative),
      silent = TRUE), "try-error"),
    inherits(try(method(matrix, x$pseudotime, missing), silent = TRUE),
      "try-error"))
  TRUE
}
