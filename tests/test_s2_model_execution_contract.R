test_dir <- if (file.exists("helper_s2_fixture.R")) "." else "tests"
source(file.path(test_dir, "helper_s2_fixture.R"), local = .GlobalEnv)

.trajectory_s2_model_params <- function(operation) {
  switch(operation,
    de = list(lineage = "L1", features = c("g1", "g2", "g3"),
      chunk_size = 2L),
    smooth = list(lineage = "L1", features = c("g1", "g2", "g3"), n_bins = 5L,
      min_cells_per_bin = 1L, chunk_size = 2L),
    signature_smooth = list(feature_set = "selected", lineage = "L1",
      n_bins = 5L, min_cells_per_bin = 1L, chunk_size = 2L),
    vector_smooth = list(lineage = "L1", signal_keys = "rising",
      n_bins = 5L, min_cells_per_bin = 1L),
    geneset_dynamics = list(namespace = "trajectory", lineage = "L1",
      feature_sets = "selected",
      adjust_depth = FALSE,
      min_set_size = 1L, n_bins = 5L, min_cells_per_bin = 1L,
      min_bins_observed = 2L, bootstrap_n = 8L, bootstrap_seed = 17L,
      n_permutations = 8L, permutation_seed = 19L, chunk_size = 2L),
    signal_dynamics = list(lineage = "L1", n_bins = 5L,
      min_cells_per_bin = 1L, min_bins_observed = 2L,
      bootstrap_n = 8L, bootstrap_seed = 17L))
}

.trajectory_s2_synthetic_slingshot <- function(cells) {
  stopifnot(requireNamespace("slingshot", quietly = TRUE))
  coordinates <- matrix(seq_len(length(cells) * 2L),
    nrow = length(cells), ncol = 2L,
    dimnames = list(cells, c("x", "y")))
  model <- slingshot::newSlingshotDataSet(coordinates,
    rep("A", length(cells)))
  curve <- function(lambda) structure(list(s = coordinates,
    ord = seq_along(cells), lambda = lambda,
    dist = rep(0, length(cells)), w = rep(1, length(cells))),
    class = "principal_curve")
  methods::slot(model, "curves") <- list(
    L1 = curve(seq(0, 1, length.out = length(cells))),
    L2 = curve(seq(1, 0, length.out = length(cells))))
  stopifnot(methods::validObject(model))
  model
}

.trajectory_s2_run_model_direct <- function(x, operation, params = NULL,
                                            model = NULL,
                                            membership = NULL,
                                            fault_terminal = NULL,
                                            fitted_cells = x$cells,
                                            fitted_source_space = NULL,
                                            include_fitted_space = TRUE,
                                            effective_cells = x$cells,
                                            model_payload_reader = NULL) {
  boundary <- trajectory_s2_mock_boundary(x, fault_terminal)
  boundary$readers$membership <- if (is.null(membership)) NULL else
    function() membership
  boundary$context$effective_cells <- as.character(effective_cells)
  boundary$readers$matrix_chunks <- function(features, cells) {
    stopifnot(identical(as.character(cells), as.character(effective_cells)),
      all(features %in% x$features))
    boundary$state$calls$matrix_chunks <- c(
      boundary$state$calls$matrix_chunks, as.character(features))
    list(values = x$values[features, cells, drop = FALSE])
  }
  if (is.null(fitted_source_space)) fitted_source_space <- matrix(
    seq_len(length(fitted_cells) * 2L), ncol = 2L,
    dimnames = list(as.character(fitted_cells), c("PC1", "PC2")))
  path <- tempfile("s2_g1_model_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  if (is.null(model)) model <- .trajectory_s2_synthetic_slingshot(x$cells)
  saveRDS(model, path)
  boundary$readers$model_payload <- function() {
    boundary$state$calls$model_payload <- as.integer(
      boundary$state$calls$model_payload %||% 0L) + 1L
    stopifnot(identical(boundary$state$calls$model_payload, 1L))
    c(list(primary_path = "model.rds"),
      if (isTRUE(include_fitted_space))
        list(fitted_source_space = fitted_source_space) else list(), list(
      runtime_private_path = function(entry) {
        stopifnot(identical(entry, "model.rds"))
        path
      }))
  }
  if (!is.null(model_payload_reader))
    boundary$readers$model_payload <- model_payload_reader
  name <- paste0(".run_trajectory_", operation, "_model")
  runner <- getFromNamespace(name, "scwbTrajectory")
  if (is.null(params)) params <- .trajectory_s2_model_params(operation)
  value <- tryCatch(runner(boundary$readers, boundary$writers,
    params, boundary$context), error = identity)
  list(value = value, outputs = boundary$state$writes,
    calls = boundary$state$calls)
}

test_trajectory_s2_model_direct_oracle_contract <- function() {
  if (!requireNamespace("slingshot", quietly = TRUE))
    return(structure(TRUE, skip = "slingshot unavailable"))
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  reference <- trajectory_s2_reference()
  mapping <- c(de = "trajectory_de", smooth = "trajectory_smooth",
    signature_smooth = "trajectory_signature_smooth",
    vector_smooth = "trajectory_vector_smooth",
    geneset_dynamics = "trajectory_geneset_dynamics_model",
    signal_dynamics = "trajectory_signal_dynamics_model")
  set.seed(20260923)
  seed <- get(".Random.seed", envir = .GlobalEnv)
  for (operation in names(mapping)) {
    observed <- .trajectory_s2_run_model_direct(x, operation)
    stopifnot(!inherits(observed$value, "error"),
      identical(observed$calls$model_payload, 1L))
    trajectory_s2_compare_outputs(x$bundle, mapping[[operation]],
      list(outputs = observed$outputs), reference)
    if (operation %in% c("de", "smooth", "signature_smooth",
      "geneset_dynamics")) stopifnot(identical(
        observed$calls$matrix_chunks, x$features))
    if (operation %in% c("geneset_dynamics", "signal_dynamics")) {
      repeated <- .trajectory_s2_run_model_direct(x, operation)
      stopifnot(identical(observed$outputs, repeated$outputs))
    }
  }
  subset <- x$cells[1:12]
  observed <- .trajectory_s2_run_model_direct(x, "smooth",
    effective_cells = subset)
  parameter_reference <- jsonlite::read_json(file.path(test_dir,
    "fixtures", "trajectory_s2_parameter_reference.json"),
    simplifyVector = FALSE)
  stopifnot(!inherits(observed$value, "error"))
  trajectory_s2_compare_outputs(x$bundle, "smooth_subset",
    list(outputs = observed$outputs), parameter_reference)
  stopifnot(identical(seed, get(".Random.seed", envir = .GlobalEnv)))
  TRUE
}

test_trajectory_s2_model_direct_rejection_contract <- function() {
  if (!requireNamespace("slingshot", quietly = TRUE))
    return(structure(TRUE, skip = "slingshot unavailable"))
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  before <- wb_art_list(x$bundle)
  for (operation in c("de", "smooth", "signature_smooth",
      "vector_smooth", "geneset_dynamics", "signal_dynamics")) {
    invalid <- .trajectory_s2_run_model_direct(x, operation,
      model = list(pseudotime_table = x$model$pseudotime_table))
    stopifnot(inherits(invalid$value, "error"),
      grepl("authenticated lineage projection",
        conditionMessage(invalid$value), fixed = TRUE),
      length(invalid$outputs) == 0L)
    wrong <- .trajectory_s2_model_params(operation)
    wrong$lineage <- "missing"
    if (identical(operation, "signature_smooth"))
      wrong$feature_set <- "selected"
    rejected <- .trajectory_s2_run_model_direct(x, operation, wrong)
    stopifnot(inherits(rejected$value, "error"),
      length(rejected$outputs) == 0L)
  }
  for (operation in c("de", "smooth", "signature_smooth")) {
    wrong <- .trajectory_s2_model_params(operation)
    wrong$feature_set <- "missing"
    rejected <- .trajectory_s2_run_model_direct(x, operation, wrong)
    stopifnot(inherits(rejected$value, "error"),
      length(rejected$outputs) == 0L)
  }
  wrong <- .trajectory_s2_model_params("vector_smooth")
  wrong$signal_keys <- c("rising", "missing")
  rejected <- .trajectory_s2_run_model_direct(x, "vector_smooth", wrong)
  stopifnot(inherits(rejected$value, "error"),
    length(rejected$outputs) == 0L)
  foreign <- .trajectory_s2_synthetic_slingshot(
    paste0("foreign_", seq_along(x$cells)))
  for (operation in c("de", "signal_dynamics")) {
    rejected <- .trajectory_s2_run_model_direct(x, operation,
      model = foreign)
    stopifnot(inherits(rejected$value, "error"),
      length(rejected$outputs) == 0L)
  }
  reordered <- .trajectory_s2_run_model_direct(x, "de",
    fitted_cells = rev(x$cells))
  stopifnot(inherits(reordered$value, "error"),
    grepl("differ from the authenticated fitted axis",
      conditionMessage(reordered$value), fixed = TRUE),
    length(reordered$outputs) == 0L)
  subset_model <- .trajectory_s2_synthetic_slingshot(x$cells[1:12])
  different_subset <- .trajectory_s2_run_model_direct(x,
    "signal_dynamics", model = subset_model,
    fitted_cells = x$cells[2:13])
  stopifnot(inherits(different_subset$value, "error"),
    grepl("differ from the authenticated fitted axis",
      conditionMessage(different_subset$value), fixed = TRUE),
    length(different_subset$outputs) == 0L)
  missing_fitted <- .trajectory_s2_run_model_direct(x, "de",
    include_fitted_space = FALSE)
  malformed_fitted <- .trajectory_s2_run_model_direct(x, "de",
    fitted_source_space = data.frame(PC1 = seq_along(x$cells),
      row.names = x$cells))
  for (rejected in list(missing_fitted, malformed_fitted))
    stopifnot(inherits(rejected$value, "error"),
      grepl("no valid authenticated fitted source space",
        conditionMessage(rejected$value), fixed = TRUE),
      length(rejected$outputs) == 0L)
  fault <- .trajectory_s2_run_model_direct(x, "smooth",
    fault_terminal = "fit_summary")
  stopifnot(inherits(fault$value, "error"),
    grepl("injected terminal fault", conditionMessage(fault$value),
      fixed = TRUE),
    length(fault$outputs) == 1L,
    identical(before, wb_art_list(x$bundle)))
  TRUE
}

test_trajectory_s2_model_parameter_equivalence_contract <- function() {
  if (!requireNamespace("slingshot", quietly = TRUE))
    return(structure(TRUE, skip = "slingshot unavailable"))
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  reference <- jsonlite::read_json(file.path(test_dir, "fixtures",
    "trajectory_s2_parameter_reference.json"), simplifyVector = FALSE)
  stopifnot(identical(reference$core_commit,
      "ba24635611756240f21f6817044883b832f12701"),
    identical(reference$comparison_policy, "exact_canonical_value"))
  variants <- list(
    de_variant = list(operation = "de", params = list(lineage = "L2",
      features = x$features[1:3], min_sd = 0.1,
      p_adjust_method = "bonferroni", chunk_size = 2L)),
    smooth_variant = list(operation = "smooth", params = list(
      lineage = "L2", features = x$features[1:3], n_bins = 5L,
      summary_fun = "median", min_cells_per_bin = 2L, chunk_size = 2L)),
    signature_variant = list(operation = "signature_smooth", params = list(
      feature_set = "selected", lineage = "L2", score_method = "sum",
      n_bins = 5L, min_cells_per_bin = 2L, chunk_size = 2L)),
    vector_variant = list(operation = "vector_smooth", params = list(
      lineage = "L2", signal_keys = "falling", n_bins = 5L,
      summary_fun = "median", min_cells_per_bin = 2L)),
    de_reordered_lineages = list(operation = "de", params = list(
      lineage = c("L2", "L1"), features = x$features[1:3],
      chunk_size = 2L)),
    de_feature_intersection = list(operation = "de", params = list(
      feature_set = "selected", features = c("g1", "g3"),
      lineage = "L1", chunk_size = 2L)),
    smooth_reordered_lineages = list(operation = "smooth", params = list(
      lineage = c("L2", "L1"), features = x$features[1:3],
      n_bins = 5L, min_cells_per_bin = 1L, chunk_size = 2L)),
    geneset_membership = list(operation = "geneset_dynamics",
      params = .trajectory_s2_model_params("geneset_dynamics"),
      membership = x$membership),
    signal_membership = list(operation = "signal_dynamics",
      params = .trajectory_s2_model_params("signal_dynamics"),
      membership = x$membership))
  for (case_name in names(variants)) {
    case <- variants[[case_name]]
    observed <- .trajectory_s2_run_model_direct(x, case$operation,
      params = case$params, membership = case$membership %||% NULL)
    stopifnot(!inherits(observed$value, "error"))
    trajectory_s2_compare_outputs(x$bundle, case_name,
      list(outputs = observed$outputs), reference)
    if (identical(case_name, "de_feature_intersection"))
      stopifnot(identical(unique(as.character(observed$outputs$
        de_table$table$feature_id)), "g1"))
    expected_params <- reference$cases[[case_name]]$params
    for (terminal in observed$outputs) {
      common <- intersect(names(terminal$meta), names(expected_params))
      for (key in common) {
        want <- expected_params[[key]]
        same <- if (is.null(want)) {
          length(terminal$meta[[key]]) == 1L &&
            is.na(terminal$meta[[key]])
        } else isTRUE(all.equal(terminal$meta[[key]],
          unlist(want, use.names = FALSE), check.attributes = FALSE))
        if (!same) stop(case_name, "/", terminal$table_type,
          ": effective parameter differs: ", key, call. = FALSE)
      }
    }
  }
  base <- .trajectory_s2_run_model_direct(x, "geneset_dynamics")$outputs
  altered <- .trajectory_s2_model_params("geneset_dynamics")
  altered$bootstrap_seed <- 18L
  changed_ci <- .trajectory_s2_run_model_direct(x, "geneset_dynamics",
    altered)$outputs
  altered <- .trajectory_s2_model_params("geneset_dynamics")
  altered$permutation_seed <- 20L
  changed_null <- .trajectory_s2_run_model_direct(x,
    "geneset_dynamics", altered)$outputs
  signal <- .trajectory_s2_run_model_direct(x,
    "signal_dynamics")$outputs
  altered <- .trajectory_s2_model_params("signal_dynamics")
  altered$bootstrap_seed <- 18L
  signal_changed <- .trajectory_s2_run_model_direct(x,
    "signal_dynamics", altered)$outputs
  stopifnot(!identical(base$signature_table$table$value_ci_low,
      changed_ci$signature_table$table$value_ci_low),
    !identical(base$dynamics_table$table$null_mean_dynamic_score,
      changed_null$dynamics_table$table$null_mean_dynamic_score),
    !identical(signal$signal_table$table$value_ci_low,
      signal_changed$signal_table$table$value_ci_low))
  TRUE
}
