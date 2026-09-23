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
                                            fault_terminal = NULL) {
  boundary <- trajectory_s2_mock_boundary(x, fault_terminal)
  boundary$readers$membership <- if (is.null(membership)) NULL else
    function() membership
  path <- tempfile("s2_g1_model_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  if (is.null(model)) model <- .trajectory_s2_synthetic_slingshot(x$cells)
  saveRDS(model, path)
  boundary$readers$model_payload <- function() {
    boundary$state$calls$model_payload <- as.integer(
      boundary$state$calls$model_payload %||% 0L) + 1L
    stopifnot(identical(boundary$state$calls$model_payload, 1L))
    list(primary_path = "model.rds",
      runtime_private_path = function(entry) {
        stopifnot(identical(entry, "model.rds"))
        path
      })
  }
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
  variants <- list(
    de = list(params = list(lineage = "L2", features = x$features[1:3],
      min_sd = 0.1, p_adjust_method = "bonferroni", chunk_size = 2L),
      core = function() wb_op_trajectory_de(x$bundle, "RNA", x$model,
        lineage = "L2", features = x$features[1:3], min_sd = 0.1,
        p_adjust_method = "bonferroni", chunk_size = 2L)),
    smooth = list(params = list(lineage = "L2", features = x$features[1:3],
      n_bins = 5L, summary_fun = "median", min_cells_per_bin = 2L,
      chunk_size = 2L),
      core = function() wb_op_trajectory_smooth(x$bundle, "RNA", x$model,
        lineage = "L2", features = x$features[1:3], n_bins = 5L,
        summary_fun = "median", min_cells_per_bin = 2L, chunk_size = 2L)),
    signature_smooth = list(params = list(feature_set = "selected",
      lineage = "L2", score_method = "sum", n_bins = 5L,
      min_cells_per_bin = 2L, chunk_size = 2L),
      core = function() wb_op_trajectory_signature_smooth(x$bundle,
        "RNA", x$model, feature_set = "selected", lineage = "L2",
        score_method = "sum", n_bins = 5L, min_cells_per_bin = 2L)),
    vector_smooth = list(params = list(lineage = "L2",
      signal_keys = "falling", n_bins = 5L, summary_fun = "median",
      min_cells_per_bin = 2L),
      core = function() wb_op_trajectory_vector_smooth(x$bundle,
        x$signal, x$model, lineage = "L2", signal_keys = "falling",
        n_bins = 5L, summary_fun = "median", min_cells_per_bin = 2L)))
  for (operation in names(variants)) {
    observed <- .trajectory_s2_run_model_direct(x, operation,
      variants[[operation]]$params)
    stopifnot(!inherits(observed$value, "error"))
    expected <- variants[[operation]]$core()$outputs
    stopifnot(identical(names(observed$outputs), names(expected)))
    for (name in names(expected)) {
      actual_identity <- trajectory_identity(x$bundle,
        observed$outputs[[name]])
      expected_identity <- trajectory_identity(x$bundle, expected[[name]])
      stopifnot(identical(actual_identity$content_sha256,
        expected_identity$content_sha256),
        identical(actual_identity$structural_fingerprint,
          expected_identity$structural_fingerprint))
    }
  }
  for (operation in c("de", "smooth")) {
    params <- .trajectory_s2_model_params(operation)
    params$lineage <- c("L2", "L1")
    observed <- .trajectory_s2_run_model_direct(x, operation, params)
    stopifnot(!inherits(observed$value, "error"))
    expected <- if (identical(operation, "de"))
      wb_op_trajectory_de(x$bundle, "RNA", x$model,
        lineage = c("L2", "L1"), features = x$features[1:3],
        chunk_size = 2L)$outputs else
      wb_op_trajectory_smooth(x$bundle, "RNA", x$model,
        lineage = c("L2", "L1"), features = x$features[1:3],
        n_bins = 5L, min_cells_per_bin = 1L,
        chunk_size = 2L)$outputs
    for (name in names(expected)) stopifnot(identical(
      trajectory_identity(x$bundle, observed$outputs[[name]])$content_sha256,
      trajectory_identity(x$bundle, expected[[name]])$content_sha256))
  }
  for (operation in c("geneset_dynamics", "signal_dynamics")) {
    observed <- .trajectory_s2_run_model_direct(x, operation,
      membership = x$membership)
    stopifnot(!inherits(observed$value, "error"))
    expected <- if (identical(operation, "geneset_dynamics"))
      wb_op_trajectory_geneset_dynamics(x$bundle, "RNA", model = x$model,
        membership = x$membership, namespace = "trajectory",
        lineage = "L1", min_set_size = 1L, n_bins = 5L,
        min_cells_per_bin = 1L, min_bins_observed = 2L,
        bootstrap_n = 8L, bootstrap_seed = 17L,
        n_permutations = 8L, permutation_seed = 19L)$outputs else
      wb_op_trajectory_signal_dynamics(x$bundle, x$signal,
        model = x$model, membership = x$membership,
        lineage = "L1", n_bins = 5L, min_cells_per_bin = 1L,
        min_bins_observed = 2L, bootstrap_n = 8L,
        bootstrap_seed = 17L)$outputs
    for (name in names(expected)) stopifnot(identical(
      trajectory_identity(x$bundle, observed$outputs[[name]])$content_sha256,
      trajectory_identity(x$bundle, expected[[name]])$content_sha256))
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
