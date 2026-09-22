source(file.path(test_dir, "helper_s2_fixture.R"), local = .GlobalEnv)

.trajectory_s2_explicit_params <- function(operation) {
  if (identical(operation, "geneset")) return(list(
    namespace = "trajectory", pseudotime_key = "L1", min_set_size = 1L,
    n_bins = 5L, min_cells_per_bin = 1L, min_bins_observed = 2L,
    bootstrap_n = 8L, bootstrap_seed = 17L, n_permutations = 8L,
    permutation_seed = 19L, chunk_size = 2L))
  list(pseudotime_key = "L1", n_bins = 5L, min_cells_per_bin = 1L,
    min_bins_observed = 2L, bootstrap_n = 8L, bootstrap_seed = 17L)
}

.trajectory_s2_run_direct <- function(x, operation, params = NULL,
                                      fault_terminal = NULL,
                                      membership = NULL) {
  boundary <- trajectory_s2_mock_boundary(x, fault_terminal)
  if (!is.null(membership)) {
    boundary$readers$membership <- function() membership
  }
  runner <- getFromNamespace(paste0(".run_trajectory_", operation,
    "_dynamics_explicit"), "scwbTrajectory")
  if (is.null(params)) params <- .trajectory_s2_explicit_params(operation)
  value <- tryCatch(runner(boundary$readers, boundary$writers, params,
    boundary$context), error = identity)
  list(value = value, boundary = boundary,
    outputs = boundary$state$writes)
}

test_trajectory_s2_explicit_runners_oracle_and_seed_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  stopifnot(!any(c("scwbtrajectory::trajectory_geneset_dynamics",
    "scwbtrajectory::trajectory_signal_dynamics") %in% wb_list_ops()))
  reference <- trajectory_s2_reference()
  set.seed(20260922)
  seed_before <- get(".Random.seed", envir = .GlobalEnv)
  for (operation in c("geneset", "signal")) {
    first <- .trajectory_s2_run_direct(x, operation)
    second <- .trajectory_s2_run_direct(x, operation)
    name <- paste0("trajectory_", operation, "_dynamics_explicit")
    stopifnot(!inherits(first$value, "error"),
      !inherits(second$value, "error"))
    trajectory_s2_compare_outputs(x$bundle, name,
      list(outputs = first$outputs), reference)
    trajectory_s2_compare_outputs(x$bundle, name,
      list(outputs = second$outputs), reference)
    stopifnot(identical(first$outputs, second$outputs),
      identical(first$boundary$state$calls$pseudotime, 1L),
      identical(first$boundary$state$calls$membership, 1L))
    if (identical(operation, "geneset")) {
      stopifnot(identical(first$boundary$state$calls$feature_sets, 1L),
        identical(first$boundary$state$calls$matrix_chunks, x$features))
    } else {
      stopifnot(identical(first$boundary$state$calls$signal, 1L))
    }
  }
  baseline <- .trajectory_s2_run_direct(x, "geneset")$outputs
  changed_bootstrap <- .trajectory_s2_explicit_params("geneset")
  changed_bootstrap$bootstrap_seed <- 18L
  changed_ci <- .trajectory_s2_run_direct(x, "geneset",
    changed_bootstrap)$outputs
  changed_permutation <- .trajectory_s2_explicit_params("geneset")
  changed_permutation$permutation_seed <- 20L
  changed_null <- .trajectory_s2_run_direct(x, "geneset",
    changed_permutation)$outputs
  signal_baseline <- .trajectory_s2_run_direct(x, "signal")$outputs
  signal_seed <- .trajectory_s2_explicit_params("signal")
  signal_seed$bootstrap_seed <- 18L
  signal_changed <- .trajectory_s2_run_direct(x, "signal", signal_seed)$outputs
  stopifnot(!identical(baseline$signature_table$table$value_ci_low,
      changed_ci$signature_table$table$value_ci_low),
    !identical(baseline$dynamics_table$table$null_mean_dynamic_score,
      changed_null$dynamics_table$table$null_mean_dynamic_score),
    !identical(signal_baseline$signal_table$table$value_ci_low,
      signal_changed$signal_table$table$value_ci_low))
  stopifnot(identical(seed_before,
    get(".Random.seed", envir = .GlobalEnv)))
  rm(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", seed_before, envir = .GlobalEnv),
    add = TRUE)
  no_prior_seed <- .trajectory_s2_run_direct(x, "geneset")
  stopifnot(!inherits(no_prior_seed$value, "error"),
    !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  TRUE
}

test_trajectory_s2_explicit_runners_rejection_and_rollback_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  reference <- trajectory_s2_reference()
  before <- wb_art_list(x$bundle)
  few <- x$membership
  few$value[9:length(few$value)] <- NA_character_
  for (operation in c("geneset", "signal")) {
    no_key <- .trajectory_s2_explicit_params(operation)
    no_key$pseudotime_key <- NULL
    absent <- .trajectory_s2_run_direct(x, operation, no_key)
    stopifnot(inherits(absent$value, "error"),
      length(absent$outputs) == 0L)
    wrong <- .trajectory_s2_explicit_params(operation)
    wrong$pseudotime_key <- "missing"
    bad_key <- .trajectory_s2_run_direct(x, operation, wrong)
    stopifnot(inherits(bad_key$value, "error"),
      length(bad_key$outputs) == 0L)
    bad_selector <- .trajectory_s2_explicit_params(operation)
    if (identical(operation, "geneset"))
      bad_selector$feature_sets <- "missing" else
      bad_selector$signal_keys <- "missing"
    missing_method_key <- .trajectory_s2_run_direct(x, operation,
      bad_selector)
    stopifnot(inherits(missing_method_key$value, "error"),
      length(missing_method_key$outputs) == 0L)
    bad_membership <- .trajectory_s2_run_direct(x, operation,
      membership = few)
    stopifnot(inherits(bad_membership$value, "error"),
      length(bad_membership$outputs) == 0L)
    bad_kind <- .trajectory_s2_run_direct(x, operation,
      membership = x$signal)
    stopifnot(inherits(bad_kind$value, "error"),
      length(bad_kind$outputs) == 0L)

    fault <- .trajectory_s2_run_direct(x, operation,
      fault_terminal = "fit_summary")
    stopifnot(inherits(fault$value, "error"),
      grepl("injected terminal fault", conditionMessage(fault$value)),
      !"fit_summary" %in% names(fault$outputs),
      length(fault$outputs) > 0L,
      identical(before, wb_art_list(x$bundle)))
    retry <- .trajectory_s2_run_direct(x, operation)
    stopifnot(!inherits(retry$value, "error"))
    trajectory_s2_compare_outputs(x$bundle,
      paste0("trajectory_", operation, "_dynamics_explicit"),
      list(outputs = retry$outputs), reference)
  }
  TRUE
}
