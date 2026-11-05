test_dir <- if (file.exists("helper_s2_fixture.R")) "." else "tests"
source(file.path(test_dir, "helper_s2_fixture.R"), local = .GlobalEnv)

.trajectory_s2_explicit_params <- function(operation) {
  if (identical(operation, "geneset")) return(list(
    namespace = "trajectory", pseudotime_key = "L1", min_set_size = 1L,
    adjust_depth = FALSE,
    n_bins = 5L, min_cells_per_bin = 1L, min_bins_observed = 2L,
    bootstrap_n = 8L, bootstrap_seed = 17L, n_permutations = 8L,
    permutation_seed = 19L, chunk_size = 2L))
  list(pseudotime_key = "L1", n_bins = 5L, min_cells_per_bin = 1L,
    min_bins_observed = 2L, bootstrap_n = 8L, bootstrap_seed = 17L)
}

.trajectory_s2_run_direct <- function(x, operation, params = NULL,
                                      fault_terminal = NULL,
                                      membership = NULL,
                                      ancestry_ids = list()) {
  boundary <- trajectory_s2_mock_boundary(x, fault_terminal)
  if (!is.null(membership)) {
    boundary$readers$membership <- function() membership
  }
  stopifnot(all(names(ancestry_ids) %in% c("pseudotime", "membership")))
  make_ancestry <- function(name, id) {
    force(name); force(id)
    function() {
      key <- paste0(name, "_ancestry")
      boundary$state$calls[[key]] <- as.integer(
        boundary$state$calls[[key]] %||% 0L) + 1L
      list(artifact_id = id)
    }
  }
  for (name in names(ancestry_ids))
    boundary$readers[[paste0(name, "_ancestry")]] <-
      make_ancestry(name, ancestry_ids[[name]])
  runner <- getFromNamespace(paste0(".run_trajectory_", operation,
    "_dynamics_explicit"), "scwbTrajectory")
  if (is.null(params)) params <- .trajectory_s2_explicit_params(operation)
  value <- tryCatch(runner(boundary$readers, boundary$writers, params,
    boundary$context), error = identity)
  list(value = value, boundary = boundary,
    outputs = boundary$state$writes)
}

test_trajectory_s2_explicit_ancestry_metadata_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  for (operation in c("geneset", "signal")) {
    observed <- .trajectory_s2_run_direct(x, operation,
      ancestry_ids = list(pseudotime = "pt_artifact",
        membership = "membership_artifact"))
    stopifnot(!inherits(observed$value, "error"),
      identical(observed$outputs$fit_summary$table$
        pseudotime_artifact_id, "pt_artifact"),
      identical(observed$outputs$fit_summary$table$
        membership_artifact_id, "membership_artifact"),
      identical(observed$boundary$state$calls$pseudotime_ancestry, 1L),
      identical(observed$boundary$state$calls$membership_ancestry, 1L))
    invalid <- .trajectory_s2_run_direct(x, operation,
      ancestry_ids = list(pseudotime = ""))
    stopifnot(inherits(invalid$value, "error"),
      length(invalid$outputs) == 0L)
  }
  TRUE
}

test_trajectory_s2_explicit_runners_oracle_and_seed_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  stopifnot(all(c("scwbtrajectory::trajectory_geneset_dynamics_explicit",
    "scwbtrajectory::trajectory_signal_dynamics_explicit") %in%
    wb_list_ops()))
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
    mixed_selector <- .trajectory_s2_explicit_params(operation)
    if (identical(operation, "geneset"))
      mixed_selector$feature_sets <- c("selected", "missing") else
      mixed_selector$signal_keys <- c("rising", "missing")
    mixed <- .trajectory_s2_run_direct(x, operation, mixed_selector)
    stopifnot(inherits(mixed$value, "error"),
      grepl("requested .* not found", conditionMessage(mixed$value)),
      length(mixed$outputs) == 0L)
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

test_trajectory_s2_authenticated_reader_boundary_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  coverage <- getFromNamespace(".wb_typed_effective_cells", "scworkbench")
  bind <- getFromNamespace(".wb_operation_pack_v2_bind_inputs", "scworkbench")
  check <- function(name, primitive, kind, pattern) {
    rejected <- tryCatch(coverage(x$bundle, "all", list(),
      stats::setNames(list(primitive), name),
      stats::setNames(list(kind), name),
      "trajectory_signal_dynamics", 1L), error = identity)
    stopifnot(inherits(rejected, "error"),
      grepl(pattern, conditionMessage(rejected), fixed = TRUE))
  }
  for (name in c("pseudotime", "signal")) {
    family <- if (identical(name, "pseudotime")) x$explicit else x$signal
    foreign <- family
    foreign$key <- c(foreign$key, foreign$key[[1L]])
    foreign$cell_id <- c(foreign$cell_id, "foreign_cell")
    foreign$value <- c(foreign$value, 0.5)
    check(name, foreign, "cell_vector_family", "unknown cell IDs")
    duplicate <- family
    duplicate$cell_id[[2L]] <- duplicate$cell_id[[1L]]
    check(name, duplicate, "cell_vector_family", "duplicate key/cell rows")
  }
  foreign_membership <- x$membership
  foreign_membership$cell_id[[1L]] <- "foreign_cell"
  check("membership", foreign_membership, "cell_vector", "unknown cell IDs")
  duplicate_membership <- x$membership
  duplicate_membership$cell_id[[2L]] <- duplicate_membership$cell_id[[1L]]
  check("membership", duplicate_membership, "cell_vector", "duplicate cell IDs")

  # A source row outside the requested session subset is valid at this stage.
  subset <- x$cells[1:12]
  accepted <- coverage(x$bundle, "all", list(cell_ids = subset),
    list(pseudotime = x$explicit, signal = x$signal,
      membership = x$membership),
    list(pseudotime = "cell_vector_family",
      signal = "cell_vector_family", membership = "cell_vector"),
    "trajectory_signal_dynamics", 1L)
  stopifnot(identical(accepted$cells, subset))
  bound <- bind(x$bundle, list(pseudotime = x$explicit,
    membership = x$membership), subset)
  stopifnot(identical(bound$pseudotime$cell_id, rep(subset, each = 2L)),
    identical(bound$pseudotime$key, rep(c("L1", "L2"), length(subset))),
    identical(bound$membership$cell_id, subset))

  TRUE
}
