test_dir <- if (file.exists("helper_s2_fixture.R")) "." else "tests"
source(file.path(test_dir, "helper_s2_fixture.R"), local = .GlobalEnv)
source(file.path(test_dir, "test_s2_explicit_execution_contract.R"),
  local = .GlobalEnv)
source(file.path(test_dir, "test_s2_model_execution_contract.R"),
  local = .GlobalEnv)

test_trajectory_s2_frozen_method_outputs_contract <- function() {
  if (!requireNamespace("slingshot", quietly = TRUE))
    return(structure(TRUE, skip = "slingshot unavailable"))
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  reference <- trajectory_s2_reference()
  stopifnot(identical(reference$core_commit,
    "ba24635611756240f21f6817044883b832f12701"),
    identical(reference$oracle, "old_core_synthetic_pseudotime_table"))
  cases <- list(
    trajectory_de = function() .trajectory_s2_run_model_direct(x, "de"),
    trajectory_smooth = function()
      .trajectory_s2_run_model_direct(x, "smooth"),
    trajectory_signature_smooth = function()
      .trajectory_s2_run_model_direct(x, "signature_smooth"),
    trajectory_geneset_dynamics_model = function()
      .trajectory_s2_run_model_direct(x, "geneset_dynamics"),
    trajectory_geneset_dynamics_explicit = function()
      .trajectory_s2_run_direct(x, "geneset"),
    trajectory_vector_smooth = function()
      .trajectory_s2_run_model_direct(x, "vector_smooth"),
    trajectory_signal_dynamics_model = function()
      .trajectory_s2_run_model_direct(x, "signal_dynamics"),
    trajectory_signal_dynamics_explicit = function()
      .trajectory_s2_run_direct(x, "signal"))
  stopifnot(identical(sort(names(cases)), sort(names(reference$cases))))
  results <- list()
  for (name in names(cases)) {
    observed <- cases[[name]]()
    stopifnot(!inherits(observed$value, "error"))
    results[[name]] <- observed
    trajectory_s2_compare_outputs(x$bundle, name, observed, reference)
    if (grepl("_dynamics_(model|explicit)$", name)) {
      expected_source <- if (grepl("_model$", name)) "model_ref" else
        "cell_vector_family"
      stopifnot(identical(observed$outputs$fit_summary$table$
        pseudotime_source_kind, expected_source),
        identical(observed$outputs$fit_summary$table$pseudotime_key, "L1"))
    }
  }
  # Equal pseudotime values give equal method tables, while the fit summary
  # still records which authenticated source was selected.
  for (operation in c("trajectory_geneset_dynamics",
      "trajectory_signal_dynamics")) {
    model <- results[[paste0(operation, "_model")]]$outputs
    explicit <- results[[paste0(operation, "_explicit")]]$outputs
    for (terminal in setdiff(names(model), "fit_summary"))
      stopifnot(isTRUE(all.equal(model[[terminal]]$table,
        explicit[[terminal]]$table, check.attributes = FALSE,
        tolerance = 0)))
  }

  # A comparator must detect both a missing sibling and changed table data.
  original <- cases$trajectory_signal_dynamics_explicit()
  missing <- original
  missing$outputs$top_signals <- NULL
  stopifnot(inherits(try(trajectory_s2_compare_outputs(x$bundle,
    "trajectory_signal_dynamics_explicit", missing, reference),
    silent = TRUE), "try-error"))
  changed <- original
  changed$outputs$signal_table$table$value[1L] <-
    changed$outputs$signal_table$table$value[1L] + 0.25
  stopifnot(inherits(try(trajectory_s2_compare_outputs(x$bundle,
    "trajectory_signal_dynamics_explicit", changed, reference),
    silent = TRUE), "try-error"))
  TRUE
}

test_trajectory_s2_explicit_pseudotime_selection_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  resolve <- getFromNamespace(".trajectory_s2_explicit_lineages",
    "scwbTrajectory")
  selected <- resolve(x$explicit, pseudotime_key = "L1",
    membership = x$membership, cells = x$cells[1:12])
  stopifnot(identical(names(selected), "L1"),
    identical(selected$L1$cell_id, x$cells[1:12]),
    isTRUE(all.equal(selected$L1$pseudotime,
      seq(0, 1, length.out = 24L)[1:12], check.attributes = FALSE)))
  stopifnot(inherits(try(resolve(x$explicit), silent = TRUE), "try-error"),
    inherits(try(resolve(x$explicit, pseudotime_key = "missing"),
      silent = TRUE), "try-error"),
    inherits(try(resolve(x$explicit, pseudotime_key = "L1",
      cells = "foreign"), silent = TRUE), "try-error"))
  TRUE
}

test_trajectory_s2_summary_equations_contract <- function() {
  bins <- getFromNamespace(".trajectory_s2_bins", "scwbTrajectory")
  aggregate <- getFromNamespace(".trajectory_s2_aggregate_one",
    "scwbTrajectory")
  reference <- jsonlite::read_json(file.path(test_dir, "fixtures",
    "trajectory_s2_parameter_reference.json"), simplifyVector = TRUE)
  stopifnot(identical(reference$core_commit,
    "ba24635611756240f21f6817044883b832f12701"))
  pt <- seq(0, 1, length.out = 24L)
  values <- c(seq(1, 12), rep(NA_real_, 2L), seq(13, 22))
  observed_bins <- bins(pt, 5L)
  stopifnot(identical(observed_bins$bins,
      as.integer(reference$summary_equations$bins$bins)),
    identical(observed_bins$n_bins,
      as.integer(reference$summary_equations$bins$n_bins)))
  for (summary_fun in c("mean", "median")) {
    observed <- aggregate("g1", "L1", pt, values,
      observed_bins$bins, 2L, summary_fun = summary_fun)
    expected <- reference$summary_equations[[summary_fun]]
    stopifnot(isTRUE(all.equal(observed, expected,
      check.attributes = FALSE, tolerance = 1e-14)))
  }
  stopifnot(is.null(aggregate("g1", "L1", pt, values,
    observed_bins$bins, 24L)),
    inherits(try(aggregate("g1", "L1", pt, values,
      observed_bins$bins, 1L, summary_fun = "unsupported"),
      silent = TRUE), "try-error"))
  TRUE
}
