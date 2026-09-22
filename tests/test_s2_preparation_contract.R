source(file.path(test_dir, "helper_s2_fixture.R"), local = .GlobalEnv)

test_trajectory_s2_frozen_method_outputs_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  reference <- trajectory_s2_reference()
  stopifnot(identical(reference$core_commit,
    "ba24635611756240f21f6817044883b832f12701"),
    identical(reference$oracle, "old_core_synthetic_pseudotime_table"))
  cases <- trajectory_s2_cases(x)
  stopifnot(identical(sort(names(cases)), sort(names(reference$cases))))
  results <- list()
  for (name in names(cases)) {
    observed <- cases[[name]]()
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
  old <- getFromNamespace(".wb_trajectory_dynamics_resolve_lineages",
    "scworkbench")(bundle = x$bundle, pseudotime = x$explicit,
      pseudotime_key = "L1", membership = x$membership,
      cells = x$cells[1:12])$lineages
  stopifnot(isTRUE(all.equal(selected$L1, old$L1,
    check.attributes = FALSE)))
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
  old_bins <- getFromNamespace(".wb_trajectory_smooth_bins", "scworkbench")
  old_aggregate <- getFromNamespace(".wb_trajectory_smooth_aggregate_one",
    "scworkbench")
  pt <- seq(0, 1, length.out = 24L)
  values <- c(seq(1, 12), rep(NA_real_, 2L), seq(13, 22))
  observed_bins <- bins(pt, 5L)
  stopifnot(identical(observed_bins, old_bins(pt, 5L)))
  for (summary_fun in c("mean", "median")) {
    observed <- aggregate("g1", "L1", pt, values,
      observed_bins$bins, 2L, summary_fun = summary_fun)
    expected <- old_aggregate("g1", "L1", pt, values,
      observed_bins$bins, 2L, summary_fun = summary_fun)
    stopifnot(isTRUE(all.equal(observed, expected, tolerance = 0)))
  }
  stopifnot(is.null(aggregate("g1", "L1", pt, values,
    observed_bins$bins, 24L)),
    inherits(try(aggregate("g1", "L1", pt, values,
      observed_bins$bins, 1L, summary_fun = "unsupported"),
      silent = TRUE), "try-error"))
  TRUE
}
