.trajectory_occupancy_inputs <- function(cells) {
  stopifnot(length(cells) >= 8L)
  cells <- cells[seq_len(8L)]
  pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("mapped_pt", 8L), cells,
    c(0, .25, .75, 1, .1, .4, .8, NA_real_),
    value_type = "continuous")
  group_id <- scworkbench::wb_prim_cell_vector(stats::setNames(
    rep(c("A", "B"), each = 4L), cells), value_type = "categorical")
  sample_id <- scworkbench::wb_prim_cell_vector(stats::setNames(
    rep(c("s1", "s2"), times = 2L, each = 2L), cells),
    value_type = "categorical")
  accepted <- c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
  reason <- c("accepted", "accepted", "annotation_not_accepted", "accepted",
    "accepted", "annotation_confidence_too_low", "accepted", "accepted")
  mapping_acceptance <- scworkbench::wb_prim_obs_table(data.frame(
    query_cell_id = cells, reference_key = "trajectory",
    output_key = "mapped_pt", accepted = accepted,
    acceptance_reason = reason,
    raw_estimate = pseudotime$value,
    accepted_value = ifelse(accepted, pseudotime$value, NA_real_),
    stringsAsFactors = FALSE),
    table_type = "reference_map_continuous_acceptance",
    entity_kind = "query_cell_vector")
  list(cells = cells, pseudotime = pseudotime, group_id = group_id,
    sample_id = sample_id, mapping_acceptance = mapping_acceptance)
}

test_trajectory_occupancy_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  values <- .trajectory_occupancy_inputs(x$cells)
  inputs <- list(pseudotime = values$pseudotime,
    group_id = values$group_id, sample_id = values$sample_id,
    mapping_acceptance = values$mapping_acceptance,
    pseudotime_key = "mapped_pt", acceptance_output_key = "mapped_pt",
    n_bins = 2L, analysis_mode = "descriptive", min_cells = 1L)
  run <- function(extra = list()) scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_occupancy",
    inputs = utils::modifyList(inputs, extra))
  first <- run(); second <- run()
  stopifnot(identical(first$run$outputs, second$run$outputs),
    setequal(names(first$run$outputs),
      c("bin_occupancy", "acceptance_accounting", "run_summary")))

  occupancy <- first$run$outputs$bin_occupancy$table
  aggregate <- occupancy[occupancy$stratum == "group", , drop = FALSE]
  group_a <- aggregate[aggregate$group_id == "A", , drop = FALSE]
  group_b <- aggregate[aggregate$group_id == "B", , drop = FALSE]
  stopifnot(identical(group_a$n_cells, c(2L, 1L)),
    identical(group_b$n_cells, c(1L, 1L)),
    isTRUE(all.equal(group_a$fraction_within_stratum, c(2 / 3, 1 / 3))),
    isTRUE(all.equal(group_b$fraction_within_stratum, c(.5, .5))),
    nrow(occupancy) == 12L,
    all(occupancy$lower_inclusive),
    identical(occupancy$upper_inclusive[seq_len(2L)], c(FALSE, TRUE)))

  accounting <- first$run$outputs$acceptance_accounting$table
  aggregate_accounting <- accounting[accounting$stratum == "group", ,
    drop = FALSE]
  stopifnot(sum(aggregate_accounting$n_cells) == 8L,
    sum(aggregate_accounting$n_cells[
      aggregate_accounting$accounting_status == "accepted_binned"]) == 5L,
    sum(aggregate_accounting$n_cells[
      aggregate_accounting$accounting_status ==
        "accepted_pseudotime_missing"]) == 1L,
    sum(aggregate_accounting$n_cells[
      aggregate_accounting$accounting_status ==
        "annotation_not_accepted"]) == 1L,
    sum(aggregate_accounting$n_cells[
      aggregate_accounting$accounting_status ==
        "annotation_confidence_too_low"]) == 1L)

  summary <- first$run$outputs$run_summary$table
  stopifnot(identical(summary$n_cells_input, 8L),
    identical(summary$n_cells_accepted, 6L),
    identical(summary$n_cells_withheld, 2L),
    identical(summary$n_cells_accepted_finite, 5L),
    identical(summary$n_cells_accepted_pseudotime_missing, 1L),
    identical(summary$n_groups, 2L), identical(summary$n_samples, 2L),
    isTRUE(summary$accounting_complete),
    identical(first$run$params$n_bins, 2L),
    identical(first$run$params$analysis_mode, "descriptive"))

  manifest <- scworkbench::wb_operation_pack_inspect(
    "scwbtrajectory")$manifest
  operation <- Filter(function(value) identical(value$name,
    "trajectory_occupancy"), manifest$operations)[[1L]]
  acceptance_input <- Filter(function(value) identical(value$name,
    "mapping_acceptance"), operation$inputs)[[1L]]
  stopifnot(identical(acceptance_input$effective_cell_coverage,
      "reference_only"),
    identical(operation$params_schema$properties$analysis_mode$enum,
      list("descriptive")),
    all(vapply(operation$promotion_closure, function(value)
      setequal(unlist(value, use.names = FALSE), names(first$run$outputs)),
      logical(1))))

  inferential <- try(run(list(analysis_mode = "inferential")), silent = TRUE)
  incomplete <- values$mapping_acceptance
  incomplete$table <- incomplete$table[-1L, , drop = FALSE]
  incomplete_run <- try(run(list(mapping_acceptance = incomplete)), silent = TRUE)
  stopifnot(inherits(inferential, "try-error"),
    inherits(incomplete_run, "try-error"))
  TRUE
}
