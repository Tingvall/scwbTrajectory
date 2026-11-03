.trajectory_adapter_state <- function(bundle) {
  list(artifacts = scworkbench::wb_art_list(bundle),
    fingerprint = getFromNamespace(".wb_session_source_identity",
      "scworkbench")(bundle)$scientific_state_fingerprint)
}

.trajectory_adapter_source <- function(bundle, cells, values,
                                       value_type = "continuous") {
  primitive <- scworkbench::wb_prim_cell_vector(
    stats::setNames(values, cells), value_type = value_type)
  run <- scworkbench::wb_new_run("metadata_cell_vector",
    params = list(vector_name = "fallback scalar pseudotime"),
    outputs = list(vector = primitive),
    notes = "Owned scalar-pseudotime adapter fixture.")
  promoted <- scworkbench::wb_promote_analysis(bundle,
    stats::setNames(list(run), run$run_id), head_run_id = run$run_id,
    outputs = "vector", created_by = "pseudotime-adapter-contract")
  list(primitive = primitive,
    artifact_id = as.character(promoted$promoted$artifact_id[[1L]]),
    snapshot_id = as.character(promoted$provenance_snapshot_id)[1L])
}

.trajectory_adapter_expect_error <- function(expression, pattern = NULL) {
  error <- tryCatch(force(expression), error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern))
    stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
  invisible(error)
}

test_trajectory_pseudotime_adapter_lifecycle_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  selected <- x$cells[c(1:7, 10:14, 17:22)]
  values <- c(-4.25, -1, 0, 0.125, pi, 1 / 3, 10.5,
    20, 20.25, 40, 80, 160, 320, 640, 641, 900, 1200, 4096)
  source <- .trajectory_adapter_source(x$bundle, selected, values)
  source_before <- scworkbench::wb_art_load(x$bundle, source$artifact_id,
    expected_kind = "cell_vector")
  source_values <- as.numeric(source_before$value)
  before <- .trajectory_adapter_state(x$bundle)

  manifest <- scworkbench::wb_operation_pack_inspect(
    "scwbtrajectory")$manifest
  operations <- stats::setNames(manifest$operations,
    vapply(manifest$operations, `[[`, character(1), "name"))
  operation <- operations$trajectory_pseudotime_from_cell_vector
  stopifnot(identical(operation$params_schema$required,
      list("lineage_key")),
    identical(vapply(operation$outputs, `[[`, character(1), "name"),
      c("pseudotime", "adaptation_summary")),
    identical(operation$provenance_reads[[1L]]$kind,
      "artifact_ancestry"))

  adapted <- scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_pseudotime_from_cell_vector",
    inputs = list(pseudotime = scworkbench::wb_ref_artifact(
      source$artifact_id, "vector"), lineage_key = "fallback_pt",
      min_cells = 1L))
  family <- adapted$run$outputs$pseudotime
  summary <- adapted$run$outputs$adaptation_summary$table
  stopifnot(
    identical(as.character(family$key), rep("fallback_pt", length(selected))),
    identical(as.character(family$cell_id), selected),
    identical(as.numeric(family$value), source_values),
    identical(as.character(family$meta$source_artifact_id),
      source$artifact_id),
    identical(as.character(family$meta$source_snapshot_id),
      source$snapshot_id),
    identical(as.character(summary$lineage_key), "fallback_pt"),
    identical(as.character(summary$source_artifact_id), source$artifact_id),
    identical(as.integer(summary$input_cell_count), 18L),
    identical(as.integer(summary$finite_cell_count), 18L),
    identical(as.integer(summary$missing_cell_count), 0L),
    identical(before, .trajectory_adapter_state(x$bundle)),
    identical(unclass(source_before), unclass(scworkbench::wb_art_load(
      x$bundle, source$artifact_id, expected_kind = "cell_vector"))))

  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    adapted$session, outputs = c("pseudotime", "adaptation_summary"),
    created_by = "pseudotime-adapter-contract")
  family_id <- as.character(promoted$promoted$artifact_id[
    promoted$promoted$output == "pseudotime"][[1L]])
  transition_id <- as.character(adapted$run$operation_pack_transition$
    transition_id)[1L]
  record_path <- file.path(getFromNamespace(
    ".wb_operation_pack_v2_store", "scworkbench")(x$bundle), transition_id)
  record <- getFromNamespace(".wb_operation_pack_v2_record_validate",
    "scworkbench")(record_path, load_payloads = TRUE)
  ancestry <- jsonlite::read_json(file.path(record_path,
    "provenance_reads", "source_ancestry", "payload.json"),
    simplifyVector = FALSE)
  stopifnot(isTRUE(record$ok),
    identical(ancestry$artifact_id, source$artifact_id),
    identical(ancestry$source_snapshot_id, source$snapshot_id))

  signal <- scworkbench::wb_prim_cell_vector_family(
    rep(c("rising", "falling"), each = length(x$cells)),
    rep(x$cells, 2L),
    c(seq(0, 1, length.out = length(x$cells)),
      seq(1, 0, length.out = length(x$cells))),
    value_type = "continuous")
  consumer <- scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_signal_dynamics_explicit",
    inputs = list(signal = signal,
      pseudotime = scworkbench::wb_ref_artifact(family_id, "pseudotime"),
      pseudotime_key = "fallback_pt", n_bins = 5L,
      min_cells_per_bin = 1L, min_bins_observed = 2L,
      bootstrap_n = 0L, min_cells = 1L))
  fit_summary <- consumer$run$outputs$fit_summary$table
  stopifnot(identical(as.character(fit_summary$pseudotime_key),
      "fallback_pt"),
    identical(as.character(fit_summary$pseudotime_artifact_id), family_id),
    identical(as.integer(fit_summary$n_pseudotime_cells), 18L),
    identical(unclass(source_before), unclass(scworkbench::wb_art_load(
      x$bundle, source$artifact_id, expected_kind = "cell_vector"))))
  invisible(scworkbench::wb_session_discard(consumer$session))
  TRUE
}

test_trajectory_pseudotime_adapter_missing_and_rejection_contract <-
    function() {
  adapt <- getFromNamespace(".trajectory_pseudotime_adaptation",
    "scwbTrajectory")
  ancestry <- list(artifact_id = "artifact_source",
    source_snapshot_id = "snapshot_source")
  source <- scworkbench::wb_prim_cell_vector(stats::setNames(
    c(-3.25, NA_real_, 9.75, NA_real_), c("c4", "c2", "c1", "c3")),
    value_type = "continuous")
  source_before <- serialize(source, NULL, version = 3L)
  adapted <- adapt(source, "fallback_pt", ancestry)
  stopifnot(
    identical(adapted$pseudotime$cell_id, c("c4", "c1")),
    identical(adapted$pseudotime$value, c(-3.25, 9.75)),
    identical(unique(adapted$pseudotime$key), "fallback_pt"),
    identical(adapted$adaptation_summary$table$input_cell_count, 4L),
    identical(adapted$adaptation_summary$table$finite_cell_count, 2L),
    identical(adapted$adaptation_summary$table$missing_cell_count, 2L),
    identical(serialize(source, NULL, version = 3L), source_before))

  categorical <- scworkbench::wb_prim_cell_vector(
    c(c1 = "early", c2 = "late"), value_type = "categorical")
  duplicate <- structure(list(kind = "cell_vector",
    cell_id = c("c1", "c1"), value = c(1, 2),
    value_type = "continuous"), class = "wb_primitive")
  no_finite <- scworkbench::wb_prim_cell_vector(
    c(c1 = NA_real_, c2 = NA_real_), value_type = "continuous")
  .trajectory_adapter_expect_error(adapt(categorical, "fallback_pt",
    ancestry), "requires continuous input")
  .trajectory_adapter_expect_error(adapt(source, " ", ancestry),
    "lineage_key must be one non-empty key")
  .trajectory_adapter_expect_error(adapt(duplicate, "fallback_pt",
    ancestry), "unique non-empty cell ID")
  .trajectory_adapter_expect_error(adapt(no_finite, "fallback_pt",
    ancestry), "contains no finite values")

  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  categorical_source <- .trajectory_adapter_source(x$bundle, x$cells,
    rep(c("early", "late"), length.out = length(x$cells)),
    value_type = "categorical")
  valid_source <- .trajectory_adapter_source(x$bundle, x$cells[1:12],
    seq(0, 1, length.out = 12L))
  before <- .trajectory_adapter_state(x$bundle)
  run <- function(input, key = "fallback_pt")
    scworkbench::wb_session_run_op(scworkbench::wb_session_new(), x$bundle,
      "scwbtrajectory::trajectory_pseudotime_from_cell_vector",
      inputs = list(pseudotime = input, lineage_key = key,
        min_cells = 1L))
  .trajectory_adapter_expect_error(run(scworkbench::wb_ref_artifact(
    categorical_source$artifact_id, "vector")))
  .trajectory_adapter_expect_error(run(x$pseudotime, "fallback_pt"),
    "requires an authenticated promoted cell-vector artifact")
  .trajectory_adapter_expect_error(run(scworkbench::wb_ref_artifact(
    valid_source$artifact_id, "vector"), " "),
    "lineage_key must be one non-empty key")
  .trajectory_adapter_expect_error(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_pseudotime_from_cell_vector",
    inputs = list(pseudotime = scworkbench::wb_ref_artifact(
      valid_source$artifact_id, "vector"), min_cells = 1L)))
  .trajectory_adapter_expect_error(run(duplicate))
  .trajectory_adapter_expect_error(run(no_finite))
  stopifnot(identical(before, .trajectory_adapter_state(x$bundle)))
  TRUE
}
