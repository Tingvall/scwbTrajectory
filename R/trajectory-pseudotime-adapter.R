.trajectory_pseudotime_adaptation <- function(source, lineage_key, ancestry) {
  .trajectory_assert(identical(as.character(source$value_type)[1L],
      "continuous"),
    "trajectory_pseudotime_from_cell_vector requires continuous input.")

  .trajectory_assert(is.character(lineage_key) && length(lineage_key) == 1L &&
      !is.na(lineage_key) && nzchar(trimws(lineage_key)) &&
      nchar(enc2utf8(lineage_key), type = "bytes") <= 1024L,
    "lineage_key must be one non-empty key of at most 1024 bytes.")

  cells <- as.character(source$cell_id)
  .trajectory_assert(length(cells) == length(source$value) && length(cells) > 0L &&
      !anyNA(cells) && all(nzchar(cells)) && !anyDuplicated(cells),
    "Input pseudotime must have one unique non-empty cell ID per value.")
  values <- suppressWarnings(as.numeric(source$value))
  invalid <- is.na(values) & !is.na(source$value)
  .trajectory_assert(!any(invalid),
    "Continuous pseudotime contains non-numeric values.")
  finite <- is.finite(values)
  .trajectory_assert(any(finite),
    "Continuous pseudotime contains no finite values.")

  source_artifact_id <- as.character(ancestry$artifact_id %||% "")[1L]
  source_snapshot_id <- as.character(ancestry$source_snapshot_id %||% "")[1L]
  .trajectory_assert(nzchar(source_artifact_id) && nzchar(source_snapshot_id),
    "Authenticated scalar pseudotime ancestry is incomplete.")

  finite_cells <- cells[finite]
  finite_values <- values[finite]
  counts <- c(input = length(cells), finite = length(finite_cells),
    missing = sum(!finite))
  meta <- list(family_kind = "trajectory_pseudotime",
    value_semantics = "preserved_scalar_pseudotime",
    adaptation = "cell_vector_to_one_key_family_no_transform",
    lineage_key = lineage_key, source_artifact_id = source_artifact_id,
    source_snapshot_id = source_snapshot_id)
  family <- scworkbench::wb_prim_cell_vector_family(
    rep(lineage_key, length(finite_cells)), finite_cells, finite_values,
    value_type = "continuous", meta = meta)
  summary <- scworkbench::wb_prim_obs_table(data.frame(
    lineage_key = lineage_key,
    source_artifact_id = source_artifact_id,
    input_cell_count = as.integer(counts[["input"]]),
    finite_cell_count = as.integer(counts[["finite"]]),
    missing_cell_count = as.integer(counts[["missing"]]),
    stringsAsFactors = FALSE),
    table_type = "trajectory_pseudotime_cell_vector_adaptation_summary",
    meta = meta)

  list(pseudotime = family, adaptation_summary = summary,
    finite_cells = finite_cells)
}

.run_trajectory_pseudotime_from_cell_vector <- function(readers, writers,
    params, context) {
  ancestry_reader <- readers$source_ancestry
  .trajectory_assert(is.function(ancestry_reader),
    paste0("trajectory_pseudotime_from_cell_vector requires an authenticated ",
      "promoted cell-vector artifact."))
  adapted <- .trajectory_pseudotime_adaptation(readers$pseudotime(),
    params$lineage_key, ancestry_reader())

  writers$write_terminal_primitive("pseudotime", adapted$pseudotime)
  writers$write_terminal_primitive("adaptation_summary",
    adapted$adaptation_summary)
  grant <- context$derived_cell_output_axis_grants[[1L]]
  declaration <- scworkbench::wb_derived_cell_output_axis_declaration(
    grant, adapted$finite_cells)
  list(status = "ok",
    derived_cell_output_axis_declarations = list(declaration))
}
