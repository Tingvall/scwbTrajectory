test_dir <- if (file.exists("helper_s2_fixture.R")) "." else "tests"
source(file.path(test_dir, "test_s2_model_execution_contract.R"),
  local = .GlobalEnv)

.trajectory_s2_stream_watch <- new.env(parent = emptyenv())

`[.trajectory_s2_stream_block` <- function(x, ..., drop = FALSE) {
  result <- NextMethod("[", drop = drop)
  if (is.matrix(result))
    class(result) <- c("trajectory_s2_stream_block", "matrix", "array")
  result
}

as.matrix.trajectory_s2_stream_block <- function(x, ...) {
  rows <- nrow(x)
  .trajectory_s2_stream_watch$rows <- c(
    .trajectory_s2_stream_watch$rows, rows)
  if (rows > .trajectory_s2_stream_watch$limit)
    stop("DE materialized more than one feature chunk")
  base::as.matrix(unclass(x), ...)
}

rbind.trajectory_s2_stream_block <- function(..., deparse.level = 1) {
  blocks <- list(...)
  if (sum(vapply(blocks, nrow, integer(1))) >
      .trajectory_s2_stream_watch$limit)
    stop("DE concatenated feature chunks")
  do.call(base::rbind, lapply(blocks, unclass))
}

test_trajectory_s2_de_streaming_block_bound_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  reference <- jsonlite::read_json(file.path(test_dir, "fixtures",
    "trajectory_s2_de_streaming_reference.json"), simplifyVector = FALSE)
  stopifnot(identical(reference$core_commit,
    "ba24635611756240f21f6817044883b832f12701"),
    identical(reference$comparison_policy, "exact_canonical_value"))
  source <- list(engine = "slingshot", axis = x$model$pseudotime_table)
  stream <- getFromNamespace(".trajectory_s2_de_stream", "scwbTrajectory")
  cases <- list(
    all_lineages_all_features = list(selected = x$features,
      params = list()),
    selected_reordered = list(selected = c("g5", "g1", "g3"),
      params = list(features = c("g5", "g1", "g3"), lineage = "L1")))
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    for (size in c(2L, 3L)) {
      requested <- list(); dense_rows <- integer()
      .trajectory_s2_stream_watch$rows <- integer()
      .trajectory_s2_stream_watch$limit <- size
      reader <- function(features, cells) {
        stopifnot(identical(as.character(cells), x$cells))
        requested[[length(requested) + 1L]] <<- as.character(features)
        block <- x$values[features, cells, drop = FALSE]
        class(block) <- c("trajectory_s2_stream_block", "matrix", "array")
        list(values = block)
      }
      observer <- function(n_features, n_cells) {
        stopifnot(identical(as.integer(n_cells),
          as.integer(length(x$cells))), n_features <= size)
        dense_rows <<- c(dense_rows, as.integer(n_features))
      }
      result <- do.call(stream, c(list(reader = reader,
        all_features = x$features,
        selected_features = case$selected, source = source,
        effective_cells = x$cells, matrix_id = "RNA",
        chunk_size = size, observe_block = observer), case$params))
      expected_chunks <- scworkbench::wb_partition_axis(x$features, size)
      stopifnot(identical(requested, expected_chunks),
        length(dense_rows) > 0L,
        identical(dense_rows, .trajectory_s2_stream_watch$rows),
        max(dense_rows) <= size,
        sum(dense_rows) == length(case$selected),
        max(dense_rows) < length(x$features))
      trajectory_s2_compare_outputs(x$bundle, case_name, result, reference)
      old_params <- reference$cases[[case_name]]$params
      meta <- result$outputs$de_table$meta
      for (key in intersect(names(meta), names(old_params)))
        stopifnot(isTRUE(all.equal(meta[[key]], old_params[[key]],
          check.attributes = FALSE)))
    }
  }
  original <- getFromNamespace(".trajectory_chunked_matrix",
    "scwbTrajectory")
  assignInNamespace(".trajectory_chunked_matrix", function(...)
    stop("DE attempted whole-matrix concatenation"),
    ns = "scwbTrajectory")
  on.exit(assignInNamespace(".trajectory_chunked_matrix", original,
    ns = "scwbTrajectory"), add = TRUE)
  for (size in c(2L, 3L)) {
    all <- .trajectory_s2_run_model_direct(x, "de",
      params = list(chunk_size = size))
    stopifnot(!inherits(all$value, "error"),
      identical(all$calls$matrix_chunks, x$features))
    trajectory_s2_compare_outputs(x$bundle,
      "all_lineages_all_features", list(outputs = all$outputs), reference)
    selected <- .trajectory_s2_run_model_direct(x, "de",
      params = list(features = c("g5", "g1", "g3"),
        lineage = "L1", chunk_size = size))
    stopifnot(!inherits(selected$value, "error"),
      identical(selected$calls$matrix_chunks, x$features))
    trajectory_s2_compare_outputs(x$bundle, "selected_reordered",
      list(outputs = selected$outputs), reference)
  }
  TRUE
}
