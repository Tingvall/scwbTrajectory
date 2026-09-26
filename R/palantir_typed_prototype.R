# Typed Phase 5 Palantir computation, initially proven as an unregistered
# prototype and now wired to the isolated candidate pack manifest. Historical
# core operations remain authoritative until a separate ownership order.
.trajectory_palantir_typed_schema <-
  "pybench.prototype.palantir_typed.v1"

.trajectory_palantir_typed_numeric <- function(value, n, label) {
  result <- as.numeric(unlist(value, use.names = FALSE))
  .trajectory_assert(length(result) == n && all(is.finite(result)),
    paste0("Private Palantir response has invalid ", label, "."))
  result
}

.trajectory_palantir_typed_rows <- function(value, fields, label) {
  .trajectory_assert(is.list(value) && length(value) > 0L &&
      all(vapply(value, function(row) is.list(row) &&
        setequal(names(row), fields), logical(1))),
    paste0("Private Palantir response has invalid ", label, " rows."))
  value
}

.trajectory_palantir_typed_matrix_rows <- function(values) {
  lapply(seq_len(nrow(values)), function(i)
    unname(as.list(as.numeric(values[i, , drop = TRUE]))))
}

.trajectory_palantir_typed_effective <- function(response, requested,
    mode) {
  actual <- response$effective_params
  if (identical(mode, "fit")) {
    fields <- c("knn", "n_diffusion_components", "n_waypoints", "seed",
      "use_early_cell_as_start", "dims")
    valid <- is.list(actual) && setequal(names(actual), fields) &&
      identical(as.integer(actual$knn), requested$knn) &&
      identical(as.integer(actual$n_diffusion_components),
        requested$n_diffusion_components) &&
      identical(as.integer(actual$n_waypoints), requested$n_waypoints) &&
      identical(as.integer(actual$seed), requested$seed) &&
      identical(isTRUE(actual$use_early_cell_as_start),
        requested$use_early_cell_as_start) &&
      identical(as.integer(unlist(actual$dims, use.names = FALSE)),
        as.integer(unlist(requested$dims, use.names = FALSE)))
  } else {
    fields <- c("prob_min", "smoothness", "sample_points",
      "min_cells_per_branch", "mask_q", "mask_eps")
    valid <- is.list(actual) && setequal(names(actual), fields) &&
      identical(as.numeric(actual$prob_min), requested$prob_min) &&
      identical(as.numeric(actual$smoothness), requested$smoothness) &&
      identical(as.integer(actual$sample_points),
        max(10L, requested$sample_points)) &&
      identical(as.integer(actual$min_cells_per_branch),
        max(3L, requested$min_cells_per_branch)) &&
      identical(as.numeric(actual$mask_q), 0.01) &&
      identical(as.numeric(actual$mask_eps), 0.01)
  }
  .trajectory_assert(isTRUE(valid),
    "Private Palantir effective parameters differ from the typed request.")
  invisible(TRUE)
}

.trajectory_palantir_typed_call <- function(request, model_bytes = NULL) {
  .trajectory_assert(requireNamespace("jsonlite", quietly = TRUE),
    "The unregistered Palantir prototype requires jsonlite.")
  stage <- tempfile("scwb-palantir-private-")
  dir.create(stage, mode = "0700")
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  if (!is.null(model_bytes)) {
    .trajectory_assert(is.raw(model_bytes) && length(model_bytes) > 0L,
      "Verified G1 model bytes are empty.")
    writeBin(model_bytes, file.path(stage, "model.pkl"))
    request$model <- list(path = "model.pkl", sha256 =
      digest::digest(model_bytes, algo = "sha256", serialize = FALSE))
  }
  request$schema <- .trajectory_palantir_typed_schema
  request_path <- file.path(stage, "request.json")
  response_path <- file.path(stage, "response.json")
  jsonlite::write_json(request, request_path, auto_unbox = TRUE,
    null = "null", digits = NA)
  python <- as.character(getOption(
    "scwbTrajectory.palantir_typed_python",
    Sys.getenv("SCWORKBENCH_PYBENCH_PYTHON", "python3")))[1L]
  .trajectory_assert(nzchar(python),
    "Private Palantir Python interpreter is not configured.")
  output <- suppressWarnings(system2(python, c("-m",
    "pybench.prototype.palantir_typed", "--request",
    shQuote(request_path), "--response", shQuote(response_path)),
    stdout = TRUE, stderr = TRUE))
  status <- as.integer(attr(output, "status") %||% 0L)
  .trajectory_assert(identical(status, 0L) && file.exists(response_path),
    paste0("Private Palantir computation failed before terminal writes: ",
      paste(output, collapse = "\n")))
  response <- jsonlite::read_json(response_path, simplifyVector = FALSE)
  .trajectory_assert(identical(as.character(response$schema),
      .trajectory_palantir_typed_schema) &&
      identical(as.character(response$mode), as.character(request$mode)),
    "Private Palantir response schema or mode changed.")
  if (identical(request$mode, "fit")) {
    descriptor <- response$model
    .trajectory_assert(is.list(descriptor) &&
        identical(as.character(descriptor$path), "model.pkl"),
      "Private Palantir model descriptor is invalid.")
    path <- file.path(stage, "model.pkl")
    .trajectory_assert(file.exists(path) && !dir.exists(path),
      "Private Palantir model bytes are absent.")
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    .trajectory_assert(identical(as.integer(descriptor$size_bytes),
        as.integer(length(bytes))) &&
        identical(as.character(descriptor$sha256),
          digest::digest(bytes, algo = "sha256", serialize = FALSE)),
      "Private Palantir model descriptor disagrees with model bytes.")
    response$model_bytes <- bytes
  }
  response
}

.trajectory_palantir_typed_fit <- function(readers, writers, params,
    context) {
  space <- .trajectory_space_matrix(readers$fit_space())
  cells <- as.character(context$effective_cells)
  index <- match(cells, rownames(space))
  .trajectory_assert(length(cells) >= 10L && !anyDuplicated(cells) &&
      !anyNA(index),
    "Palantir fitted cells must be an ordered subset of the typed space.")
  dims <- as.integer(params$dims %||% seq_len(ncol(space)))
  .trajectory_assert(length(dims) >= 2L && !anyNA(dims) &&
      !anyDuplicated(dims) && all(dims >= 1L & dims <= ncol(space)),
    "Palantir dimensions are invalid for the typed space.")
  root <- as.character(params$root_cell)[1L]
  terminals <- as.character(params$terminal_cells %||% character())
  .trajectory_assert(length(root) == 1L && !is.na(root) &&
      root %in% cells && length(terminals) > 0L &&
      !anyDuplicated(terminals) && all(terminals %in% cells),
    "Palantir root or terminal cells are invalid.")
  effective <- list(dims = unname(as.list(dims)),
    knn = as.integer(params$knn %||% 30L),
    n_diffusion_components = as.integer(
      params$n_diffusion_components %||% 10L),
    n_waypoints = as.integer(params$n_waypoints %||% 300L),
    use_early_cell_as_start = isTRUE(
      params$use_early_cell_as_start %||% TRUE),
    seed = as.integer(params$seed %||% 20L))
  request <- list(mode = "fit", cell_ids = unname(as.list(cells)),
    coordinates = .trajectory_palantir_typed_matrix_rows(
      space[index, , drop = FALSE]), root_cell = root,
    terminal_cells = unname(as.list(terminals)), params = effective)
  response <- .trajectory_palantir_typed_call(request)
  .trajectory_palantir_typed_effective(response, effective, "fit")
  .trajectory_assert(identical(as.character(unlist(response$cell_ids,
      use.names = FALSE)), cells),
    "Private Palantir fitted cell order changed.")
  pseudotime <- .trajectory_palantir_typed_numeric(response$pseudotime,
    length(cells), "pseudotime")
  entropy <- .trajectory_palantir_typed_numeric(response$entropy,
    length(cells), "entropy")
  relation_rows <- .trajectory_palantir_typed_rows(
    response$branch_membership, c("key", "cell_id", "weight"),
    "branch membership")
  keys <- vapply(relation_rows, function(row) as.character(row$key),
    character(1))
  relation_cells <- vapply(relation_rows,
    function(row) as.character(row$cell_id), character(1))
  weight <- .trajectory_palantir_typed_numeric(lapply(relation_rows, `[[`,
    "weight"), length(relation_rows), "branch membership weight")
  .trajectory_assert(all(relation_cells %in% cells) &&
      all(weight >= 0 & weight <= 1),
    "Private Palantir relation contains foreign cells or invalid weights.")
  summary_rows <- .trajectory_palantir_typed_rows(response$branch_summary,
    c("relation_key", "branch", "n_cells_prob_gt_0", "mean_prob",
      "max_prob", "top_cell_id"), "branch summary")
  summary <- data.frame(
    relation_key = vapply(summary_rows, function(row)
      as.character(row$relation_key), character(1)),
    branch = vapply(summary_rows, function(row)
      as.character(row$branch), character(1)),
    n_cells_prob_gt_0 = vapply(summary_rows, function(row)
      as.integer(row$n_cells_prob_gt_0), integer(1)),
    mean_prob = .trajectory_palantir_typed_numeric(lapply(summary_rows,
      `[[`, "mean_prob"), length(summary_rows), "branch mean"),
    max_prob = .trajectory_palantir_typed_numeric(lapply(summary_rows,
      `[[`, "max_prob"), length(summary_rows), "branch maximum"),
    top_cell_id = vapply(summary_rows, function(row)
      as.character(row$top_cell_id), character(1)),
    stringsAsFactors = FALSE)
  .trajectory_assert(identical(as.character(summary$relation_key),
      unique(keys)) && all(summary$top_cell_id %in% cells) &&
      all(!is.na(summary$n_cells_prob_gt_0)) &&
      all(summary$n_cells_prob_gt_0 >= 0L &
        summary$n_cells_prob_gt_0 <= length(cells)),
    "Private Palantir branch summary contradicts relation keys or cells.")
  fitted_space_id <- as.character(context$input_ids$fit_space)[1L]
  .trajectory_assert(length(fitted_space_id) == 1L &&
      nzchar(fitted_space_id), "Authenticated fit-space ID is unavailable.")
  relation <- scworkbench::wb_prim_relation_family(keys, relation_cells,
    weight = weight)
  pseudotime <- scworkbench::wb_prim_cell_vector(
    stats::setNames(pseudotime, cells), value_type = "continuous")
  entropy <- scworkbench::wb_prim_cell_vector(
    stats::setNames(entropy, cells), value_type = "continuous")
  summary <- scworkbench::wb_prim_obs_table(summary,
    table_type = "palantir_branch_summary", entity_kind = "branch")
  pairs <- lapply(seq_along(keys), function(i)
    list(key = keys[[i]], cell_id = relation_cells[[i]]))
  writers$declare_derived_fitted_axis(cells)
  writers$declare_relation_family_rows("branch_membership", unique(keys),
    list(pairs), c(weight = TRUE, value = FALSE, distance = FALSE))
  writers$write_terminal_primitive("pseudotime", pseudotime)
  writers$write_terminal_primitive("branch_membership", relation)
  writers$write_terminal_primitive("branch_summary", summary)
  writers$write_terminal_primitive("entropy", entropy)
  tree <- writers$open_model_tree("model")
  writeBin(response$model_bytes, tree$path(tree$primary_path))
  model <- scworkbench::wb_prim_model_ref(fitted_space_id,
    fitted_space_id, engine = "palantir")
  model$model_type <- "palantir_fit"
  model$cells <- cells
  model$fit_space_id <- fitted_space_id
  tree$finish(model)
  invisible(list(status = "ok"))
}

.trajectory_palantir_typed_project <- function(readers, writers, params,
    context) {
  # A matching digest proves byte integrity, not that a pickle is safe to
  # execute. The operator must explicitly attest trust in the model source.
  .trajectory_assert(identical(as.character(
      params$pickle_execution_trust %||% "")[1L],
      "operator_verified_source_v1"),
    "Palantir pickle execution requires an explicitly verified model source.")
  verified <- readers$fitted_model()
  .trajectory_assert(identical(as.character(verified$primary_path),
      "model.pkl"), "Authenticated Palantir model primary path changed.")
  .trajectory_assert(is.character(verified$source_space_id) &&
      length(verified$source_space_id) == 1L &&
      nzchar(verified$source_space_id) &&
      is.matrix(verified$fitted_source_space) &&
      identical(as.character(rownames(verified$fitted_source_space)),
        as.character(context$effective_cells)),
    "Authenticated Palantir model source linkage or fitted axis changed.")
  entries <- verified$entries
  .trajectory_assert(is.list(entries) && length(entries) >= 1L,
    "Authenticated Palantir model tree is empty.")
  primary <- Filter(function(entry) is.list(entry) &&
      identical(as.character(entry$path), "model.pkl"), entries)
  .trajectory_assert(length(primary) == 1L,
    "Authenticated Palantir model descriptor is ambiguous.")
  model_bytes <- verified$read_bytes(verified$primary_path)
  .trajectory_assert(is.raw(model_bytes) && length(model_bytes) > 0L &&
      identical(as.integer(primary[[1L]]$size),
        as.integer(length(model_bytes))) &&
      identical(as.character(primary[[1L]]$sha256),
        digest::digest(model_bytes, algo = "sha256", serialize = FALSE)),
    "Authenticated Palantir model checksum changed before private execution.")
  display <- .trajectory_space_matrix(readers$plot_space())
  .trajectory_assert(ncol(display) == 2L,
    "Palantir display space is not two-dimensional.")
  cells <- as.character(context$effective_cells)
  index <- match(cells, rownames(display))
  .trajectory_assert(length(cells) > 0L && !anyDuplicated(cells) &&
      !anyNA(index),
    "Palantir display space lacks the exact ordered fitted cell axis.")
  effective <- list(prob_min = as.numeric(params$prob_min %||% 0.05),
    smoothness = as.numeric(params$smoothness %||% 1),
    min_cells_per_branch = as.integer(
      params$min_cells_per_branch %||% 50L),
    sample_points = as.integer(params$approx_points %||% 200L))
  response <- .trajectory_palantir_typed_call(list(mode = "project",
    cell_ids = unname(as.list(cells)),
    coordinates = .trajectory_palantir_typed_matrix_rows(
      display[index, , drop = FALSE]), params = effective), model_bytes)
  .trajectory_palantir_typed_effective(response, effective, "project")
  rows <- .trajectory_palantir_typed_rows(response$curves,
    c("key", "point", "x", "y"), "projected curve")
  keys <- vapply(rows, function(row) as.character(row$key), character(1))
  points <- vapply(rows, function(row) as.integer(row$point), integer(1))
  x <- .trajectory_palantir_typed_numeric(lapply(rows, `[[`, "x"),
    length(rows), "curve x")
  y <- .trajectory_palantir_typed_numeric(lapply(rows, `[[`, "y"),
    length(rows), "curve y")
  counts <- response$branch_counts
  .trajectory_assert(is.list(counts) &&
      identical(names(counts), unique(keys)) &&
      all(vapply(counts, function(value) length(value) == 1L &&
        is.finite(as.numeric(value)) && as.numeric(value) >= 3,
        logical(1))) &&
      all(vapply(unique(keys), function(key) identical(
        as.integer(points[keys == key]),
        seq_len(sum(keys == key))), logical(1))),
    "Private Palantir curve points or branch counts are inconsistent.")
  curves <- scworkbench::wb_prim_path_2d_family(keys, points, x, y,
    meta = list(plot_space_id = as.character(context$input_ids$plot_space)[1L],
      fit_space_id = as.character(verified$source_space_id)[1L],
      curve_method = "palantir_mellon_interpolation"))
  writers$write_path_2d_family("curves", curves)
  invisible(list(status = "ok"))
}

# These two production entry points remain pack-owned; core retains its
# historical names until a separately authorized ownership migration.
.run_palantir_fit <- .trajectory_palantir_typed_fit
.run_palantir_project <- .trajectory_palantir_typed_project
