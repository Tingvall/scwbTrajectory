# Single change-point detection along a pseudotime axis. The numeric core lives
# in inst/scworkbench/changepoint_core.py and runs in private Python compute;
# this file owns axis ordering, region semantics and every terminal write.

.trajectory_cp_schema <- "scwbtrajectory.changepoint_core.v1"
.trajectory_cp_prefix <- "trajectory_single_changepoint_regions"
.trajectory_cp_codes <- c("2" = "RUNTIME_IMPORT_UNAVAILABLE",
  "3" = "RUNTIME_API_UNAVAILABLE", "4" = "REQUEST_INVALID",
  "5" = "COMPUTE_FAILED")

.trajectory_cp_stop <- function(message)
  stop(paste0(.trajectory_cp_prefix, ": ", message), call. = FALSE)

.trajectory_cp_assert <- function(ok, message)
  .trajectory_assert(isTRUE(ok), paste0(.trajectory_cp_prefix, ": ", message))

.trajectory_cp_python <- function() {
  value <- as.character(getOption("scwbTrajectory.changepoint_python",
    Sys.getenv("SCWORKBENCH_CHANGEPPOINT_PYTHON", "python3")))[1L]
  .trajectory_cp_assert(length(value) == 1L && !is.na(value) && nzchar(value),
    "the private-compute interpreter is not configured.")
  value
}

.trajectory_cp_script <- function() {
  path <- system.file("scworkbench", "changepoint_core.py",
    package = "scwbTrajectory")
  .trajectory_cp_assert(length(path) == 1L && nzchar(path) &&
      file.exists(path) && !dir.exists(path),
    "the private numeric core is not installed with this pack.")
  normalizePath(path, mustWork = TRUE)
}

.trajectory_cp_call <- function(phase, request, python, script) {
  .trajectory_cp_assert(requireNamespace("jsonlite", quietly = TRUE),
    "the private-compute bridge requires jsonlite.")
  stage <- tempfile("scwb-changepoint-")
  dir.create(stage, mode = "0700")
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  request_path <- file.path(stage, "request.json")
  response_path <- file.path(stage, "response.json")
  jsonlite::write_json(request, request_path, auto_unbox = TRUE,
    na = "null", null = "null", digits = NA)
  # A missing interpreter surfaces either as status 127 from the shell or as an
  # R-level "error in running command"; both mean the same thing here.
  launched <- tryCatch(list(output = suppressWarnings(system2(python,
      c(shQuote(script), "--phase", phase, "--request", shQuote(request_path),
        "--response", shQuote(response_path)), stdout = TRUE, stderr = TRUE)),
      launch_error = NA_character_),
    error = function(error) list(output = character(0),
      launch_error = conditionMessage(error)))
  output <- launched$output
  if (!is.na(launched$launch_error))
    .trajectory_cp_stop(paste0("private numeric core phase '", phase,
      "' is not runnable: the configured interpreter '", python,
      "' could not be launched (", launched$launch_error, ")."))
  status <- as.integer(attr(output, "status") %||% 0L)
  transcript <- paste(substr(output[seq_len(min(length(output), 12L))], 1L,
    400L), collapse = "\n")
  if (!file.exists(response_path))
    .trajectory_cp_stop(paste0("private numeric core phase '", phase,
      if (identical(status, 127L))
        "' is not runnable: the configured interpreter was not found."
      else paste0("' exited with status ", status,
        " and produced no response. "),
      if (nzchar(transcript)) paste0(" Output:\n", transcript) else ""))
  response <- tryCatch(jsonlite::read_json(response_path,
    simplifyVector = FALSE), error = function(error)
      .trajectory_cp_stop(paste0("private numeric core phase '", phase,
        "' returned unparseable JSON: ", conditionMessage(error))))
  .trajectory_cp_assert(is.list(response) &&
      identical(as.character(response$schema %||% ""), .trajectory_cp_schema) &&
      identical(as.character(response$phase %||% ""), phase),
    paste0("private numeric core phase '", phase,
      "' changed the request schema or phase."))
  if (!identical(response$ok, TRUE)) {
    code <- as.character(as.integer(response$error %||% 5L)[1L])
    .trajectory_cp_stop(paste0("private numeric core phase '", phase,
      "' failed [", .trajectory_cp_codes[[code]] %||% paste0("EXIT_", code),
      "]: ", as.character(response$message %||% "no message")[1L]))
  }
  response
}

# The sort exists in exactly one place: Python never re-orders a vector. Ties
# break on ascending cell_id in the C locale in both directions, so a
# descending axis is rev(ascending) only when no pseudotime tie group exists.
.trajectory_cp_axis <- function(pseudotime, membership, cells, direction,
    compartment) {
  .trajectory_kind(pseudotime, "cell_vector")
  time <- suppressWarnings(as.numeric(pseudotime$value))
  names(time) <- as.character(pseudotime$cell_id)
  selected <- cells
  .trajectory_cp_assert(length(selected) > 0L && !anyDuplicated(selected) &&
      !anyNA(selected), "the effective cell axis is not unique and non-empty.")
  if (!is.null(membership)) {
    .trajectory_kind(membership, "cell_vector")
    label <- as.character(membership$value)
    label[is.na(label)] <- ""
    labels <- unique(sort(label[nzchar(label)]))
    included <- as.character(membership$cell_id)[label == compartment]
    .trajectory_cp_assert(length(included) > 0L, paste0(
      "membership compartment '", compartment, "' names no cell; the labelled",
      " compartments are: ",
      if (length(labels)) paste0("'", paste(labels, collapse = "', '"), "'")
      else "none", "."))
    selected <- intersect(included, selected)
    .trajectory_cp_assert(length(selected) > 0L, paste0(
      "membership compartment '", compartment,
      "' shares no cells with the effective cell axis."))
  }
  value <- suppressWarnings(as.numeric(unname(time[match(selected, names(time))])))
  keep <- !is.na(value)
  dropped <- as.integer(sum(!keep))
  .trajectory_cp_assert(any(keep),
    "no selected cell carries finite pseudotime.")
  cell_id <- selected[keep]
  ordered <- order(cell_id, method = "radix")
  key <- value[ordered]
  if (identical(direction, "descending")) key <- -key
  ranked <- order(key, method = "radix")
  cell_id <- cell_id[ordered][ranked]
  value <- value[ordered][ranked]
  .trajectory_cp_assert(length(cell_id) >= 2L, paste0(
    "the analysed pseudotime axis holds ", length(cell_id),
    " cells; a change point needs at least two."))
  list(cell_id = cell_id, pseudotime = value, n_cells = length(cell_id),
    n_input_cells = length(cells), n_dropped_cells = dropped,
    direction = direction, compartment = compartment)
}

# A declared rank window [rank_start, rank_end) is half-open: rank_start is the
# first included rank and rank_end - 1 the last. The notebook's strict rank
# predicates translate exactly, so `index > 6000 & index < 6500` becomes
# [6001, 6500). An empty declaration analyses the whole axis in both directions.
.trajectory_cp_windows <- function(value, n) {
  if (is.null(value)) value <- list()
  if (is.data.frame(value))
    value <- lapply(seq_len(nrow(value)), function(i)
      stats::setNames(as.list(value[i, , drop = FALSE]), names(value)))
  if (!is.list(value))
    .trajectory_cp_stop("region_windows must be an array of objects.")
  if (length(value) == 4L && all(c("region_id", "direction", "rank_start",
      "rank_end") %in% names(value))) value <- list(value)
  if (!length(value))
    return(lapply(c("up", "down"), function(key)
      list(region_id = "whole_axis", direction = key, rank_start = 0L,
        rank_end = as.integer(n), region_source = "default_whole_axis")))
  fields <- c("region_id", "direction", "rank_start", "rank_end")
  rows <- lapply(seq_along(value), function(i) {
    entry <- value[[i]]
    .trajectory_cp_assert(is.list(entry) && all(fields %in% names(entry)),
      paste0("region_windows[", i, "] must declare region_id, direction, ",
        "rank_start and rank_end."))
    region_id <- as.character(entry$region_id %||% "")[1L]
    direction <- as.character(entry$direction %||% "")[1L]
    start <- suppressWarnings(as.integer(entry$rank_start)[1L])
    end <- suppressWarnings(as.integer(entry$rank_end)[1L])
    .trajectory_cp_assert(!is.na(region_id) && nzchar(region_id),
      paste0("region_windows[", i, "] needs a non-empty region_id."))
    .trajectory_cp_assert(grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$", region_id),
      paste0("region_windows[", i, "] region_id must match ",
        "^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$, got '", region_id, "'."))
    .trajectory_cp_assert(direction %in% c("up", "down"), paste0(
      "region_windows[", i, "] direction must be 'up' or 'down', got '",
      direction, "'."))
    .trajectory_cp_assert(!is.na(start) && !is.na(end) && start >= 0L &&
        end >= start && end <= n, paste0("region_windows[", i,
      "] declares [", start, ", ", end, ") outside the half-open rank ",
      "domain [0, ", n, "]."))
    list(region_id = region_id, direction = direction, rank_start = start,
      rank_end = end, region_source = "declared")
  })
  pairs <- paste(vapply(rows, function(row) row$region_id, character(1)),
    vapply(rows, function(row) row$direction, character(1)), sep = "\034")
  .trajectory_cp_assert(!anyDuplicated(pairs),
    "region_windows repeats a (region_id, direction) pair.")
  rows
}

.trajectory_cp_numeric <- function(value, length_expected, label) {
  result <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  .trajectory_cp_assert(length(result) == length_expected &&
      all(is.finite(result)),
    paste0("the private numeric core returned an invalid ", label, "."))
  result
}

.trajectory_cp_rows <- function(value, fields, label) {
  .trajectory_cp_assert(is.list(value) &&
      all(vapply(value, function(row) is.list(row) &&
        all(fields %in% names(row)), logical(1))),
    paste0("the private numeric core returned invalid ", label, "."))
  value
}

.trajectory_cp_scalar <- function(value, label) {
  .trajectory_cp_assert(!is.null(value) && length(unlist(value)) == 1L,
    paste0("the private numeric core omitted ", label, "."))
  as.character(unlist(value, use.names = FALSE)[1L])
}

.trajectory_cp_text <- function(value) {
  value <- unlist(value, use.names = FALSE)
  if (is.null(value) || !length(value)) return(NA_character_)
  as.character(value)[1L]
}

.trajectory_cp_number <- function(value) {
  value <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  if (!length(value)) return(NA_real_)
  value[1L]
}

# A region's member ids travel as one JSON array inside a character column, so
# the encoding belongs to a JSON writer: hand-joined quoting escaped the double
# quote only, leaving a backslash, a newline or a tab in a feature id able to
# publish an array no parser accepts. auto_unbox stays off so a single member is
# still an array.
.trajectory_cp_ids_json <- function(ids) {
  as.character(jsonlite::toJSON(as.character(unname(ids)), auto_unbox = FALSE))
}

# One wide row, not a long key/value table: retained_parameters addresses
# run_summary by column name, so every recorded fact must be a field.
.trajectory_cp_row <- function(pairs) {
  fields <- stats::setNames(lapply(pairs, function(pair) pair[[2L]][1L]),
    vapply(pairs, function(pair) pair[[1L]], character(1)))
  .trajectory_cp_assert(!anyDuplicated(names(fields)) &&
      all(lengths(fields) == 1L),
    "the run summary row has duplicate or empty fields.")
  as.data.frame(fields, stringsAsFactors = FALSE)
}

.run_trajectory_single_changepoint_regions <- function(readers, writers,
    params, context) {
  direction <- as.character(params$pseudotime_direction %||% "ascending")[1L]
  .trajectory_cp_assert(direction %in% c("ascending", "descending"),
    "pseudotime_direction must be 'ascending' or 'descending'.")
  compartment <- as.character(params$membership_compartment %||% "")[1L]
  smooth_window <- as.integer(params$smooth_window %||% 200L)[1L]
  window_width <- as.integer(params$window_width %||% 100L)[1L]
  window_jump <- as.integer(params$window_jump %||% 5L)[1L]
  window_min_size <- as.integer(params$window_min_size %||% 2L)[1L]
  direction_window <- as.integer(params$direction_window %||% 5L)[1L]
  density_smooth <- as.integer(params$density_smooth_window %||% 100L)[1L]
  density_scale <- as.numeric(params$density_scale %||% 20)[1L]
  change_threshold <- as.numeric(params$change_threshold %||% 30)[1L]
  chunk_size <- as.integer(params$chunk_size %||% 256L)[1L]
  non_finite_policy <- as.character(params$non_finite_policy %||% "error")[1L]
  degenerate_policy <- as.character(params$degenerate_bkp_policy %||%
    "error")[1L]
  for (spec in list(list(smooth_window, "smooth_window", 1L),
      list(window_width, "window_width", 2L),
      list(window_jump, "window_jump", 1L),
      list(window_min_size, "window_min_size", 1L),
      list(direction_window, "direction_window", 1L),
      list(density_smooth, "density_smooth_window", 1L),
      list(chunk_size, "chunk_size", 1L)))
    .trajectory_cp_assert(!is.na(spec[[1L]]) && spec[[1L]] >= spec[[3L]],
      paste0(spec[[2L]], " must be an integer >= ", spec[[3L]], "."))
  .trajectory_cp_assert(window_min_size < window_width,
    "window_min_size must be smaller than window_width.")
  .trajectory_cp_assert(non_finite_policy %in% c("error", "drop"),
    "non_finite_policy must be 'error' or 'drop'.")
  .trajectory_cp_assert(degenerate_policy %in% c("error", "drop"),
    "degenerate_bkp_policy must be 'error' or 'drop'.")
  .trajectory_cp_assert(is.finite(density_scale) && density_scale > 0 &&
      is.finite(change_threshold),
    "density_scale must be positive and change_threshold finite.")

  all_features <- as.character(
    context$typed_slice_reads$matrix_chunks$features)
  features <- .trajectory_feature_sets(readers, params, all_features,
    required = TRUE, prefix = .trajectory_cp_prefix)
  python <- .trajectory_cp_python()
  script <- .trajectory_cp_script()
  script_sha256 <- digest::digest(script, algo = "sha256", file = TRUE)
  matrix_id <- as.character(context$input_ids$matrix)[1L]
  .trajectory_cp_assert(length(matrix_id) == 1L && nzchar(matrix_id),
    "the authenticated matrix ID is unavailable.")
  effective_cells <- as.character(context$effective_cells)
  bound_membership <- is.function(readers[["membership"]])
  .trajectory_cp_assert(nzchar(compartment) == bound_membership, paste0(
    "membership and membership_compartment must be supplied together: ",
    if (bound_membership) "membership is bound but membership_compartment "
      else "membership_compartment is set but membership is ", "is absent."))

  axis <- .trajectory_cp_axis(readers$pseudotime(),
    if (bound_membership) readers$membership() else NULL, effective_cells,
    direction, compartment)
  cells <- axis$cell_id
  n <- axis$n_cells
  regions <- .trajectory_cp_windows(params$region_windows, n)

  probe <- .trajectory_cp_call("probe", list(schema = .trajectory_cp_schema,
    phase = "probe", params = list(window_width = window_width,
      window_jump = window_jump, window_min_size = window_min_size)),
    python, script)
  .trajectory_cp_assert(
    identical(.trajectory_cp_scalar(probe$changepoint_core_sha256,
      "the probed core checksum"), script_sha256),
    paste0("the private numeric core on disk (", script_sha256, ") is not the ",
      "core the interpreter executed (",
      .trajectory_cp_scalar(probe$changepoint_core_sha256, "checksum"), ")."))

  method <- list(smooth_window = smooth_window, window_width = window_width,
    window_jump = window_jump, window_min_size = window_min_size,
    direction_window = direction_window, non_finite_policy = non_finite_policy,
    degenerate_bkp_policy = degenerate_policy)
  partitions <- scworkbench::wb_partition_axis(all_features, chunk_size)
  records <- list()
  for (partition in partitions) {
    block <- readers$matrix_chunks(features = partition,
      cells = effective_cells)$values
    keep <- partition %in% features
    block <- block[keep, , drop = FALSE]
    names <- as.character(rownames(block) %||% partition[keep])
    .trajectory_cp_assert(identical(names, as.character(partition[keep])),
      "an authenticated matrix chunk did not preserve the partition feature axis.")
    index <- match(cells, colnames(block))
    .trajectory_cp_assert(!anyNA(index),
      "the analysed cell axis exceeds the authenticated matrix cell axis.")
    payload <- lapply(seq_along(names), function(i) list(
      feature_id = names[i],
      values = I(as.list(as.numeric(block[i, index, drop = FALSE])))))
    response <- .trajectory_cp_call("features", list(
      schema = .trajectory_cp_schema, phase = "features",
      cells = I(as.list(cells)), features = I(payload), params = method),
      python, script)
    part <- .trajectory_cp_rows(response$records,
      c("feature_id", "b", "delta", "direction", "analysed", "reason",
        "missing_values"), "feature records")
    .trajectory_cp_assert(identical(vapply(part,
      function(row) as.character(row$feature_id), character(1)), names),
      "the private numeric core did not return one record per submitted feature.")
    records <- c(records, part)
  }
  .trajectory_cp_assert(length(records) == length(features), paste0(
    "collected ", length(records), " feature records for ", length(features),
    " selected features."))

  counts <- list(all = integer(n + 1L), up = integer(n + 1L),
    down = integer(n + 1L))
  analysed <- vapply(records, function(row) identical(row$analysed, TRUE),
    logical(1))
  boundary <- integer(length(records))
  for (i in which(analysed)) {
    b <- suppressWarnings(as.integer(records[[i]]$b)[1L])
    .trajectory_cp_assert(!is.na(b) && b > 0L && b < n, paste0("feature '",
      records[[i]]$feature_id, "' returned boundary rank ", b,
      " outside the open interval (0, ", n, ")."))
    boundary[[i]] <- b
    counts$all[[b + 1L]] <- counts$all[[b + 1L]] + 1L
    key <- as.character(records[[i]]$direction %||% "")
    .trajectory_cp_assert(key %in% c("up", "down"), paste0("feature '",
      records[[i]]$feature_id, "' returned direction '", key, "'."))
    counts[[key]][[b + 1L]] <- counts[[key]][[b + 1L]] + 1L
  }
  .trajectory_cp_assert(sum(counts$up) + sum(counts$down) == sum(counts$all),
    "the directional breakpoint counts do not sum to the all-direction count.")

  globals <- .trajectory_cp_call("global", list(schema = .trajectory_cp_schema,
    phase = "global", cell_ids = I(as.list(cells)), counts = list(
      all = I(counts$all), up = I(counts$up), down = I(counts$down)),
    params = list(density_smooth_window = density_smooth,
      density_scale = density_scale, change_threshold = change_threshold)),
    python, script)
  density <- list(all = .trajectory_cp_numeric(globals$density_all, n,
      "all-direction density"),
    up = .trajectory_cp_numeric(globals$density_up, n, "up density"),
    down = .trajectory_cp_numeric(globals$density_down, n, "down density"))
  # Region significance is direction-matched: an up region is gated by the
  # up-direction density and a down region by the down-direction density, as in
  # the notebook. The combined density mixes the two directions, so a strong up
  # population would otherwise carry a lone down breakpoint over the threshold;
  # it stays published, but as a diagnostic.
  .trajectory_cp_assert(identical(.trajectory_cp_scalar(
    globals$significance$region_selection, "the region significance basis"),
    "direction_matched_density_v1"),
    paste0("the private numeric core selected regions on '",
      .trajectory_cp_scalar(globals$significance$region_selection,
        "the region significance basis"),
      "', not on the direction-matched density."))
  significant_ranks <- list()
  significant_set <- list()
  for (key in c("all", "up", "down")) {
    field <- if (identical(key, "all")) "significant_ranks" else
      paste0("significant_ranks_", key)
    ranks <- unique(sort(as.integer(unlist(globals[[field]],
      use.names = FALSE))))
    .trajectory_cp_assert(!anyNA(ranks) && all(ranks >= 0L & ranks < n),
      paste0("the private numeric core returned a ", key,
        "-direction significant rank outside the axis."))
    set <- rep(FALSE, n)
    if (length(ranks)) set[ranks + 1L] <- TRUE
    significant_ranks[[key]] <- ranks
    significant_set[[key]] <- set
  }
  # Diagnostic aliases: the combined density's significance for a feature, and
  # the combined rank set recorded in run_summary.
  significant <- significant_ranks$all
  combined_set <- significant_set$all

  feature_ids <- vapply(records, function(row) as.character(row$feature_id),
    character(1))
  delta <- vapply(records, function(row)
    .trajectory_cp_number(row$delta), numeric(1))
  directions <- vapply(records, function(row)
    .trajectory_cp_text(row$direction), character(1))
  reasons <- vapply(records, function(row)
    .trajectory_cp_text(row$reason), character(1))
  missing_counts <- vapply(records, function(row)
    suppressWarnings(as.integer(unlist(row$missing_values,
      use.names = FALSE)[1L])), integer(1))
  .trajectory_cp_assert(!anyNA(missing_counts) && all(missing_counts >= 0L) &&
      all(missing_counts <= n),
    paste0("the private numeric core returned a missing-value count outside [0, ",
      n, "]."))
  # Significance of a feature on its own terms: the density of the direction the
  # breakpoint actually moved in, at that feature's own boundary rank. This is the
  # gate the regions apply, published per feature so a caller can see why a
  # feature is absent from its region row.
  own_significant <- vapply(seq_along(records), function(i) {
    if (!analysed[[i]]) return(NA)
    significant_set[[directions[[i]]]][[boundary[[i]] + 1L]]
  }, logical(1))
  included <- which(analysed)
  window_first <- ifelse(analysed, pmax(0L, boundary - direction_window),
    NA_integer_)
  window_last <- ifelse(analysed, pmin(n, boundary + direction_window) - 1L,
    NA_integer_)
  at <- function(rank) ifelse(rank >= 1L & rank <= n,
    cells[pmin(pmax(rank, 1L), n)], NA_character_)
  feature_table <- data.frame(
    feature_id = feature_ids,
    analysed = analysed,
    reason = reasons,
    missing_value_count = missing_counts,
    boundary_rank = ifelse(analysed, boundary, NA_integer_),
    boundary_before_cell_id = ifelse(analysed, at(boundary), NA_character_),
    boundary_after_cell_id = ifelse(analysed, at(boundary + 1L), NA_character_),
    change_point_pseudotime = ifelse(analysed,
      axis$pseudotime[pmin(boundary + 1L, n)], NA_real_),
    delta = delta,
    direction = directions,
    density_significant = ifelse(analysed,
      combined_set[pmin(boundary + 1L, n)], NA),
    direction_density_significant = own_significant,
    delta_window_first_rank = window_first,
    delta_window_last_rank = window_last,
    delta_window_first_cell_id = ifelse(analysed, at(window_first + 1L),
      NA_character_),
    delta_window_last_cell_id = ifelse(analysed, at(window_last + 1L),
      NA_character_),
    delta_window_clipped = ifelse(analysed,
      (boundary - direction_window) < 0L |
        (boundary + direction_window) > n, NA),
    stringsAsFactors = FALSE)
  .trajectory_cp_assert(all(feature_table$boundary_rank[included] > 0L),
    "an analysed feature lost its boundary rank in the published table.")

  region_table <- do.call(rbind, c(lapply(regions, function(window) {
    mask <- analysed & directions == window$direction
    members <- character(0)
    if (any(mask)) {
      rank <- boundary[mask]
      member_index <- which(mask)[rank >= window$rank_start &
        rank < window$rank_end &
        significant_set[[window$direction]][rank + 1L]]
      members <- feature_ids[member_index]
    }
    empty <- window$rank_end <= window$rank_start
    data.frame(region_id = window$region_id, direction = window$direction,
      rank_start = window$rank_start, rank_end = window$rank_end,
      rank_count = window$rank_end - window$rank_start,
      first_included_rank = if (empty) NA_integer_ else window$rank_start,
      last_included_rank = if (empty) NA_integer_ else window$rank_end - 1L,
      first_included_cell_id = if (empty) NA_character_ else
        cells[window$rank_start + 1L],
      last_included_cell_id = if (empty) NA_character_ else
        cells[window$rank_end],
      first_included_pseudotime = if (empty) NA_real_ else
        axis$pseudotime[window$rank_start + 1L],
      last_included_pseudotime = if (empty) NA_real_ else
        axis$pseudotime[window$rank_end],
      feature_count = length(members),
      feature_ids_json = .trajectory_cp_ids_json(members),
      significance_basis = paste0(window$direction, "_density"),
      region_source = window$region_source, stringsAsFactors = FALSE)
  }), list()))

  axis_count <- length(effective_cells)
  slots <- match(cells, effective_cells)
  .trajectory_cp_assert(!anyNA(slots),
    "the analysed cell axis is not a subset of the effective cell axis.")
  value <- rep(NA_real_, 3L * axis_count)
  for (block in list(list(key = "all", values = density$all, position = 1L),
      list(key = "up", values = density$up, position = 2L),
      list(key = "down", values = density$down, position = 3L)))
    value[(block$position - 1L) * axis_count + slots] <- block$values
  family <- scworkbench::wb_prim_cell_vector_family(
    rep(c("all", "up", "down"), each = axis_count),
    rep(effective_cells, 3L), value, value_type = "continuous",
    meta = list(matrix_id = matrix_id,
      rank_axis = "analysed pseudotime order over effective_cells_v1"))

  alignment <- globals$alignment %||% list()
  # The boundary-pseudotime statistics describe the detected change points, so
  # they are the min/median/max of the per-feature boundary pseudotimes - the same
  # values the feature table publishes. The cell axis is not a boundary list:
  # indexing it with the feature mask recycles that mask over the axis and reports
  # axis statistics, which is a different quantity wearing the same column name.
  pseudotime_at_boundary <- feature_table$change_point_pseudotime[which(analysed)]
  .trajectory_cp_assert(length(pseudotime_at_boundary) == sum(analysed) &&
      !anyNA(pseudotime_at_boundary),
      "an analysed feature has no boundary pseudotime to summarise.")
  summary <- .trajectory_cp_row(list(
    list("operation_contract", .trajectory_cp_schema),
    list("matrix_id", matrix_id),
    list("n_cells_input", as.integer(axis$n_input_cells)),
    list("n_cells_analysed", as.integer(n)),
    list("n_cells_dropped_non_finite_pseudotime",
      as.integer(axis$n_dropped_cells)),
    list("n_features_selected", as.integer(length(features))),
    list("n_features_analysed", as.integer(sum(analysed))),
    list("n_features_dropped", as.integer(sum(!analysed))),
    list("n_features_with_missing_values", as.integer(sum(missing_counts >
      0L))),
    list("pseudotime_direction", direction),
    list("membership_compartment", if (bound_membership) compartment
      else "all_effective_cells"),
    list("tie_break", "ascending_cell_id_both_directions"),
    list("smooth_window", smooth_window),
    list("window_width", window_width),
    list("window_jump", window_jump),
    list("window_min_size", window_min_size),
    list("direction_window", direction_window),
    list("density_smooth_window", density_smooth),
    list("density_scale", density_scale),
    list("change_threshold", change_threshold),
    list("change_threshold_units",
      "absolute units of a density_scale-scaled breakpoint count; not scale free"),
    list("non_finite_policy", non_finite_policy),
    list("degenerate_bkp_policy", degenerate_policy),
    list("chunk_size", chunk_size),
    list("partition_count", as.integer(length(partitions))),
    list("region_window_count", as.integer(length(regions))),
    list("rank_convention", paste0("boundary rank b sits before the cell at ",
      "0-based rank b; region ranks are half-open [rank_start, rank_end)")),
    list("density_alignment", .trajectory_cp_scalar(alignment$published,
      "the published density convention")),
    list("notebook_raw_density_convention",
      .trajectory_cp_scalar(alignment$notebook_raw,
        "the notebook raw density convention")),
    list("density_shift_formula",
      .trajectory_cp_scalar(alignment$shift, "the density shift formula")),
    list("raw_peak_rank_all", as.integer(.trajectory_cp_scalar(
      alignment$raw_peak_rank_all, "the raw density peak rank"))),
    list("aligned_peak_rank_all", as.integer(.trajectory_cp_scalar(
      alignment$aligned_peak_rank_all, "the aligned density peak rank"))),
    # significant_rank_count/first/last are the combined density's, which is a
    # diagnostic; the two directional counts are the gate the region rows applied.
    list("significant_rank_count", as.integer(length(significant))),
    list("significant_rank_role", "combined_density_diagnostic_only"),
    list("significant_rank_count_up",
      as.integer(length(significant_ranks$up))),
    list("significant_rank_count_down",
      as.integer(length(significant_ranks$down))),
    list("significant_rank_first", if (length(significant))
      as.integer(significant[[1L]]) else NA_integer_),
    list("significant_rank_last", if (length(significant))
      as.integer(significant[[length(significant)]]) else NA_integer_),
    list("region_significance_basis", .trajectory_cp_scalar(
      globals$significance$region_selection, "the region significance basis")),
    list("region_significance_rule", .trajectory_cp_scalar(
      globals$significance$rule, "the region significance rule")),
    list("combined_density_role", .trajectory_cp_scalar(
      globals$significance$combined_density_role,
      "the role of the combined density")),
    list("python_interpreter_option", python),
    list("python_executable", .trajectory_cp_scalar(probe$executable,
      "the interpreter realpath")),
    list("python_executable_input", .trajectory_cp_scalar(
      probe$executable_input, "the interpreter argv path")),
    list("python_version", .trajectory_cp_scalar(probe$python_version,
      "the Python version")),
    list("ruptures_version", .trajectory_cp_scalar(probe$ruptures_version,
      "the ruptures distribution version")),
    list("ruptures_dunder_version",
      .trajectory_cp_scalar(probe$ruptures_dunder_version,
        "the ruptures module version")),
    list("numpy_version", .trajectory_cp_scalar(probe$numpy_version,
      "the numpy version")),
    list("scipy_version", .trajectory_cp_scalar(probe$scipy_version,
      "the scipy version")),
    list("pandas_version", .trajectory_cp_scalar(probe$pandas_version,
      "the pandas version")),
    list("changepoint_core_path", script),
    list("changepoint_core_sha256", script_sha256),
    list("probe_signal", .trajectory_cp_scalar(probe$probe$signal,
      "the probe signal")),
    list("probe_n", as.integer(.trajectory_cp_scalar(probe$probe$probe_n,
      "the probe length"))),
    list("probe_planted_boundary", as.integer(.trajectory_cp_scalar(
      probe$probe$planted_boundary, "the planted step"))),
    list("probe_recovered_boundary", as.integer(.trajectory_cp_scalar(
      probe$probe$recovered_boundary, "the recovered step"))),
    list("probe_terminal_sample_count", as.integer(.trajectory_cp_scalar(
      probe$probe$terminal_sample_count, "the probe terminal sample count"))),
    list("change_point_pseudotime_min", if (length(pseudotime_at_boundary))
      min(pseudotime_at_boundary) else NA_real_),
    list("change_point_pseudotime_median", if (length(pseudotime_at_boundary))
      stats::median(pseudotime_at_boundary) else NA_real_),
    list("change_point_pseudotime_max", if (length(pseudotime_at_boundary))
      max(pseudotime_at_boundary) else NA_real_)))

  writers$write_terminal_primitive("feature_single_changepoints",
    scworkbench::wb_prim_feature_table(feature_table,
      table_type = "trajectory_single_changepoints",
      meta = list(matrix_id = matrix_id, pseudotime_direction = direction,
        n_cells = n, smooth_window = smooth_window, width = window_width,
        jump = window_jump, min_size = window_min_size,
        change_threshold = change_threshold, density_scale = density_scale)))
  writers$write_terminal_primitive("trend_regions",
    scworkbench::wb_prim_obs_table(region_table,
      table_type = "trajectory_single_changepoint_regions",
      entity_kind = "changepoint_region",
      meta = list(matrix_id = matrix_id, n_cells = n,
        rank_convention = "half_open_start_included_end_excluded")))
  writers$write_terminal_primitive("trend_density", family)
  writers$write_terminal_primitive("run_summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "trajectory_single_changepoint_summary",
      entity_kind = "analysis_run"))
  invisible(list(status = "ok"))
}
