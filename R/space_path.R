.trajectory_score_linear <- function(cells, path, anchors = NULL) {
  cells <- as.matrix(cells); path <- as.matrix(path)
  .trajectory_assert(nrow(path) >= 2L && ncol(cells) == ncol(path),
    "Path scoring requires compatible cell and path dimensions.")
  a <- path[-nrow(path), , drop = FALSE]
  b <- path[-1L, , drop = FALSE]
  vector <- b - a
  length <- sqrt(rowSums(vector^2))
  .trajectory_assert(all(length > 0), "Path contains zero-length segments.")
  cumulative <- c(0, cumsum(length)); total <- sum(length)
  best_distance <- rep(Inf, nrow(cells)); best_segment <- integer(nrow(cells))
  best_fraction <- numeric(nrow(cells))
  for (i in seq_len(nrow(a))) {
    delta <- sweep(cells, 2L, a[i, ], "-")
    fraction <- as.numeric((delta %*% vector[i, ]) /
      sum(vector[i, ]^2))
    fraction <- pmin(1, pmax(0, fraction))
    projected <- sweep(matrix(fraction, nrow(cells), ncol(cells)), 2L,
      vector[i, ], "*")
    projected <- sweep(projected, 2L, a[i, ], "+")
    distance <- rowSums((cells - projected)^2)
    take <- distance < best_distance
    best_distance[take] <- distance[take]
    best_segment[take] <- i; best_fraction[take] <- fraction[take]
  }
  progress <- (cumulative[best_segment] +
    best_fraction * length[best_segment]) / total
  distance <- sqrt(best_distance)
  if (is.null(anchors)) anchors <- path[c(1L, nrow(path)), , drop = FALSE]
  anchors <- as.matrix(anchors)
  nearest_anchor <- vapply(seq_len(nrow(cells)), function(i)
    which.min(rowSums((anchors - matrix(cells[i, ], nrow(anchors),
      ncol(cells), byrow = TRUE))^2)), integer(1))
  list(progress = progress, distance = distance,
    nearest_segment = as.numeric(best_segment),
    nearest_anchor = as.numeric(nearest_anchor))
}

.run_space_path_score <- function(readers, writers, params, context) {
  space <- readers$space(); path <- readers$path()
  .trajectory_kind(space, "space"); .trajectory_kind(path,
    "space_path_family")
  path_matrix <- .trajectory_path_matrix(path)
  path_space <- unique(as.character(path$space_id))
  .trajectory_assert(length(path_space) == 1L && identical(path_space,
      as.character(context$input_ids$space)[1L]),
    "space_path_score requires path and scoring space to match.")
  values <- .trajectory_space_matrix(space)
  cells <- as.character(context$effective_cells)
  index <- match(cells, rownames(values))
  .trajectory_assert(!anyNA(index),
    "space_path_score effective cells are absent from the scoring space.")
  values <- values[index, , drop = FALSE]
  scored <- .trajectory_score_linear(values, path_matrix$points,
    .trajectory_anchor_matrix(path))
  prefix <- as.character(params$score_prefix %||% "")[1L]
  key <- function(value) if (nzchar(prefix)) paste0(prefix, value) else value
  output <- scworkbench::wb_prim_cell_vector_family(
    key = rep(vapply(c("progress", "distance", "nearest_segment",
      "nearest_anchor"), key, character(1)), each = length(cells)),
    cell_id = rep(cells, 4L), value = c(scored$progress, scored$distance,
      scored$nearest_segment, scored$nearest_anchor),
    value_type = "continuous")
  summary <- data.frame(path_id = path_matrix$path_id,
    space_id = path_space, fit_mode = as.character(
      path$meta$fit_mode %||% "linear")[1L],
    n_checkpoints = as.integer(path$meta$checkpoint_n %||% NA_integer_)[1L],
    n_path_points = as.integer(nrow(path_matrix$points)),
    cells_scored = as.integer(nrow(values)),
    mean_distance = mean(scored$distance),
    median_distance = stats::median(scored$distance),
    p90_distance = as.numeric(stats::quantile(scored$distance, .9,
      names = FALSE)), stringsAsFactors = FALSE)
  writers$write_terminal_primitive("path_scores", output)
  writers$write_terminal_primitive("fit_summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "space_path_score_summary", entity_kind = "path",
      meta = list(path_id = path_matrix$path_id, space_id = path_space)))
  invisible(list(status = "ok"))
}

.trajectory_local_projection <- function(source, target, path) {
  if (!requireNamespace("RANN", quietly = TRUE)) stop(
    "Package 'RANN' is required for local path interpolation.", call. = FALSE)
  .trajectory_assert(nrow(source) == nrow(target) && nrow(source) >= 2L &&
      ncol(target) == 2L && ncol(source) == ncol(path),
    "Source and target spaces are not aligned for path projection.")
  neighbor <- RANN::nn2(source, path, k = min(20L, nrow(source)))
  output <- matrix(NA_real_, nrow(path), 2L,
    dimnames = list(NULL, c("x", "y")))
  for (i in seq_len(nrow(path))) {
    index <- neighbor$nn.idx[i, ]; distance <- neighbor$nn.dists[i, ]
    if (any(distance == 0)) {
      output[i, ] <- colMeans(target[index[distance == 0], , drop = FALSE])
    } else {
      weight <- 1 / (distance + 1e-8); weight <- weight / sum(weight)
      output[i, ] <- colSums(target[index, , drop = FALSE] * weight)
    }
  }
  output
}

.run_space_path_project_2d <- function(readers, writers, params, context) {
  path <- readers$path(); source <- readers$source_space()
  target <- readers$plot_space()
  .trajectory_kind(path, "space_path_family")
  source_id <- unique(as.character(path$space_id))
  .trajectory_assert(length(source_id) == 1L && identical(source_id,
      as.character(context$input_ids$source_space)[1L]),
    "space_path_project_2d requires source_space to match path$space_id.")
  source <- .trajectory_space_matrix(source)
  target <- .trajectory_space_matrix(target)
  .trajectory_assert(ncol(target) == 2L,
    "space_path_project_2d requires an exactly two-dimensional plot_space.")
  shared <- intersect(rownames(source), rownames(target))
  .trajectory_assert(length(shared) >= 2L,
    "Need at least two cells shared between source_space and plot_space.")
  source <- source[shared, , drop = FALSE]
  target <- target[shared, , drop = FALSE]
  path_matrix <- .trajectory_path_matrix(path)
  coordinates <- .trajectory_local_projection(source, target,
    path_matrix$points)
  if (!is.null(params$approx_points)) {
    retain <- min(as.integer(params$approx_points), nrow(coordinates))
    .trajectory_assert(retain >= 2L,
      "approx_points must be >= 2 when provided.")
    index <- unique(round(seq(1, nrow(coordinates), length.out = retain)))
    coordinates <- coordinates[index, , drop = FALSE]
  }
  curves <- scworkbench::wb_prim_path_2d_family(
    key = rep(path_matrix$path_id, nrow(coordinates)),
    point = seq_len(nrow(coordinates)), x = coordinates[, 1L],
    y = coordinates[, 2L])
  writers$write_path_2d_family("curves", curves)
  invisible(list(status = "ok"))
}
