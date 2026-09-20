.trajectory_route <- function(pseudotime, membership, route_key, cells,
                              weight_mode, prob_min) {
  .trajectory_kind(pseudotime, "cell_vector")
  .trajectory_kind(membership, "relation_family")
  time <- suppressWarnings(as.numeric(pseudotime$value))
  names(time) <- as.character(pseudotime$cell_id)
  time <- time[is.finite(time) & nzchar(names(time))]
  table <- scworkbench::wb_primitive_to_table(membership,
    include_context = FALSE)
  .trajectory_assert(all(c("key", "cell_id") %in% names(table)),
    "trajectory_route_smooth: branch_membership is malformed.")
  weight <- if ("weight" %in% names(table))
    suppressWarnings(as.numeric(table$weight)) else rep(1, nrow(table))
  table <- data.frame(key = as.character(table$key),
    cell_id = as.character(table$cell_id), weight_use = weight,
    stringsAsFactors = FALSE)
  table$weight_use[!is.finite(table$weight_use)] <- 0
  table <- table[table$key == route_key & table$cell_id %in% cells,
    c("cell_id", "weight_use"), drop = FALSE]
  .trajectory_assert(nrow(table) > 0L, paste0(
    "trajectory_route_smooth: no branch_membership rows found for route_key '",
    route_key, "'."))
  aggregate <- stats::aggregate(weight_use ~ cell_id, table, max)
  names(aggregate)[2L] <- "route_weight"
  aggregate <- aggregate[aggregate$cell_id %in% names(time), , drop = FALSE]
  prob_min <- as.numeric(prob_min)[1L]
  if (!is.finite(prob_min)) prob_min <- 0
  prob_min <- max(0, min(1, prob_min))
  if (identical(weight_mode, "binary_threshold")) {
    aggregate <- aggregate[aggregate$route_weight >= prob_min &
      is.finite(aggregate$route_weight), , drop = FALSE]
    aggregate$route_weight <- 1
  } else {
    aggregate <- aggregate[aggregate$route_weight > 0 &
      aggregate$route_weight >= prob_min & is.finite(aggregate$route_weight),
      , drop = FALSE]
  }
  .trajectory_assert(nrow(aggregate) >= 10L, paste0(
    "trajectory_route_smooth: route '", route_key,
    "' has fewer than 10 cells after weighting/filtering."))
  aggregate$pseudotime <- unname(time[aggregate$cell_id])
  aggregate <- aggregate[is.finite(aggregate$pseudotime), , drop = FALSE]
  .trajectory_assert(nrow(aggregate) >= 10L &&
      stats::sd(aggregate$pseudotime) > 0,
    "trajectory_route_smooth: route has insufficient finite pseudotime variation.")
  aggregate[order(aggregate$pseudotime, aggregate$cell_id),
    c("cell_id", "pseudotime", "route_weight"), drop = FALSE]
}

.trajectory_bins <- function(time, n_bins) {
  n_bins <- as.integer(n_bins)
  if (!is.finite(n_bins) || n_bins < 5L) n_bins <- 50L
  n_bins <- min(n_bins, max(5L, length(unique(time))))
  bins <- as.integer(cut(time, breaks = n_bins, include.lowest = TRUE,
    labels = FALSE))
  .trajectory_assert(all(is.finite(bins)),
    "trajectory_route_smooth: failed to build pseudotime bins.")
  list(bins = bins, n_bins = as.integer(n_bins))
}

.trajectory_route_aggregate <- function(feature, route_key, time, value,
                                        bins, weights, min_cells, summary_fun,
                                        weight_mode) {
  data <- data.frame(entity_id = feature, route_key = route_key,
    pseudotime = as.numeric(time), value = as.numeric(value),
    bin = as.integer(bins), route_weight = as.numeric(weights),
    stringsAsFactors = FALSE)
  data <- data[is.finite(data$pseudotime) & is.finite(data$value) &
    is.finite(data$bin) & is.finite(data$route_weight), , drop = FALSE]
  if (!nrow(data)) return(NULL)
  weighted <- identical(weight_mode, "probability")
  if (weighted && identical(summary_fun, "median")) stop(
    "trajectory_route_smooth: summary_fun='median' is not supported for weight_mode='probability'.",
    call. = FALSE)
  chunks <- split(data, interaction(data$entity_id, data$route_key,
    data$bin, drop = TRUE, lex.order = TRUE))
  rows <- lapply(chunks, function(item) {
    if (nrow(item) < min_cells) return(NULL)
    weight_sum <- sum(item$route_weight)
    if (weighted && (!is.finite(weight_sum) || weight_sum <= 0)) return(NULL)
    fun <- if (identical(summary_fun, "median")) stats::median else mean
    value <- if (weighted) stats::weighted.mean(item$value,
      item$route_weight) else fun(item$value)
    time <- if (weighted) stats::weighted.mean(item$pseudotime,
      item$route_weight) else fun(item$pseudotime)
    data.frame(entity_id = item$entity_id[[1L]],
      route_key = item$route_key[[1L]], bin = as.integer(item$bin[[1L]]),
      pseudotime = as.numeric(time),
      pseudotime_min = min(item$pseudotime),
      pseudotime_max = max(item$pseudotime), expr = as.numeric(value),
      n_cells = as.integer(nrow(item)), route_weight_sum = weight_sum,
      stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  output <- do.call(rbind, rows)
  output[order(output$bin), , drop = FALSE]
}

.run_trajectory_route_smooth <- function(readers, writers, params, context) {
  route_key <- as.character(params$route_key %||% "")[1L]
  .trajectory_assert(nzchar(route_key),
    "trajectory_route_smooth: route_key is required.")
  summary_fun <- as.character(params$summary_fun %||% "mean")[1L]
  .trajectory_assert(summary_fun %in% c("mean", "median"),
    "trajectory_route_smooth: summary_fun must be 'mean' or 'median'.")
  weight_mode <- as.character(params$weight_mode %||% "probability")[1L]
  all_features <- as.character(
    context$typed_slice_reads$matrix_chunks$features)
  features <- .trajectory_feature_sets(readers, params, all_features,
    required = TRUE, prefix = "trajectory_route_smooth")
  route <- .trajectory_route(readers$pseudotime(),
    readers$branch_membership(), route_key,
    as.character(context$effective_cells), weight_mode,
    params$prob_min %||% 0)
  cells <- as.character(route$cell_id)
  matrix <- .trajectory_chunked_matrix(readers$matrix_chunks, all_features,
    features, as.character(context$effective_cells), cells,
    params$chunk_size %||% 128L)
  bin <- .trajectory_bins(route$pseudotime, params$n_bins %||% 50L)
  rows <- lapply(seq_len(nrow(matrix)), function(i)
    .trajectory_route_aggregate(rownames(matrix)[i], route_key,
      route$pseudotime, matrix[i, ], bin$bins, route$route_weight,
      as.integer(params$min_cells_per_bin %||% 5L), summary_fun,
      weight_mode))
  rows <- Filter(Negate(is.null), rows)
  table <- if (length(rows)) do.call(rbind, rows) else data.frame(
    feature_id = character(), route_key = character(), bin = integer(),
    pseudotime = numeric(), pseudotime_min = numeric(),
    pseudotime_max = numeric(), expr = numeric(), n_cells = integer(),
    route_weight_sum = numeric(), label = character(),
    stringsAsFactors = FALSE)
  if (length(rows)) {
    table$feature_id <- as.character(table$entity_id)
    table$label <- table$feature_id
    table <- table[, c("feature_id", "route_key", "bin", "pseudotime",
      "pseudotime_min", "pseudotime_max", "expr", "n_cells",
      "route_weight_sum", "label"), drop = FALSE]
    rownames(table) <- NULL
  }
  matrix_id <- as.character(context$input_ids$matrix)[1L]
  summary <- data.frame(route_key = route_key,
    n_cells_route = as.integer(length(cells)),
    route_weight_sum = sum(route$route_weight),
    n_features_input = as.integer(length(features)),
    n_features_tested = as.integer(length(unique(table$feature_id))),
    n_features_with_rows = as.integer(length(unique(table$feature_id))),
    n_rows = as.integer(nrow(table)), n_bins = bin$n_bins,
    summary_fun = summary_fun,
    min_cells_per_bin = as.integer(params$min_cells_per_bin %||% 5L),
    weight_mode = weight_mode, prob_min = as.numeric(params$prob_min %||% 0),
    matrix_id = matrix_id, stringsAsFactors = FALSE)
  writers$write_terminal_primitive("trend_table",
    scworkbench::wb_prim_obs_table(table,
      table_type = "trajectory_route_smooth",
      entity_kind = "trajectory_bin", meta = list(matrix_id = matrix_id,
        route_key = route_key, summary_fun = summary_fun,
        min_cells_per_bin = as.integer(params$min_cells_per_bin %||% 5L),
        n_bins = bin$n_bins, weight_mode = weight_mode,
        prob_min = as.numeric(params$prob_min %||% 0))))
  writers$write_terminal_primitive("fit_summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "trajectory_route_smooth_summary", entity_kind = "route"))
  invisible(list(status = "ok"))
}
