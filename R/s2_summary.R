# Shared S2 method equations. These consume ordinary vectors reconstructed
# from authenticated typed inputs; no bundle or model-state handle is needed.
.trajectory_s2_bins <- function(pseudotime, n_bins) {
  n_bins <- suppressWarnings(as.integer(n_bins))
  if (!is.finite(n_bins) || n_bins < 5L) n_bins <- 50L
  n_bins <- min(n_bins, max(5L, length(unique(pseudotime))))
  bins <- as.integer(cut(pseudotime, breaks = n_bins,
    include.lowest = TRUE, labels = FALSE))
  .trajectory_assert(all(is.finite(bins)),
    "trajectory_smooth: failed to build pseudotime bins.")
  list(bins = bins, n_bins = as.integer(n_bins))
}

.trajectory_s2_aggregate_one <- function(entity_id, lineage_id, pseudotime,
                                          values, bins, min_cells_per_bin,
                                          value_col = "expr",
                                          summary_fun = "mean") {
  .trajectory_assert(summary_fun %in% c("mean", "median"),
    "trajectory_smooth: summary_fun must be 'mean' or 'median'.")
  data <- data.frame(entity_id = rep(as.character(entity_id),
    length(pseudotime)), lineage_id = rep(as.character(lineage_id),
    length(pseudotime)), pseudotime = as.numeric(pseudotime),
    value = as.numeric(values), bin = as.integer(bins),
    stringsAsFactors = FALSE)
  data <- data[is.finite(data$pseudotime) & is.finite(data$value) &
    is.finite(data$bin), , drop = FALSE]
  if (!nrow(data)) return(NULL)
  fun <- if (identical(summary_fun, "median")) stats::median else base::mean
  group <- list(entity_id = data$entity_id,
    lineage_id = data$lineage_id, bin = data$bin)
  aggregate <- stats::aggregate(data[, c("pseudotime", "value")],
    by = group, FUN = function(z) fun(z, na.rm = TRUE))
  count <- stats::aggregate(data$value, by = group,
    FUN = function(z) as.integer(sum(is.finite(z))))
  names(count)[names(count) == "x"] <- "n_cells"
  by <- c("entity_id", "lineage_id", "bin")
  output <- merge(aggregate, count, by = by, all.x = TRUE, sort = FALSE)
  minimum <- stats::aggregate(data$pseudotime, by = group,
    FUN = function(z) suppressWarnings(min(z, na.rm = TRUE)))
  maximum <- stats::aggregate(data$pseudotime, by = group,
    FUN = function(z) suppressWarnings(max(z, na.rm = TRUE)))
  names(minimum)[names(minimum) == "x"] <- "pseudotime_min"
  names(maximum)[names(maximum) == "x"] <- "pseudotime_max"
  output <- merge(output, minimum, by = by, all.x = TRUE, sort = FALSE)
  output <- merge(output, maximum, by = by, all.x = TRUE, sort = FALSE)
  names(output)[names(output) == "value"] <- value_col
  output <- output[is.finite(output$n_cells) &
    output$n_cells >= as.integer(min_cells_per_bin), , drop = FALSE]
  if (!nrow(output)) return(NULL)
  output[order(output$bin), , drop = FALSE]
}
