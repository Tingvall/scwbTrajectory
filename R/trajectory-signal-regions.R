# Permutation-qualified local signal regions over explicit fitted pseudotime.
.trajectory_regions_optional_reader <- function(readers, name) {
  reader <- readers[[name]]
  if (is.function(reader)) reader() else NULL
}

.trajectory_regions_keys <- function(value, available, label) {
  if (is.null(value)) return(as.character(unique(available)))
  keys <- unique(as.character(unlist(value, use.names = FALSE)))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  .trajectory_assert(length(keys) > 0L,
    paste0("trajectory_signal_regions: ", label, " is empty."))
  .trajectory_assert(!length(setdiff(keys, available)),
    paste0("trajectory_signal_regions: requested ", label,
      " was not found."))
  keys
}

.trajectory_regions_continuous <- function(value, label) {
  .trajectory_kind(value, "cell_vector_family")
  value_type <- as.character(value$value_type %||% "continuous")[1L]
  .trajectory_assert(!identical(value_type, "categorical"),
    paste0("trajectory_signal_regions: ", label,
      " must contain continuous values."))
  table <- scworkbench::wb_primitive_to_table(value,
    include_context = FALSE)
  .trajectory_assert(all(c("key", "cell_id", "value") %in% names(table)),
    paste0("trajectory_signal_regions: ", label,
      " table is malformed."))
  table <- data.frame(key = as.character(table$key),
    cell_id = as.character(table$cell_id),
    value = suppressWarnings(as.numeric(table$value)),
    stringsAsFactors = FALSE)
  .trajectory_assert(!anyNA(table$key) && !anyNA(table$cell_id) &&
      all(nzchar(table$key)) && all(nzchar(table$cell_id)),
    paste0("trajectory_signal_regions: ", label,
      " has an empty key or cell ID."))
  pair <- paste(table$key, table$cell_id, sep = "\r")
  .trajectory_assert(!anyDuplicated(pair),
    paste0("trajectory_signal_regions: ", label,
      " has duplicate key/cell rows."))
  table
}

.trajectory_regions_summary <- function(x, method) {
  if (identical(method, "median")) stats::median(x) else mean(x)
}

.trajectory_regions_max_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

.trajectory_regions_seed <- function(seed, lineage, signal) {
  bytes <- utf8ToInt(enc2utf8(paste(lineage, signal, sep = "\r")))
  offset <- if (length(bytes)) sum((as.double(bytes) *
    seq_along(bytes)) %% 1000003) else 0
  as.integer((as.double(seed) + offset) %% 2147483646 + 1)
}

.trajectory_regions_with_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv,
    inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

.trajectory_regions_permute <- function(values, blocks = NULL) {
  if (is.null(blocks)) return(sample(values, length(values), replace = FALSE))
  output <- values
  groups <- split(seq_along(values), as.character(blocks), drop = TRUE)
  for (index in groups) {
    if (length(index) > 1L)
      output[index] <- sample(values[index], length(index), replace = FALSE)
  }
  output
}

.trajectory_regions_extrema <- function(values, direction,
                                        include_endpoints = TRUE) {
  values <- as.numeric(values)
  if (length(values) < 2L || length(unique(values)) < 2L)
    return(integer(0))
  runs <- rle(values)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  selected <- logical(length(runs$values))
  for (ii in seq_along(runs$values)) {
    if (!include_endpoints && ii %in% c(1L, length(runs$values))) next
    left <- if (ii == 1L) NA_real_ else runs$values[ii - 1L]
    right <- if (ii == length(runs$values)) NA_real_ else
      runs$values[ii + 1L]
    if (identical(direction, "peak")) {
      selected[ii] <- (is.na(left) || runs$values[ii] > left) &&
        (is.na(right) || runs$values[ii] > right)
    } else {
      selected[ii] <- (is.na(left) || runs$values[ii] < left) &&
        (is.na(right) || runs$values[ii] < right)
    }
  }
  as.integer(floor((starts[selected] + ends[selected]) / 2))
}

.trajectory_regions_windows <- function(pseudotime, samples,
                                        width_fraction, step_fraction,
                                        min_cells, min_samples) {
  pt_min <- min(pseudotime)
  pt_max <- max(pseudotime)
  pt_range <- pt_max - pt_min
  .trajectory_assert(is.finite(pt_range) && pt_range > 0,
    "trajectory_signal_regions: pseudotime has zero range.")
  width <- width_fraction * pt_range
  step <- step_fraction * pt_range
  centers <- seq(pt_min, pt_max, by = step)
  if (utils::tail(centers, 1L) < pt_max - sqrt(.Machine$double.eps) *
      max(1, abs(pt_max))) centers <- c(centers, pt_max)
  centers <- unique(pmin(centers, pt_max))
  left <- pmax(pt_min, centers - width / 2)
  right <- pmin(pt_max, centers + width / 2)
  index <- lapply(seq_along(centers), function(ii)
    which(pseudotime >= left[ii] & pseudotime <= right[ii]))
  n_cells <- lengths(index)
  n_samples <- if (is.null(samples)) rep(NA_integer_, length(index)) else
    vapply(index, function(ii) length(unique(samples[ii])), integer(1))
  keep <- n_cells >= min_cells
  if (!is.null(samples)) keep <- keep & n_samples >= min_samples
  list(center = centers[keep], left = left[keep], right = right[keep],
    index = index[keep], n_cells = as.integer(n_cells[keep]),
    n_samples = as.integer(n_samples[keep]), range = pt_range,
    width = width, step = step)
}

.trajectory_regions_adjust <- function(values, covariates, method,
                                       spline_df) {
  if (identical(method, "none")) return(list(values = values,
    rank = 1L, columns = 1L, r_squared = NA_real_,
    raw_abs_spearman = NA_real_, adjusted_abs_spearman = NA_real_,
    model = "none"))
  pieces <- list(`(Intercept)` = rep(1, length(values)))
  for (key in colnames(covariates)) {
    x <- as.numeric(covariates[, key])
    if (identical(method, "linear")) {
      piece <- matrix(x, ncol = 1L,
        dimnames = list(NULL, key))
    } else {
      .trajectory_assert(length(unique(x)) > spline_df,
        paste0("trajectory_signal_regions: covariate '", key,
          "' has too few unique values for natural_spline."))
      piece <- splines::ns(x, df = spline_df)
      colnames(piece) <- paste0(key, "_ns", seq_len(ncol(piece)))
    }
    pieces[[length(pieces) + 1L]] <- piece
  }
  design <- do.call(cbind, pieces)
  fit <- stats::lm.fit(design, values)
  .trajectory_assert(identical(as.integer(fit$rank), ncol(design)),
    "trajectory_signal_regions: covariate design is rank deficient.")
  adjusted <- as.numeric(fit$residuals)
  total <- sum((values - mean(values))^2)
  r_squared <- if (total > 0) 1 - sum(adjusted^2) / total else 0
  raw_cor <- vapply(seq_len(ncol(covariates)), function(ii)
    suppressWarnings(abs(stats::cor(values, covariates[, ii],
      method = "spearman"))), numeric(1))
  adjusted_cor <- vapply(seq_len(ncol(covariates)), function(ii)
    suppressWarnings(abs(stats::cor(adjusted, covariates[, ii],
      method = "spearman"))), numeric(1))
  list(values = adjusted, rank = as.integer(fit$rank),
    columns = as.integer(ncol(design)), r_squared = as.numeric(r_squared),
    raw_abs_spearman = .trajectory_regions_max_finite(raw_cor),
    adjusted_abs_spearman = .trajectory_regions_max_finite(adjusted_cor),
    model = if (identical(method, "linear"))
      "signal ~ covariates" else paste0(
        "signal ~ ns(each covariate, df=", spline_df, ")"))
}

.trajectory_regions_dynamic_row <- function(signal_key, lineage_id,
                                            status, n_cells, n_samples,
                                            n_windows, seed,
                                            permutation_mode) {
  data.frame(signal_key = signal_key, feature_set = signal_key,
    lineage_id = lineage_id, status = status,
    n_cells = as.integer(n_cells), n_samples = as.integer(n_samples),
    n_windows = as.integer(n_windows), baseline = NA_real_,
    observed_dynamic_range = NA_real_, peak_pseudotime = NA_real_,
    trough_pseudotime = NA_real_, peak_value = NA_real_,
    trough_value = NA_real_, dynamic_exceedances = NA_integer_,
    permutations_used = NA_integer_, p_dynamic = NA_real_,
    q_dynamic = NA_real_, dynamic_significant = FALSE,
    null_range_mean = NA_real_, null_range_sd = NA_real_,
    null_range_q95 = NA_real_, null_range_q99 = NA_real_,
    permutation_seed = as.integer(seed),
    permutation_mode = permutation_mode, n_blocks = NA_integer_,
    n_singleton_blocks = NA_integer_, stringsAsFactors = FALSE)
}

.trajectory_regions_empty_candidates <- function() data.frame(
  candidate_id = character(), signal_key = character(),
  feature_set = character(), lineage_id = character(),
  direction = character(), window_index = integer(),
  peak_pseudotime = numeric(), interval_start = numeric(),
  interval_end = numeric(), window_value = numeric(), baseline = numeric(),
  effect_size = numeric(), raw_effect_size = numeric(), n_cells = integer(),
  n_samples = integer(), region_exceedances = integer(),
  permutations_used = integer(), p_empirical = numeric(),
  q_adjusted = numeric(), p_dynamic = numeric(), q_dynamic = numeric(),
  dynamic_significant = logical(), region_supported = logical(),
  merged_region_id = character(), permutation_seed = integer(),
  permutation_mode = character(), stringsAsFactors = FALSE)

.trajectory_regions_empty_regions <- function() data.frame(
  region_id = character(), signal_key = character(),
  feature_set = character(), lineage_id = character(),
  direction = character(), interval_start = numeric(),
  interval_end = numeric(), peak_pseudotime = numeric(),
  effect_size = numeric(), raw_effect_size = numeric(),
  p_empirical = numeric(), q_adjusted = numeric(),
  p_dynamic = numeric(), q_dynamic = numeric(), n_cells = integer(),
  n_samples = integer(), peak_n_cells = integer(),
  peak_n_samples = integer(), n_constituent_windows = integer(),
  constituent_candidate_ids = character(), permutations_used = integer(),
  permutation_seed = integer(), permutation_mode = character(),
  stringsAsFactors = FALSE)

.trajectory_signal_regions <- function(effective_cells, signal, pseudotime,
                                      pseudotime_key, covariates = NULL,
                                      sample_id = NULL, membership = NULL,
                                      depth = NULL,
                                      signal_keys = NULL,
                                      covariate_keys = NULL,
                                      pseudotime_direction = "as_is",
                                      window_width_fraction = 0.05,
                                      window_step_fraction = 0.01,
                                      min_cells_per_window = 20L,
                                      min_samples_per_window = 1L,
                                      summary_fun = "median",
                                      n_permutations = 1000L,
                                      seed = 1L,
                                      covariate_adjustment = "none",
                                      covariate_spline_df = 4L,
                                      adjust_depth = TRUE,
                                      depth_field = "nCount_RNA",
                                      permutation_mode = "auto",
                                      include_endpoint_extrema = TRUE,
                                      dynamic_q_cutoff = 0.05,
                                      region_q_cutoff = 0.05,
                                      merge_gap_fraction = NULL,
                                      signal_artifact_id = "",
                                      pseudotime_artifact_id = "",
                                      covariates_artifact_id = "",
                                      depth_artifact_id = "",
                                      sample_id_artifact_id = "",
                                      membership_artifact_id = "") {
  pseudotime_key <- as.character(pseudotime_key %||% "")[1L]
  .trajectory_assert(!is.na(pseudotime_key) && nzchar(trimws(
    pseudotime_key)),
    "trajectory_signal_regions: pseudotime_key is required.")
  pseudotime_key <- trimws(pseudotime_key)
  pseudotime_direction <- as.character(pseudotime_direction)[1L]
  summary_fun <- as.character(summary_fun)[1L]
  covariate_adjustment <- as.character(covariate_adjustment)[1L]
  adjust_depth <- isTRUE(adjust_depth)
  depth_field <- trimws(as.character(depth_field %||% "nCount_RNA")[1L])
  permutation_mode <- as.character(permutation_mode)[1L]
  .trajectory_assert(pseudotime_direction %in% c("as_is", "invert"),
    "trajectory_signal_regions: invalid pseudotime_direction.")
  .trajectory_assert(summary_fun %in% c("mean", "median"),
    "trajectory_signal_regions: summary_fun must be mean or median.")
  .trajectory_assert(covariate_adjustment %in%
      c("none", "linear", "natural_spline"),
    "trajectory_signal_regions: invalid covariate_adjustment.")
  .trajectory_assert(!is.na(depth_field) && nzchar(depth_field),
    "trajectory_signal_regions: depth_field is required.")
  .trajectory_assert(permutation_mode %in%
      c("auto", "blocked", "unblocked"),
    "trajectory_signal_regions: invalid permutation_mode.")

  width_fraction <- suppressWarnings(as.numeric(window_width_fraction)[1L])
  step_fraction <- suppressWarnings(as.numeric(window_step_fraction)[1L])
  min_cells <- suppressWarnings(as.integer(min_cells_per_window)[1L])
  min_samples <- suppressWarnings(as.integer(min_samples_per_window)[1L])
  permutations <- suppressWarnings(as.integer(n_permutations)[1L])
  seed <- suppressWarnings(as.integer(seed)[1L])
  spline_df <- suppressWarnings(as.integer(covariate_spline_df)[1L])
  dynamic_q_cutoff <- suppressWarnings(as.numeric(dynamic_q_cutoff)[1L])
  region_q_cutoff <- suppressWarnings(as.numeric(region_q_cutoff)[1L])
  include_endpoint_extrema <- isTRUE(include_endpoint_extrema)
  .trajectory_assert(is.finite(width_fraction) && width_fraction > 0 &&
      width_fraction <= 1, "trajectory_signal_regions: invalid window width.")
  .trajectory_assert(is.finite(step_fraction) && step_fraction > 0 &&
      step_fraction <= 1, "trajectory_signal_regions: invalid window step.")
  .trajectory_assert(is.finite(min_cells) && min_cells >= 2L,
    "trajectory_signal_regions: min_cells_per_window must be >= 2.")
  .trajectory_assert(is.finite(min_samples) && min_samples >= 1L,
    "trajectory_signal_regions: min_samples_per_window must be >= 1.")
  .trajectory_assert(is.finite(permutations) && permutations >= 1L,
    "trajectory_signal_regions: n_permutations must be >= 1.")
  .trajectory_assert(is.finite(seed) && seed >= 1L,
    "trajectory_signal_regions: seed must be >= 1.")
  .trajectory_assert(is.finite(spline_df) && spline_df >= 2L,
    "trajectory_signal_regions: covariate_spline_df must be >= 2.")
  .trajectory_assert(is.finite(dynamic_q_cutoff) &&
      dynamic_q_cutoff > 0 && dynamic_q_cutoff <= 1 &&
      is.finite(region_q_cutoff) && region_q_cutoff > 0 &&
      region_q_cutoff <= 1,
    "trajectory_signal_regions: q-value cutoffs must lie in (0, 1].")
  if (is.null(merge_gap_fraction)) merge_gap_fraction <- width_fraction
  merge_gap_fraction <- suppressWarnings(as.numeric(
    merge_gap_fraction)[1L])
  .trajectory_assert(is.finite(merge_gap_fraction) &&
      merge_gap_fraction >= 0 && merge_gap_fraction <= 1,
    "trajectory_signal_regions: invalid merge_gap_fraction.")

  signal_table <- .trajectory_regions_continuous(signal, "signal")
  already_depth_adjusted <- isTRUE(signal$meta$depth_adjusted)
  if (adjust_depth && already_depth_adjusted) {
    prior_field <- as.character(signal$meta$depth_field %||% "")[1L]
    .trajectory_assert(!nzchar(prior_field) || identical(prior_field,
        depth_field),
      paste0("trajectory_signal_regions: signal was depth-adjusted using ",
        "a different depth field."))
  }
  pseudotime_table <- .trajectory_regions_continuous(pseudotime,
    "pseudotime")
  available_signal_keys <- unique(signal_table$key)
  signal_keys <- .trajectory_regions_keys(signal_keys,
    available_signal_keys, "signal_keys")
  .trajectory_assert(pseudotime_key %in% pseudotime_table$key,
    "trajectory_signal_regions: pseudotime_key was not found.")

  cells <- unique(as.character(effective_cells))
  .trajectory_assert(length(cells) > 0L && !anyNA(cells) &&
      all(nzchar(cells)) && !anyDuplicated(cells),
    "trajectory_signal_regions: effective cell axis is invalid.")
  if (!is.null(membership)) {
    .trajectory_kind(membership, "cell_vector")
    member_cells <- as.character(membership$cell_id)
    member_value <- as.character(membership$value)
    .trajectory_assert(!anyDuplicated(member_cells),
      "trajectory_signal_regions: membership has duplicate cells.")
    member_cells <- member_cells[!is.na(member_value) &
      nzchar(member_value)]
    cells <- cells[cells %in% member_cells]
  }
  depth_values <- NULL
  if (adjust_depth && !already_depth_adjusted) {
    .trajectory_assert(!is.null(depth), paste0(
      "trajectory_signal_regions: adjust_depth=TRUE requires an authenticated ",
      "depth input unless the signal is already depth-adjusted."))
    depth_values <- .trajectory_depth_vector(depth, effective_cells,
      "trajectory_signal_regions")
  }

  pt_all <- pseudotime_table[pseudotime_table$key == pseudotime_key &
    pseudotime_table$cell_id %in% cells, , drop = FALSE]
  pt_finite <- pt_all[is.finite(pt_all$value), , drop = FALSE]
  .trajectory_assert(nrow(pt_finite) >= 3L,
    "trajectory_signal_regions: selected pseudotime has too few finite cells.")
  pt <- pt_finite$value
  names(pt) <- pt_finite$cell_id
  if (identical(pseudotime_direction, "invert"))
    pt <- max(pt) - pt

  covariate_table <- NULL
  if (!is.null(covariates)) {
    covariate_table <- .trajectory_regions_continuous(covariates,
      "covariates")
    covariate_keys <- .trajectory_regions_keys(covariate_keys,
      unique(covariate_table$key), "covariate_keys")
  }
  .trajectory_assert((identical(covariate_adjustment, "none") &&
      is.null(covariate_table)) || (!identical(covariate_adjustment,
      "none") && !is.null(covariate_table)),
    paste0("trajectory_signal_regions: covariates and ",
      "covariate_adjustment must be supplied together."))

  sample_values <- NULL
  if (!is.null(sample_id)) {
    .trajectory_kind(sample_id, "cell_vector")
    sample_cells <- as.character(sample_id$cell_id)
    sample_values_all <- as.character(sample_id$value)
    .trajectory_assert(!anyDuplicated(sample_cells) &&
        !anyNA(sample_cells) && all(nzchar(sample_cells)),
      "trajectory_signal_regions: sample_id has invalid cell IDs.")
    .trajectory_assert(!anyNA(sample_values_all) &&
        all(nzchar(sample_values_all)),
      "trajectory_signal_regions: sample_id has missing block labels.")
    sample_values <- stats::setNames(sample_values_all, sample_cells)
  }
  actual_permutation_mode <- if (identical(permutation_mode, "auto")) {
    if (is.null(sample_values)) "unblocked" else "blocked"
  } else permutation_mode
  .trajectory_assert(!identical(actual_permutation_mode, "blocked") ||
      !is.null(sample_values),
    "trajectory_signal_regions: blocked permutations require sample_id.")
  .trajectory_assert(!is.null(sample_values) || min_samples == 1L,
    paste0("trajectory_signal_regions: min_samples_per_window > 1 ",
      "requires sample_id."))

  curve_rows <- list()
  dynamic_rows <- list()
  candidate_rows <- list()
  adjustment_rows <- list()
  analysis_data <- list()
  curve_n <- dynamic_n <- candidate_n <- adjustment_n <- 0L
  analyzable <- 0L

  for (signal_key in signal_keys) {
    signal_sub <- signal_table[signal_table$key == signal_key &
      signal_table$cell_id %in% names(pt), , drop = FALSE]
    signal_value <- stats::setNames(signal_sub$value, signal_sub$cell_id)
    shared <- names(pt)[names(pt) %in% names(signal_value)]
    raw <- as.numeric(signal_value[shared])
    time <- as.numeric(pt[shared])
    samples <- if (is.null(sample_values)) NULL else
      as.character(sample_values[shared])
    signal_depth <- if (is.null(depth_values)) NULL else
      as.numeric(depth_values[shared])
    covariate_matrix <- NULL
    if (!is.null(covariate_table)) {
      covariate_matrix <- vapply(covariate_keys, function(key) {
        sub <- covariate_table[covariate_table$key == key, , drop = FALSE]
        value <- stats::setNames(sub$value, sub$cell_id)
        as.numeric(value[shared])
      }, numeric(length(shared)))
      if (is.null(dim(covariate_matrix))) covariate_matrix <- matrix(
        covariate_matrix, ncol = 1L)
      colnames(covariate_matrix) <- covariate_keys
    }
    complete <- is.finite(raw) & is.finite(time)
    if (!is.null(covariate_matrix)) complete <- complete &
      apply(covariate_matrix, 1L, function(x) all(is.finite(x)))
    if (!is.null(signal_depth)) complete <- complete &
      is.finite(signal_depth) & signal_depth >= 0
    if (!is.null(samples)) complete <- complete & !is.na(samples) &
      nzchar(samples)
    shared <- shared[complete]
    raw <- raw[complete]
    time <- time[complete]
    if (!is.null(samples)) samples <- samples[complete]
    if (!is.null(signal_depth)) signal_depth <- signal_depth[complete]
    if (!is.null(covariate_matrix)) covariate_matrix <-
      covariate_matrix[complete, , drop = FALSE]
    effective_seed <- .trajectory_regions_seed(seed, pseudotime_key,
      signal_key)

    if (length(raw) < min_cells) {
      dynamic_n <- dynamic_n + 1L
      dynamic_rows[[dynamic_n]] <- .trajectory_regions_dynamic_row(
        signal_key, pseudotime_key, "insufficient_cells", length(raw),
        if (is.null(samples)) NA_integer_ else length(unique(samples)),
        0L, effective_seed, actual_permutation_mode)
      adjustment_n <- adjustment_n + 1L
      adjustment_rows[[adjustment_n]] <- data.frame(
        signal_key = signal_key, feature_set = signal_key,
        lineage_id = pseudotime_key, status = "insufficient_cells",
        adjustment_method = covariate_adjustment,
        adjustment_model = if (identical(covariate_adjustment, "none"))
          "none" else "not_fitted", n_cells_complete = length(raw),
        n_covariates = length(covariate_keys %||% character()),
        design_rank = NA_integer_, design_columns = NA_integer_,
        model_r_squared = NA_real_, raw_max_abs_spearman = NA_real_,
        adjusted_max_abs_spearman = NA_real_,
        depth_adjustment_status = if (!adjust_depth) "disabled" else
          if (already_depth_adjusted) "already_adjusted" else "not_fitted",
        depth_field = depth_field, depth_artifact_id = depth_artifact_id,
        depth_coefficient = NA_real_, depth_partial_r_squared = NA_real_,
        depth_spline_df = if (adjust_depth) 4L else NA_integer_,
        stringsAsFactors = FALSE)
      next
    }

    depth_adjustment <- NULL
    depth_adjusted <- raw
    if (adjust_depth && !already_depth_adjusted) {
      depth_adjustment <- .trajectory_depth_adjust_vector(raw, time,
        signal_depth, spline_df = 4L,
        label = paste0("trajectory_signal_regions signal '", signal_key,
          "'"))
      depth_adjusted <- depth_adjustment$values
    }
    adjusted <- .trajectory_regions_adjust(depth_adjusted, covariate_matrix,
      covariate_adjustment, spline_df)
    values <- adjusted$values
    windows <- .trajectory_regions_windows(time, samples,
      width_fraction, step_fraction, min_cells, min_samples)
    adjustment_n <- adjustment_n + 1L
    adjustment_rows[[adjustment_n]] <- data.frame(
      signal_key = signal_key, feature_set = signal_key,
      lineage_id = pseudotime_key,
      status = if (length(windows$index) >= 2L) "ok" else
        "insufficient_windows", adjustment_method = covariate_adjustment,
      adjustment_model = adjusted$model, n_cells_complete = length(raw),
      n_covariates = length(covariate_keys %||% character()),
      design_rank = adjusted$rank, design_columns = adjusted$columns,
      model_r_squared = adjusted$r_squared,
      raw_max_abs_spearman = adjusted$raw_abs_spearman,
      adjusted_max_abs_spearman = adjusted$adjusted_abs_spearman,
      depth_adjustment_status = if (!adjust_depth) "disabled" else
        if (already_depth_adjusted) "already_adjusted" else "applied",
      depth_field = depth_field, depth_artifact_id = depth_artifact_id,
      depth_coefficient = if (is.null(depth_adjustment)) NA_real_ else
        depth_adjustment$diagnostics$depth_coefficient[[1L]],
      depth_partial_r_squared = if (is.null(depth_adjustment)) NA_real_ else
        depth_adjustment$diagnostics$depth_partial_r_squared[[1L]],
      depth_spline_df = if (adjust_depth) 4L else NA_integer_,
      stringsAsFactors = FALSE)
    if (length(windows$index) < 2L) {
      dynamic_n <- dynamic_n + 1L
      dynamic_rows[[dynamic_n]] <- .trajectory_regions_dynamic_row(
        signal_key, pseudotime_key, "insufficient_windows", length(raw),
        if (is.null(samples)) NA_integer_ else length(unique(samples)),
        length(windows$index), effective_seed, actual_permutation_mode)
      next
    }

    if (identical(actual_permutation_mode, "blocked")) {
      sizes <- table(samples)
      .trajectory_assert(any(sizes >= 2L), paste0(
        "trajectory_signal_regions: signal '", signal_key,
        "' has no permutable sample block."))
      blocks <- samples
      n_blocks <- length(sizes)
      n_singletons <- sum(sizes == 1L)
    } else {
      blocks <- NULL
      n_blocks <- if (is.null(samples)) NA_integer_ else
        length(unique(samples))
      n_singletons <- if (is.null(samples)) NA_integer_ else
        sum(table(samples) == 1L)
    }

    window_value <- vapply(windows$index, function(ii)
      .trajectory_regions_summary(values[ii], summary_fun), numeric(1))
    raw_window_value <- vapply(windows$index, function(ii)
      .trajectory_regions_summary(raw[ii], summary_fun), numeric(1))
    baseline <- .trajectory_regions_summary(values, summary_fun)
    raw_baseline <- .trajectory_regions_summary(raw, summary_fun)
    deviation <- window_value - baseline
    observed_range <- max(window_value) - min(window_value)
    peak_index <- .trajectory_regions_extrema(window_value, "peak",
      include_endpoint_extrema)
    trough_index <- .trajectory_regions_extrema(window_value, "trough",
      include_endpoint_extrema)
    peak_index <- peak_index[deviation[peak_index] > 0]
    trough_index <- trough_index[deviation[trough_index] < 0]

    null <- .trajectory_regions_with_seed(effective_seed, {
      range_null <- positive_null <- negative_null <-
        numeric(permutations)
      for (permutation in seq_len(permutations)) {
        permuted <- .trajectory_regions_permute(values, blocks)
        perm_window <- vapply(windows$index, function(ii)
          .trajectory_regions_summary(permuted[ii], summary_fun),
          numeric(1))
        range_null[permutation] <- max(perm_window) - min(perm_window)
        perm_deviation <- perm_window - baseline
        positive_null[permutation] <- max(perm_deviation)
        negative_null[permutation] <- max(-perm_deviation)
      }
      list(range = range_null, positive = positive_null,
        negative = negative_null)
    })
    dynamic_exceedances <- sum(null$range >= observed_range)
    p_dynamic <- (dynamic_exceedances + 1) / (permutations + 1)
    analyzable <- analyzable + 1L
    dynamic_n <- dynamic_n + 1L
    dynamic_rows[[dynamic_n]] <- data.frame(
      signal_key = signal_key, feature_set = signal_key,
      lineage_id = pseudotime_key, status = "ok",
      n_cells = length(raw), n_samples = if (is.null(samples))
        NA_integer_ else length(unique(samples)),
      n_windows = length(windows$index), baseline = baseline,
      observed_dynamic_range = observed_range,
      peak_pseudotime = windows$center[which.max(window_value)],
      trough_pseudotime = windows$center[which.min(window_value)],
      peak_value = max(window_value), trough_value = min(window_value),
      dynamic_exceedances = dynamic_exceedances,
      permutations_used = permutations, p_dynamic = p_dynamic,
      q_dynamic = NA_real_, dynamic_significant = FALSE,
      null_range_mean = mean(null$range),
      null_range_sd = stats::sd(null$range),
      null_range_q95 = as.numeric(stats::quantile(null$range, 0.95,
        names = FALSE)), null_range_q99 = as.numeric(stats::quantile(
        null$range, 0.99, names = FALSE)),
      permutation_seed = effective_seed,
      permutation_mode = actual_permutation_mode,
      n_blocks = as.integer(n_blocks),
      n_singleton_blocks = as.integer(n_singletons),
      stringsAsFactors = FALSE)

    local_direction <- rep("", length(window_value))
    local_direction[peak_index] <- "elevated"
    local_direction[trough_index] <- "reduced"
    curve_n <- curve_n + 1L
    curve_rows[[curve_n]] <- data.frame(
      signal_key = signal_key, feature_set = signal_key,
      lineage_id = pseudotime_key,
      bin = seq_along(window_value), window_index = seq_along(window_value),
      pseudotime = windows$center, pseudotime_min = windows$left,
      pseudotime_max = windows$right, value = window_value,
      value_raw = raw_window_value, value_smooth = window_value,
      value_se = NA_real_, value_ci_low = NA_real_,
      value_ci_high = NA_real_, n_cells = windows$n_cells,
      n_samples = windows$n_samples, baseline = baseline,
      raw_baseline = raw_baseline, deviation = deviation,
      local_direction = local_direction, p_empirical = NA_real_,
      q_adjusted = NA_real_, region_supported = FALSE,
      region_id = "", label = signal_key, stringsAsFactors = FALSE)

    local_index <- c(peak_index, trough_index)
    local_type <- c(rep("elevated", length(peak_index)),
      rep("reduced", length(trough_index)))
    if (length(local_index)) {
      for (jj in seq_along(local_index)) {
        ii <- local_index[jj]
        direction <- local_type[jj]
        effect <- if (identical(direction, "elevated"))
          deviation[ii] else -abs(deviation[ii])
        null_deviation <- if (identical(direction, "elevated"))
          null$positive else null$negative
        observed_deviation <- abs(deviation[ii])
        exceedances <- sum(null_deviation >= observed_deviation)
        candidate_n <- candidate_n + 1L
        candidate_rows[[candidate_n]] <- data.frame(
          candidate_id = sprintf("candidate_%06d", candidate_n),
          signal_key = signal_key, feature_set = signal_key,
          lineage_id = pseudotime_key, direction = direction,
          window_index = ii, peak_pseudotime = windows$center[ii],
          interval_start = windows$left[ii],
          interval_end = windows$right[ii],
          window_value = window_value[ii], baseline = baseline,
          effect_size = effect,
          raw_effect_size = raw_window_value[ii] - raw_baseline,
          n_cells = windows$n_cells[ii],
          n_samples = windows$n_samples[ii],
          region_exceedances = exceedances,
          permutations_used = permutations,
          p_empirical = (exceedances + 1) / (permutations + 1),
          q_adjusted = NA_real_, p_dynamic = p_dynamic,
          q_dynamic = NA_real_, dynamic_significant = FALSE,
          region_supported = FALSE, merged_region_id = "",
          permutation_seed = effective_seed,
          permutation_mode = actual_permutation_mode,
          stringsAsFactors = FALSE)
      }
    }
    analysis_data[[signal_key]] <- list(pseudotime = time,
      samples = samples, range = windows$range)
  }

  .trajectory_assert(analyzable > 0L,
    "trajectory_signal_regions: no signal had enough window support.")
  curve_table <- do.call(rbind, curve_rows)
  rownames(curve_table) <- NULL
  dynamic_table <- do.call(rbind, dynamic_rows)
  rownames(dynamic_table) <- NULL
  valid_dynamic <- dynamic_table$status == "ok" &
    is.finite(dynamic_table$p_dynamic)
  dynamic_table$q_dynamic[valid_dynamic] <- stats::p.adjust(
    dynamic_table$p_dynamic[valid_dynamic], method = "BH")
  dynamic_table$dynamic_significant[valid_dynamic] <-
    dynamic_table$q_dynamic[valid_dynamic] <= dynamic_q_cutoff

  candidate_table <- if (candidate_n) do.call(rbind, candidate_rows) else
    .trajectory_regions_empty_candidates()
  region_table <- .trajectory_regions_empty_regions()
  if (nrow(candidate_table)) {
    candidate_table$q_adjusted <- stats::p.adjust(
      candidate_table$p_empirical, method = "BH")
    dynamic_match <- match(candidate_table$signal_key,
      dynamic_table$signal_key)
    candidate_table$q_dynamic <- dynamic_table$q_dynamic[dynamic_match]
    candidate_table$dynamic_significant <-
      dynamic_table$dynamic_significant[dynamic_match]
    candidate_table$region_supported <-
      candidate_table$dynamic_significant &
      candidate_table$q_adjusted <= region_q_cutoff

    supported <- which(candidate_table$region_supported)
    region_rows <- list()
    region_n <- 0L
    if (length(supported)) {
      groups <- split(supported, paste(candidate_table$signal_key[supported],
        candidate_table$direction[supported], sep = "\r"))
      for (group in groups) {
        group <- group[order(candidate_table$interval_start[group],
          candidate_table$interval_end[group])]
        datum <- analysis_data[[candidate_table$signal_key[group[1L]]]]
        gap <- merge_gap_fraction * datum$range
        clusters <- list(group[1L])
        for (index in group[-1L]) {
          current <- clusters[[length(clusters)]]
          current_right <- max(candidate_table$interval_end[current])
          if (candidate_table$interval_start[index] - current_right < gap) {
            clusters[[length(clusters)]] <- c(current, index)
          } else {
            clusters[[length(clusters) + 1L]] <- index
          }
        }
        for (cluster in clusters) {
          region_n <- region_n + 1L
          region_id <- sprintf("region_%06d", region_n)
          primary <- cluster[which.max(abs(
            candidate_table$effect_size[cluster]))]
          interval_start <- min(candidate_table$interval_start[cluster])
          interval_end <- max(candidate_table$interval_end[cluster])
          cells_in_region <- datum$pseudotime >= interval_start &
            datum$pseudotime <= interval_end
          region_rows[[region_n]] <- data.frame(region_id = region_id,
            signal_key = candidate_table$signal_key[primary],
            feature_set = candidate_table$feature_set[primary],
            lineage_id = candidate_table$lineage_id[primary],
            direction = candidate_table$direction[primary],
            interval_start = interval_start, interval_end = interval_end,
            peak_pseudotime = candidate_table$peak_pseudotime[primary],
            effect_size = candidate_table$effect_size[primary],
            raw_effect_size = candidate_table$raw_effect_size[primary],
            p_empirical = candidate_table$p_empirical[primary],
            q_adjusted = candidate_table$q_adjusted[primary],
            p_dynamic = candidate_table$p_dynamic[primary],
            q_dynamic = candidate_table$q_dynamic[primary],
            n_cells = sum(cells_in_region),
            n_samples = if (is.null(datum$samples)) NA_integer_ else
              length(unique(datum$samples[cells_in_region])),
            peak_n_cells = candidate_table$n_cells[primary],
            peak_n_samples = candidate_table$n_samples[primary],
            n_constituent_windows = length(cluster),
            constituent_candidate_ids = paste(
              candidate_table$candidate_id[cluster], collapse = ";"),
            permutations_used = candidate_table$permutations_used[primary],
            permutation_seed = candidate_table$permutation_seed[primary],
            permutation_mode = candidate_table$permutation_mode[primary],
            stringsAsFactors = FALSE)
          candidate_table$merged_region_id[cluster] <- region_id
        }
      }
      region_table <- do.call(rbind, region_rows)
      rownames(region_table) <- NULL
    }

    curve_key <- paste(curve_table$signal_key, curve_table$window_index,
      curve_table$local_direction, sep = "\r")
    candidate_key <- paste(candidate_table$signal_key,
      candidate_table$window_index, candidate_table$direction, sep = "\r")
    matched <- match(curve_key, candidate_key)
    has_candidate <- !is.na(matched)
    curve_table$p_empirical[has_candidate] <-
      candidate_table$p_empirical[matched[has_candidate]]
    curve_table$q_adjusted[has_candidate] <-
      candidate_table$q_adjusted[matched[has_candidate]]
    curve_table$region_supported[has_candidate] <-
      candidate_table$region_supported[matched[has_candidate]]
    curve_table$region_id[has_candidate] <-
      candidate_table$merged_region_id[matched[has_candidate]]
  }

  adjustment_table <- do.call(rbind, adjustment_rows)
  rownames(adjustment_table) <- NULL
  fit_summary <- data.frame(
    lineage_id = pseudotime_key, pseudotime_key = pseudotime_key,
    signal_artifact_id = as.character(signal_artifact_id),
    pseudotime_artifact_id = as.character(pseudotime_artifact_id),
    covariates_artifact_id = as.character(covariates_artifact_id),
    depth_artifact_id = as.character(depth_artifact_id),
    sample_id_artifact_id = as.character(sample_id_artifact_id),
    membership_artifact_id = as.character(membership_artifact_id),
    n_effective_cells = length(cells),
    n_pseudotime_finite = nrow(pt_finite),
    n_pseudotime_missing = length(cells) - nrow(pt_finite),
    n_signal_keys_requested = length(signal_keys),
    n_signal_keys_analyzed = analyzable,
    n_candidate_regions = nrow(candidate_table),
    n_supported_regions = nrow(region_table),
    window_width_fraction = width_fraction,
    window_step_fraction = step_fraction,
    min_cells_per_window = min_cells,
    min_samples_per_window = min_samples,
    summary_fun = summary_fun, n_permutations = permutations,
    seed = seed, permutation_mode_requested = permutation_mode,
    permutation_mode_used = actual_permutation_mode,
    sample_blocks_available = !is.null(sample_values),
    covariate_adjustment = covariate_adjustment,
    covariate_spline_df = spline_df,
    adjust_depth = adjust_depth,
    depth_field = depth_field,
    depth_adjustment_status = if (!adjust_depth) "disabled" else
      if (already_depth_adjusted) "already_adjusted" else "applied",
    depth_spline_df = if (adjust_depth) 4L else NA_integer_,
    covariate_keys = paste(covariate_keys %||% character(),
      collapse = ";"), pseudotime_direction = pseudotime_direction,
    include_endpoint_extrema = include_endpoint_extrema,
    dynamic_q_cutoff = dynamic_q_cutoff,
    region_q_cutoff = region_q_cutoff,
    merge_gap_fraction = merge_gap_fraction,
    region_null = "direction_specific_permutation_maximum_deviation",
    dynamic_null = "permutation_peak_to_trough",
    stringsAsFactors = FALSE)

  outputs <- list(
    signal_table = scworkbench::wb_prim_obs_table(curve_table,
      table_type = "trajectory_signal_dynamics_curve",
      entity_kind = "trajectory_window", meta = list(
        operation = "trajectory_signal_regions",
        lineage_id = pseudotime_key,
        pseudotime_artifact_id = as.character(pseudotime_artifact_id),
        overlapping_windows = TRUE)),
    dynamic_tests = scworkbench::wb_prim_obs_table(dynamic_table,
      table_type = "trajectory_signal_regions_dynamic",
      entity_kind = "signal"),
    candidate_regions = scworkbench::wb_prim_obs_table(candidate_table,
      table_type = "trajectory_signal_regions_candidate",
      entity_kind = "signal_region_candidate"),
    regions = scworkbench::wb_prim_obs_table(region_table,
      table_type = "trajectory_signal_regions",
      entity_kind = "signal_region"),
    adjustment_summary = scworkbench::wb_prim_obs_table(adjustment_table,
      table_type = "trajectory_signal_regions_adjustment",
      entity_kind = "signal"),
    fit_summary = scworkbench::wb_prim_obs_table(fit_summary,
      table_type = "trajectory_signal_regions_summary",
      entity_kind = "lineage"))
  list(outputs = outputs, params = list(
    pseudotime_key = pseudotime_key, signal_keys = signal_keys,
    covariate_keys = covariate_keys %||% character(),
    n_permutations = permutations, seed = seed,
    permutation_mode = actual_permutation_mode,
    covariate_adjustment = covariate_adjustment,
    adjust_depth = adjust_depth, depth_field = depth_field,
    n_supported_regions = nrow(region_table)))
}

.run_trajectory_signal_regions <- function(readers, writers, params,
                                           context) {
  result <- .trajectory_s2_call(.trajectory_signal_regions,
    list(effective_cells = as.character(context$effective_cells),
      signal = readers$signal(), pseudotime = readers$pseudotime(),
      covariates = .trajectory_regions_optional_reader(readers,
        "covariates"),
      depth = .trajectory_regions_optional_reader(readers, "depth"),
      sample_id = .trajectory_regions_optional_reader(readers,
        "sample_id"),
      membership = .trajectory_regions_optional_reader(readers,
        "membership"),
      signal_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "signal_ancestry"),
      pseudotime_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "pseudotime_ancestry"),
      covariates_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "covariates_ancestry"),
      depth_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "depth_ancestry"),
      sample_id_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "sample_id_ancestry"),
      membership_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "membership_ancestry")), params)
  .trajectory_s2_publish(result, writers, c("signal_table",
    "dynamic_tests", "candidate_regions", "regions",
    "adjustment_summary", "fit_summary"))
}
