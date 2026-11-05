# Global and trajectory-window preranked GSEA over authenticated matrix and
# feature-set registry inputs.
.trajectory_gsea_empty_enrichment <- function() data.frame(
  analysis_scope = character(), namespace = character(),
  feature_set = character(), lineage_id = character(),
  ranking_id = character(), status = character(),
  window_index = integer(), window_start_fraction = numeric(),
  window_end_fraction = numeric(), window_start = numeric(),
  window_end = numeric(), n_cells = integer(), n_samples = integer(),
  raw_set_size = integer(), overlapping_set_size = integer(),
  eligible_set_size = integer(), tested_set_size = integer(),
  enrichment_score = numeric(), normalized_enrichment_score = numeric(),
  p_value = numeric(), q_within_ranking = numeric(),
  q_global_family = numeric(), log2_error = numeric(),
  enrichment_direction = character(), leading_edge_size = integer(),
  supported = logical(), stringsAsFactors = FALSE)

.trajectory_gsea_empty_intervals <- function() data.frame(
  interval_id = character(), namespace = character(),
  feature_set = character(), lineage_id = character(),
  enrichment_direction = character(), interval_start_fraction = numeric(),
  interval_end_fraction = numeric(), interval_start = numeric(),
  interval_end = numeric(), peak_ranking_id = character(),
  peak_window_index = integer(), peak_normalized_enrichment_score = numeric(),
  peak_p_value = numeric(), peak_q_within_ranking = numeric(),
  peak_q_global_family = numeric(), minimum_q_global_family = numeric(),
  n_cells = integer(), n_samples = integer(), peak_n_cells = integer(),
  peak_n_samples = integer(), n_constituent_windows = integer(),
  constituent_ranking_ids = character(),
  consensus_leading_edge_fraction = numeric(), stringsAsFactors = FALSE)

.trajectory_gsea_empty_leading <- function() data.frame(
  analysis_scope = character(), namespace = character(),
  feature_set = character(), lineage_id = character(),
  ranking_id = character(), feature_id = character(), gene_rank = integer(),
  raw_ranking_statistic = numeric(),
  standardized_ranking_statistic = numeric(),
  enrichment_direction = character(), leading_edge_order = integer(),
  leading_edge_boundary_tied = logical(), interval_id = character(),
  interval_start_fraction = numeric(), interval_end_fraction = numeric(),
  interval_start = numeric(), interval_end = numeric(),
  peak_window = logical(), consensus_leading_edge = logical(),
  stringsAsFactors = FALSE)

.trajectory_gsea_empty_ranked <- function() data.frame(
  analysis_scope = character(), namespace = character(),
  lineage_id = character(), ranking_id = character(),
  feature_id = character(), feature_index = integer(), gene_rank = integer(),
  raw_ranking_statistic = numeric(),
  standardized_ranking_statistic = numeric(), eligible = logical(),
  exclusion_reason = character(), tie_group_id = integer(),
  tie_group_size = integer(), stringsAsFactors = FALSE)

.trajectory_gsea_empty_summary <- function() data.frame(
  analysis_scope = character(), namespace = character(),
  lineage_id = character(), ranking_id = character(), status = character(),
  rank_method = character(), window_index = integer(),
  window_start_fraction = numeric(), window_end_fraction = numeric(),
  window_start = numeric(), window_end = numeric(), n_cells = integer(),
  n_samples = integer(), early_fraction = numeric(), late_fraction = numeric(),
  early_n_cells = integer(), late_n_cells = integer(),
  early_n_samples = integer(), late_n_samples = integer(),
  matrix_id = character(),
  pseudotime_artifact_id = character(), membership_artifact_id = character(),
  sample_id_artifact_id = character(), covariates_artifact_id = character(),
  pseudotime_direction = character(), permutation_mode = character(),
  n_permutations = integer(), ranking_seed = integer(),
  fgsea_seed = integer(), covariate_adjustment = character(),
  covariate_spline_df = integer(), feature_universe_size = integer(),
  eligible_feature_count = integer(), excluded_feature_count = integer(),
  tied_feature_count = integer(), tied_group_count = integer(),
  maximum_tie_size = integer(), tested_set_count = integer(),
  excluded_set_count = integer(), fgsea_version = character(),
  score_type = character(), gsea_param = numeric(), sample_size = integer(),
  n_perm_simple = integer(), eps = numeric(), min_set_size = integer(),
  max_set_size = integer(), p_adjust_method = character(),
  stringsAsFactors = FALSE)

.trajectory_gsea_continuous <- function(value, label) {
  .trajectory_kind(value, "cell_vector_family")
  .trajectory_assert(!identical(as.character(
      value$value_type %||% "continuous")[1L], "categorical"),
    paste0("trajectory_preranked_gsea: ", label,
      " must contain continuous values."))
  table <- scworkbench::wb_primitive_to_table(value,
    include_context = FALSE)
  .trajectory_assert(all(c("key", "cell_id", "value") %in% names(table)),
    paste0("trajectory_preranked_gsea: ", label, " table is malformed."))
  table <- data.frame(key = as.character(table$key),
    cell_id = as.character(table$cell_id),
    value = suppressWarnings(as.numeric(table$value)),
    stringsAsFactors = FALSE)
  pair <- paste(table$key, table$cell_id, sep = "\r")
  .trajectory_assert(!anyNA(table$key) && !anyNA(table$cell_id) &&
      all(nzchar(table$key)) && all(nzchar(table$cell_id)) &&
      !anyDuplicated(pair), paste0("trajectory_preranked_gsea: ", label,
        " has invalid or duplicate key/cell rows."))
  table
}

.trajectory_gsea_scalar <- function(value, label) {
  .trajectory_kind(value, "cell_vector")
  table <- scworkbench::wb_primitive_to_table(value,
    include_context = FALSE)
  .trajectory_assert(all(c("cell_id", "value") %in% names(table)),
    paste0("trajectory_preranked_gsea: ", label, " table is malformed."))
  table <- data.frame(cell_id = as.character(table$cell_id),
    value = as.character(table$value), stringsAsFactors = FALSE)
  .trajectory_assert(!anyNA(table$cell_id) && all(nzchar(table$cell_id)) &&
      !anyDuplicated(table$cell_id), paste0("trajectory_preranked_gsea: ",
        label, " has invalid or duplicate cell rows."))
  table
}

.trajectory_gsea_seed <- function(seed, ranking_id, stage) {
  bytes <- utf8ToInt(enc2utf8(paste(stage, ranking_id, sep = "\r")))
  offset <- if (length(bytes)) sum((as.double(bytes) *
    seq_along(bytes)) %% 1000003) else 0
  as.integer((as.double(seed) + offset) %% 2147483646 + 1)
}

.trajectory_gsea_registry <- function(snapshot, namespace, universe) {
  .trajectory_assert(is.list(snapshot) && is.list(snapshot$sets),
    "trajectory_preranked_gsea: feature-set snapshot is unavailable.")
  snapshot_universe <- as.character(unlist(snapshot$feature_universe,
    use.names = FALSE))
  .trajectory_assert(identical(snapshot_universe, universe),
    paste0("trajectory_preranked_gsea: registry and authenticated matrix ",
      "feature universes differ."))
  matches <- Filter(function(item) identical(
    as.character(item$namespace %||% "")[1L], namespace), snapshot$sets)
  .trajectory_assert(length(matches) > 0L,
    paste0("trajectory_preranked_gsea: no feature sets found for namespace='",
      namespace, "'."))
  set_names <- vapply(matches, function(item)
    as.character(item$feature_set %||% "")[1L], character(1))
  .trajectory_assert(all(nzchar(set_names)) && !anyDuplicated(set_names),
    "trajectory_preranked_gsea: feature-set names are empty or ambiguous.")
  stats::setNames(lapply(seq_along(matches), function(ii) {
    raw <- vapply(matches[[ii]]$memberships %||% list(), function(row)
      as.character(row$feature_id %||% "")[1L], character(1))
    raw <- unique(raw[!is.na(raw) & nzchar(raw)])
    list(feature_set = set_names[[ii]], raw = raw,
      overlap = raw[raw %in% universe])
  }), set_names)
}

.trajectory_gsea_adjust <- function(matrix, covariates, method, spline_df) {
  if (identical(method, "none")) return(list(values = matrix,
    design_rank = 1L, design_columns = 1L, model = "none"))
  terms <- lapply(seq_len(ncol(covariates)), function(ii) {
    value <- covariates[, ii]
    if (identical(method, "linear")) matrix(value, ncol = 1L,
      dimnames = list(NULL, colnames(covariates)[ii])) else {
      out <- splines::ns(value, df = spline_df)
      colnames(out) <- paste0(colnames(covariates)[ii], "_ns",
        seq_len(ncol(out)))
      out
    }
  })
  design <- cbind(`(Intercept)` = 1, do.call(cbind, terms))
  fit <- qr(design)
  .trajectory_assert(fit$rank == ncol(design),
    "trajectory_preranked_gsea: covariate design is rank deficient.")
  valid <- rowSums(!is.finite(matrix)) == 0L
  adjusted <- matrix(NA_real_, nrow(matrix), ncol(matrix),
    dimnames = dimnames(matrix))
  if (any(valid)) adjusted[valid, ] <- t(qr.resid(fit,
    t(matrix[valid, , drop = FALSE])))
  list(values = adjusted, design_rank = fit$rank,
    design_columns = ncol(design), model = method)
}

.trajectory_gsea_monotonic <- function(matrix, pseudotime) {
  y <- rank(pseudotime, ties.method = "average")
  y <- y - mean(y)
  y_norm <- sqrt(sum(y^2))
  ranks <- t(apply(matrix, 1L, rank, ties.method = "average"))
  ranks <- ranks - rowMeans(ranks)
  denom <- sqrt(rowSums(ranks^2)) * y_norm
  value <- as.numeric(ranks %*% y) / denom
  value[!is.finite(value)] <- NA_real_
  value
}

.trajectory_gsea_contrast <- function(matrix, pseudotime,
                                      early_fraction, late_fraction) {
  early <- pseudotime <= early_fraction
  late <- pseudotime >= 1 - late_fraction
  list(value = rowMeans(matrix[, late, drop = FALSE]) -
      rowMeans(matrix[, early, drop = FALSE]),
    early = early, late = late)
}

.trajectory_gsea_windows <- function(pseudotime, samples, width, step,
                                     min_cells, min_samples) {
  final_start <- 1 - width
  starts <- seq(0, final_start, by = step)
  if (!length(starts) || utils::tail(starts, 1L) < final_start - 1e-12)
    starts <- c(starts, final_start)
  starts <- unique(pmin(starts, final_start))
  lapply(seq_along(starts), function(ii) {
    inside <- pseudotime >= starts[ii] & pseudotime <= starts[ii] + width
    outside <- !inside
    window_samples <- if (is.null(samples)) NA_integer_ else
      length(unique(samples[inside]))
    background_samples <- if (is.null(samples)) NA_integer_ else
      length(unique(samples[outside]))
    sample_ok <- if (is.null(samples)) min_samples == 1L else
      window_samples >= min_samples && background_samples >= min_samples
    supported <- sum(inside) >= min_cells && sum(outside) >= min_cells &&
      sample_ok
    list(ranking_id = sprintf("local_w%04d", ii), window_index = ii,
      start_fraction = starts[ii], end_fraction = starts[ii] + width,
      inside = inside, outside = outside, n_cells = sum(inside),
      n_background = sum(outside), n_samples = window_samples,
      n_background_samples = background_samples, supported = supported,
      status = if (supported) "ok" else "insufficient_support")
  })
}

.trajectory_gsea_local_effects <- function(matrix, windows) {
  vapply(windows, function(window)
    rowMeans(matrix[, window$inside, drop = FALSE]) -
      rowMeans(matrix[, window$outside, drop = FALSE]),
    numeric(nrow(matrix)))
}

.trajectory_gsea_base_exclusions <- function(raw, adjusted) {
  reason <- rep("", nrow(raw))
  non_finite <- rowSums(!is.finite(raw)) > 0L
  reason[non_finite] <- "non_finite_expression"
  undetected <- !non_finite & rowSums(raw != 0) == 0L
  reason[undetected] <- "undetected"
  variance <- apply(adjusted, 1L, stats::var)
  zero_information <- reason == "" & (!is.finite(variance) | variance <= 0)
  reason[zero_information] <- "zero_information"
  reason
}

.trajectory_gsea_finalize_ranking <- function(raw_statistic, null_sum,
                                              null_sumsq, permutations,
                                              base_reason, feature_ids,
                                              scope, namespace, lineage_id,
                                              ranking_id) {
  null_mean <- null_sum / permutations
  variance <- (null_sumsq - permutations * null_mean^2) /
    (permutations - 1L)
  variance[variance < 0 & variance > -sqrt(.Machine$double.eps)] <- 0
  null_sd <- sqrt(variance)
  reason <- base_reason
  invalid_observed <- reason == "" & !is.finite(raw_statistic)
  reason[invalid_observed] <- "non_finite_ranking_statistic"
  invalid_null <- reason == "" & (!is.finite(null_sd) | null_sd <= 0)
  reason[invalid_null] <- "zero_permutation_variance"
  eligible <- reason == ""
  statistic <- rep(NA_real_, length(feature_ids))
  statistic[eligible] <- (raw_statistic[eligible] -
    null_mean[eligible]) / null_sd[eligible]
  .trajectory_assert(sum(eligible) >= 2L,
    paste0("trajectory_preranked_gsea: ranking '", ranking_id,
      "' has fewer than two eligible features."))
  order_index <- which(eligible)[order(-statistic[eligible],
    which(eligible), method = "radix")]
  gene_rank <- rep(NA_integer_, length(feature_ids))
  gene_rank[order_index] <- seq_along(order_index)
  tie_group <- tie_size <- rep(NA_integer_, length(feature_ids))
  runs <- rle(statistic[order_index])
  group <- rep(seq_along(runs$lengths), runs$lengths)
  sizes <- rep(runs$lengths, runs$lengths)
  tie_group[order_index] <- group
  tie_size[order_index] <- sizes
  tied <- runs$lengths > 1L
  table <- data.frame(analysis_scope = scope, namespace = namespace,
    lineage_id = lineage_id, ranking_id = ranking_id,
    feature_id = feature_ids, feature_index = seq_along(feature_ids),
    gene_rank = gene_rank, raw_ranking_statistic = raw_statistic,
    standardized_ranking_statistic = statistic, eligible = eligible,
    exclusion_reason = reason, tie_group_id = tie_group,
    tie_group_size = tie_size, stringsAsFactors = FALSE)
  list(table = table,
    stats = stats::setNames(statistic[order_index], feature_ids[order_index]),
    diagnostics = list(eligible = sum(eligible),
      excluded = sum(!eligible), tied_features = sum(runs$lengths[tied]),
      tied_groups = sum(tied), maximum_tie_size = if (length(runs$lengths))
        max(runs$lengths) else 0L))
}

.trajectory_gsea_fgsea <- function(ranking, sets, scope, namespace,
                                   lineage_id, window, score_type,
                                   gsea_param, sample_size, n_perm_simple,
                                   eps, min_set_size, max_set_size,
                                   p_adjust_method, seed) {
  rows <- vector("list", length(sets))
  tested <- list()
  tested_names <- character()
  for (ii in seq_along(sets)) {
    item <- sets[[ii]]
    eligible <- item$overlap[item$overlap %in% names(ranking$stats)]
    status <- if (!length(item$overlap)) "no_matrix_overlap" else
      if (!length(eligible)) "no_eligible_features" else
      if (length(eligible) < min_set_size) "below_min_set_size" else
      if (length(eligible) > max_set_size) "above_max_set_size" else "ok"
    rows[[ii]] <- data.frame(analysis_scope = scope,
      namespace = namespace, feature_set = item$feature_set,
      lineage_id = lineage_id, ranking_id = unique(ranking$table$ranking_id),
      status = status,
      window_index = as.integer(window$window_index %||% NA_integer_),
      window_start_fraction = as.numeric(window$start_fraction %||% NA_real_),
      window_end_fraction = as.numeric(window$end_fraction %||% NA_real_),
      window_start = as.numeric(window$start %||% NA_real_),
      window_end = as.numeric(window$end %||% NA_real_),
      n_cells = as.integer(window$n_cells %||% NA_integer_),
      n_samples = as.integer(window$n_samples %||% NA_integer_),
      raw_set_size = length(item$raw),
      overlapping_set_size = length(item$overlap),
      eligible_set_size = length(eligible), tested_set_size = 0L,
      enrichment_score = NA_real_, normalized_enrichment_score = NA_real_,
      p_value = NA_real_, q_within_ranking = NA_real_,
      q_global_family = NA_real_, log2_error = NA_real_,
      enrichment_direction = "", leading_edge_size = 0L,
      supported = FALSE, stringsAsFactors = FALSE)
    if (identical(status, "ok")) {
      tested[[item$feature_set]] <- eligible
      tested_names <- c(tested_names, item$feature_set)
    }
  }
  result <- do.call(rbind, rows)
  leading <- .trajectory_gsea_empty_leading()
  if (!length(tested)) return(list(enrichment = result, leading = leading))
  fgsea_result <- .trajectory_regions_with_seed(seed,
    fgsea::fgseaMultilevel(pathways = tested, stats = ranking$stats,
      sampleSize = sample_size, minSize = min_set_size,
      maxSize = max_set_size, eps = eps, scoreType = score_type,
      nproc = 1L, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE),
      gseaParam = gsea_param, nPermSimple = n_perm_simple))
  fgsea_result <- as.data.frame(fgsea_result, stringsAsFactors = FALSE)
  .trajectory_assert(setequal(as.character(fgsea_result$pathway),
      tested_names),
    "trajectory_preranked_gsea: fgsea returned an incomplete set family.")
  leading_rows <- list()
  leading_n <- 0L
  for (name in tested_names) {
    target <- match(name, result$feature_set)
    source <- match(name, as.character(fgsea_result$pathway))
    result$tested_set_size[target] <- as.integer(fgsea_result$size[source])
    result$enrichment_score[target] <- as.numeric(fgsea_result$ES[source])
    result$normalized_enrichment_score[target] <-
      as.numeric(fgsea_result$NES[source])
    result$p_value[target] <- as.numeric(fgsea_result$pval[source])
    if ("log2err" %in% names(fgsea_result)) result$log2_error[target] <-
      as.numeric(fgsea_result$log2err[source])
    direction <- if (result$normalized_enrichment_score[target] >= 0)
      "elevated" else "reduced"
    result$enrichment_direction[target] <- direction
    genes <- as.character(fgsea_result$leadingEdge[[source]])
    result$leading_edge_size[target] <- length(genes)
    if (!length(genes)) next
    ranked <- ranking$table[ranking$table$eligible, , drop = FALSE]
    ranked <- ranked[order(ranked$gene_rank), , drop = FALSE]
    ranks <- match(genes, ranked$feature_id)
    boundary <- if (identical(direction, "elevated")) max(ranks) else
      min(ranks)
    boundary_stat <- ranked$standardized_ranking_statistic[boundary]
    tie_positions <- which(ranked$standardized_ranking_statistic ==
      boundary_stat)
    boundary_tied <- if (identical(direction, "elevated"))
      any(tie_positions > boundary) else any(tie_positions < boundary)
    for (jj in seq_along(genes)) {
      row <- ranked[match(genes[jj], ranked$feature_id), , drop = FALSE]
      leading_n <- leading_n + 1L
      leading_rows[[leading_n]] <- data.frame(analysis_scope = scope,
        namespace = namespace, feature_set = name,
        lineage_id = lineage_id, ranking_id = row$ranking_id,
        feature_id = genes[jj], gene_rank = row$gene_rank,
        raw_ranking_statistic = row$raw_ranking_statistic,
        standardized_ranking_statistic =
          row$standardized_ranking_statistic,
        enrichment_direction = direction, leading_edge_order = jj,
        leading_edge_boundary_tied = boundary_tied, interval_id = "",
        interval_start_fraction = NA_real_, interval_end_fraction = NA_real_,
        interval_start = NA_real_, interval_end = NA_real_,
        peak_window = FALSE, consensus_leading_edge = FALSE,
        stringsAsFactors = FALSE)
    }
  }
  tested_index <- result$status == "ok" & is.finite(result$p_value)
  result$q_within_ranking[tested_index] <- stats::p.adjust(
    result$p_value[tested_index], method = p_adjust_method)
  if (leading_n) leading <- do.call(rbind, leading_rows)
  list(enrichment = result, leading = leading)
}

.trajectory_gsea_merge_intervals <- function(enrichment, leading, windows,
                                             pseudotime_min,
                                             pseudotime_range, samples,
                                             merge_gap_fraction,
                                             consensus_fraction) {
  intervals <- .trajectory_gsea_empty_intervals()
  if (!nrow(enrichment) || !any(enrichment$supported))
    return(list(intervals = intervals, leading = leading))
  supported <- which(enrichment$supported)
  groups <- split(supported, paste(enrichment$feature_set[supported],
    enrichment$enrichment_direction[supported], sep = "\r"))
  interval_rows <- list()
  interval_n <- 0L
  for (group in groups) {
    group <- group[order(enrichment$window_start_fraction[group],
      enrichment$window_end_fraction[group])]
    clusters <- list(group[1L])
    for (index in group[-1L]) {
      current <- clusters[[length(clusters)]]
      gap <- enrichment$window_start_fraction[index] -
        max(enrichment$window_end_fraction[current])
      if (gap <= merge_gap_fraction) clusters[[length(clusters)]] <-
        c(current, index) else clusters[[length(clusters) + 1L]] <- index
    }
    for (cluster in clusters) {
      interval_n <- interval_n + 1L
      interval_id <- sprintf("enrichment_interval_%06d", interval_n)
      peak <- cluster[which.max(abs(
        enrichment$normalized_enrichment_score[cluster]))]
      start_fraction <- min(enrichment$window_start_fraction[cluster])
      end_fraction <- max(enrichment$window_end_fraction[cluster])
      start <- pseudotime_min + start_fraction * pseudotime_range
      end <- pseudotime_min + end_fraction * pseudotime_range
      ids <- enrichment$ranking_id[cluster]
      cell_index <- which(vapply(seq_along(windows), function(ii)
        windows[[ii]]$ranking_id %in% ids, logical(1)))
      included <- Reduce(`|`, lapply(windows[cell_index], `[[`, "inside"))
      interval_rows[[interval_n]] <- data.frame(interval_id = interval_id,
        namespace = enrichment$namespace[peak],
        feature_set = enrichment$feature_set[peak],
        lineage_id = enrichment$lineage_id[peak],
        enrichment_direction = enrichment$enrichment_direction[peak],
        interval_start_fraction = start_fraction,
        interval_end_fraction = end_fraction, interval_start = start,
        interval_end = end, peak_ranking_id = enrichment$ranking_id[peak],
        peak_window_index = enrichment$window_index[peak],
        peak_normalized_enrichment_score =
          enrichment$normalized_enrichment_score[peak],
        peak_p_value = enrichment$p_value[peak],
        peak_q_within_ranking = enrichment$q_within_ranking[peak],
        peak_q_global_family = enrichment$q_global_family[peak],
        minimum_q_global_family = min(
          enrichment$q_global_family[cluster]),
        n_cells = sum(included),
        n_samples = if (is.null(samples)) NA_integer_ else
          length(unique(samples[included])),
        peak_n_cells = enrichment$n_cells[peak],
        peak_n_samples = enrichment$n_samples[peak],
        n_constituent_windows = length(cluster),
        constituent_ranking_ids = paste(ids, collapse = ";"),
        consensus_leading_edge_fraction = consensus_fraction,
        stringsAsFactors = FALSE)
      lead_index <- which(leading$analysis_scope == "local" &
        leading$feature_set == enrichment$feature_set[peak] &
        leading$ranking_id %in% ids)
      if (length(lead_index)) {
        leading$interval_id[lead_index] <- interval_id
        leading$interval_start_fraction[lead_index] <- start_fraction
        leading$interval_end_fraction[lead_index] <- end_fraction
        leading$interval_start[lead_index] <- start
        leading$interval_end[lead_index] <- end
        leading$peak_window[lead_index] <-
          leading$ranking_id[lead_index] == enrichment$ranking_id[peak]
        appearances <- tapply(leading$ranking_id[lead_index],
          leading$feature_id[lead_index], function(value)
            length(unique(value)))
        consensus <- names(appearances)[appearances / length(ids) >=
          consensus_fraction]
        leading$consensus_leading_edge[lead_index] <-
          leading$feature_id[lead_index] %in% consensus
      }
    }
  }
  intervals <- do.call(rbind, interval_rows)
  rownames(intervals) <- NULL
  list(intervals = intervals, leading = leading)
}

.trajectory_preranked_gsea <- function(matrix, registry, effective_cells,
                                       pseudotime, namespace,
                                       pseudotime_key, membership = NULL,
                                       sample_id = NULL, covariates = NULL,
                                       analysis_scope = "both",
                                       pseudotime_direction = "as_is",
                                       global_rank_method =
                                         "monotonic_association",
                                       early_fraction = 0.2,
                                       late_fraction = 0.2,
                                       window_width_fraction = 0.15,
                                       window_step_fraction = 0.05,
                                       min_cells_per_window = 20L,
                                       min_samples_per_window = 1L,
                                       permutation_mode = "auto",
                                       n_permutations = 1000L, seed = 1L,
                                       covariate_adjustment = "none",
                                       covariate_spline_df = 4L,
                                       min_set_size = 5L,
                                       max_set_size = 5000L,
                                       score_type = "std", gsea_param = 1,
                                       sample_size = 101L,
                                       n_perm_simple = 1000L, eps = 1e-10,
                                       p_adjust_method = "BH",
                                       global_q_cutoff = 0.05,
                                       local_q_cutoff = 0.05,
                                       merge_gap_fraction = 0.05,
                                       consensus_leading_edge_fraction = 0.5,
                                       matrix_id = "",
                                       pseudotime_artifact_id = "",
                                       membership_artifact_id = "",
                                       sample_id_artifact_id = "",
                                       covariates_artifact_id = "") {
  analysis_scope <- as.character(analysis_scope)[1L]
  pseudotime_direction <- as.character(pseudotime_direction)[1L]
  global_rank_method <- as.character(global_rank_method)[1L]
  permutation_mode <- as.character(permutation_mode)[1L]
  covariate_adjustment <- as.character(covariate_adjustment)[1L]
  score_type <- as.character(score_type)[1L]
  namespace <- trimws(as.character(namespace)[1L])
  pseudotime_key <- trimws(as.character(pseudotime_key)[1L])
  .trajectory_assert(analysis_scope %in% c("global", "local", "both") &&
      pseudotime_direction %in% c("as_is", "invert") &&
      global_rank_method %in% c("monotonic_association",
        "early_late_contrast") &&
      permutation_mode %in% c("auto", "blocked", "unblocked") &&
      covariate_adjustment %in% c("none", "linear", "natural_spline") &&
      score_type %in% c("std", "pos", "neg") && nzchar(namespace) &&
      nzchar(pseudotime_key),
    "trajectory_preranked_gsea: invalid categorical parameter.")
  early_fraction <- as.numeric(early_fraction)[1L]
  late_fraction <- as.numeric(late_fraction)[1L]
  width <- as.numeric(window_width_fraction)[1L]
  step <- as.numeric(window_step_fraction)[1L]
  min_cells <- as.integer(min_cells_per_window)[1L]
  min_samples <- as.integer(min_samples_per_window)[1L]
  permutations <- as.integer(n_permutations)[1L]
  seed <- as.integer(seed)[1L]
  spline_df <- as.integer(covariate_spline_df)[1L]
  min_size <- as.integer(min_set_size)[1L]
  max_size <- as.integer(max_set_size)[1L]
  sample_size <- as.integer(sample_size)[1L]
  n_perm_simple <- as.integer(n_perm_simple)[1L]
  .trajectory_assert(all(is.finite(c(early_fraction, late_fraction, width,
      step, min_cells, min_samples, permutations, seed, spline_df,
      min_size, max_size, gsea_param, sample_size, n_perm_simple, eps,
      global_q_cutoff, local_q_cutoff, merge_gap_fraction,
      consensus_leading_edge_fraction))) && early_fraction > 0 &&
      late_fraction > 0 && early_fraction + late_fraction < 1 &&
      width > 0 && width < 1 && step > 0 && step <= 1 && min_cells >= 2L &&
      min_samples >= 1L && permutations >= 2L && seed >= 1L &&
      spline_df >= 2L && min_size >= 1L && max_size >= min_size &&
      gsea_param > 0 && sample_size >= 3L && n_perm_simple >= 1L && eps >= 0 &&
      global_q_cutoff > 0 && global_q_cutoff <= 1 && local_q_cutoff > 0 &&
      local_q_cutoff <= 1 && merge_gap_fraction >= 0 &&
      merge_gap_fraction <= 1 && consensus_leading_edge_fraction > 0 &&
      consensus_leading_edge_fraction <= 1,
    "trajectory_preranked_gsea: invalid numeric parameter.")
  .trajectory_assert(p_adjust_method %in% stats::p.adjust.methods,
    "trajectory_preranked_gsea: invalid p_adjust_method.")
  .trajectory_assert(requireNamespace("fgsea", quietly = TRUE) &&
      identical(as.character(utils::packageVersion("fgsea")), "1.34.2"),
    "trajectory_preranked_gsea: fgsea 1.34.2 is required.")

  matrix <- as.matrix(matrix)
  storage.mode(matrix) <- "double"
  feature_ids <- rownames(matrix)
  matrix_cells <- colnames(matrix)
  .trajectory_assert(length(feature_ids) > 1L && !anyNA(feature_ids) &&
      all(nzchar(feature_ids)) && !anyDuplicated(feature_ids),
    "trajectory_preranked_gsea: matrix feature IDs must be unique.")
  .trajectory_assert(length(matrix_cells) > 0L && !anyDuplicated(matrix_cells),
    "trajectory_preranked_gsea: matrix cell axis is invalid.")
  cells <- as.character(effective_cells)
  .trajectory_assert(!anyDuplicated(cells) && all(cells %in% matrix_cells),
    "trajectory_preranked_gsea: effective cell axis is invalid.")
  if (!is.null(membership)) {
    member <- .trajectory_gsea_scalar(membership, "membership")
    member_cells <- member$cell_id[!is.na(member$value) & nzchar(member$value)]
    cells <- cells[cells %in% member_cells]
  }
  pt_table <- .trajectory_gsea_continuous(pseudotime, "pseudotime")
  .trajectory_assert(pseudotime_key %in% pt_table$key,
    "trajectory_preranked_gsea: pseudotime_key was not found.")
  pt_rows <- pt_table[pt_table$key == pseudotime_key, , drop = FALSE]
  pt <- stats::setNames(pt_rows$value, pt_rows$cell_id)
  cells <- cells[cells %in% names(pt) & is.finite(pt[cells])]
  .trajectory_assert(length(cells) >= 3L,
    "trajectory_preranked_gsea: insufficient finite pseudotime support.")
  pt <- as.numeric(pt[cells])
  if (identical(pseudotime_direction, "invert")) pt <- max(pt) - pt
  pt_min <- min(pt); pt_range <- max(pt) - pt_min
  .trajectory_assert(is.finite(pt_range) && pt_range > 0,
    "trajectory_preranked_gsea: pseudotime has zero range.")
  pt_normalized <- (pt - pt_min) / pt_range

  sample_values <- NULL
  if (!is.null(sample_id)) {
    sample_table <- .trajectory_gsea_scalar(sample_id, "sample_id")
    sample_map <- stats::setNames(sample_table$value, sample_table$cell_id)
    .trajectory_assert(all(cells %in% names(sample_map)) &&
        !anyNA(sample_map[cells]) && all(nzchar(sample_map[cells])),
      "trajectory_preranked_gsea: sample_id is incomplete.")
    sample_values <- as.character(sample_map[cells])
  }
  actual_permutation_mode <- if (identical(permutation_mode, "auto")) {
    if (is.null(sample_values)) "unblocked" else "blocked"
  } else permutation_mode
  .trajectory_assert(!identical(actual_permutation_mode, "blocked") ||
      !is.null(sample_values),
    "trajectory_preranked_gsea: blocked permutations require sample_id.")
  .trajectory_assert(!is.null(sample_values) || min_samples == 1L,
    paste0("trajectory_preranked_gsea: min_samples_per_window > 1 requires ",
      "sample_id."))
  blocks <- if (identical(actual_permutation_mode, "blocked"))
    sample_values else NULL
  if (!is.null(blocks)) .trajectory_assert(any(table(blocks) >= 2L),
    "trajectory_preranked_gsea: no sample block can be permuted.")

  covariate_matrix <- NULL
  if (!is.null(covariates)) {
    covariate_table <- .trajectory_gsea_continuous(covariates,
      "covariates")
    covariate_keys <- unique(covariate_table$key)
    covariate_matrix <- vapply(covariate_keys, function(key) {
      rows <- covariate_table[covariate_table$key == key, , drop = FALSE]
      map <- stats::setNames(rows$value, rows$cell_id)
      as.numeric(map[cells])
    }, numeric(length(cells)))
    if (is.null(dim(covariate_matrix))) covariate_matrix <- matrix(
      covariate_matrix, ncol = 1L)
    colnames(covariate_matrix) <- covariate_keys
    complete <- apply(covariate_matrix, 1L, function(value)
      all(is.finite(value)))
    cells <- cells[complete]
    pt <- pt[complete]; pt_normalized <- pt_normalized[complete]
    matrix <- matrix[, match(cells, matrix_cells), drop = FALSE]
    covariate_matrix <- covariate_matrix[complete, , drop = FALSE]
    if (!is.null(sample_values)) sample_values <- sample_values[complete]
    blocks <- if (identical(actual_permutation_mode, "blocked"))
      sample_values else NULL
  } else matrix <- matrix[, match(cells, matrix_cells), drop = FALSE]
  .trajectory_assert((identical(covariate_adjustment, "none") &&
      is.null(covariate_matrix)) || (!identical(covariate_adjustment,
      "none") && !is.null(covariate_matrix)),
    paste0("trajectory_preranked_gsea: covariates and ",
      "covariate_adjustment must be supplied together."))
  .trajectory_assert(length(cells) >= 3L,
    "trajectory_preranked_gsea: insufficient complete cells.")
  pt_min <- min(pt)
  pt_range <- max(pt) - pt_min
  .trajectory_assert(is.finite(pt_range) && pt_range > 0,
    "trajectory_preranked_gsea: complete pseudotime has zero range.")
  pt_normalized <- (pt - pt_min) / pt_range
  if (!is.null(blocks)) .trajectory_assert(any(table(blocks) >= 2L),
    "trajectory_preranked_gsea: no complete sample block can be permuted.")
  adjusted <- .trajectory_gsea_adjust(matrix, covariate_matrix,
    covariate_adjustment, spline_df)
  base_reason <- .trajectory_gsea_base_exclusions(matrix, adjusted$values)
  sets <- .trajectory_gsea_registry(registry, namespace, feature_ids)

  want_global <- analysis_scope %in% c("global", "both")
  want_local <- analysis_scope %in% c("local", "both")
  rankings <- list()
  ranking_windows <- list()
  endpoint_support <- list()
  ranking_seeds <- integer()
  fgsea_seeds <- integer()

  if (want_global) {
    ranking_id <- paste0("global_", global_rank_method)
    if (identical(global_rank_method, "early_late_contrast")) {
      observed <- .trajectory_gsea_contrast(adjusted$values,
        pt_normalized, early_fraction, late_fraction)
      .trajectory_assert(sum(observed$early) >= min_cells &&
          sum(observed$late) >= min_cells,
        "trajectory_preranked_gsea: early/late contrast lacks cell support.")
      if (!is.null(sample_values)) .trajectory_assert(
        length(unique(sample_values[observed$early])) >= min_samples &&
          length(unique(sample_values[observed$late])) >= min_samples,
        "trajectory_preranked_gsea: early/late contrast lacks sample support.")
      raw_statistic <- observed$value
      endpoint_support[[ranking_id]] <- list(
        early_fraction = early_fraction, late_fraction = late_fraction,
        early_n_cells = sum(observed$early), late_n_cells = sum(observed$late),
        early_n_samples = if (is.null(sample_values)) NA_integer_ else
          length(unique(sample_values[observed$early])),
        late_n_samples = if (is.null(sample_values)) NA_integer_ else
          length(unique(sample_values[observed$late])))
    } else raw_statistic <- .trajectory_gsea_monotonic(adjusted$values,
      pt_normalized)
    if (is.null(endpoint_support[[ranking_id]])) endpoint_support[[ranking_id]] <-
      list(early_fraction = NA_real_, late_fraction = NA_real_,
        early_n_cells = NA_integer_, late_n_cells = NA_integer_,
        early_n_samples = NA_integer_, late_n_samples = NA_integer_)
    null_sum <- null_sumsq <- numeric(nrow(matrix))
    ranking_seed <- .trajectory_gsea_seed(seed, ranking_id, "ranking")
    null <- .trajectory_regions_with_seed(ranking_seed, {
      for (permutation in seq_len(permutations)) {
        permuted <- .trajectory_regions_permute(pt_normalized, blocks)
        value <- if (identical(global_rank_method,
          "early_late_contrast")) .trajectory_gsea_contrast(
            adjusted$values, permuted, early_fraction,
            late_fraction)$value else .trajectory_gsea_monotonic(
              adjusted$values, permuted)
        value[!is.finite(value)] <- 0
        null_sum <- null_sum + value
        null_sumsq <- null_sumsq + value^2
      }
      list(sum = null_sum, sumsq = null_sumsq)
    })
    rankings[[ranking_id]] <- .trajectory_gsea_finalize_ranking(
      raw_statistic, null$sum, null$sumsq, permutations, base_reason,
      feature_ids, "global", namespace, pseudotime_key, ranking_id)
    ranking_windows[[ranking_id]] <- list(window_index = NA_integer_,
      start_fraction = NA_real_, end_fraction = NA_real_,
      start = NA_real_, end = NA_real_, n_cells = length(cells),
      n_samples = if (is.null(sample_values)) NA_integer_ else
        length(unique(sample_values)), supported = TRUE, status = "ok")
    ranking_seeds[ranking_id] <- ranking_seed
    fgsea_seeds[ranking_id] <- .trajectory_gsea_seed(seed, ranking_id,
      "fgsea")
  }

  all_windows <- .trajectory_gsea_windows(pt_normalized, sample_values,
    width, step, min_cells, min_samples)
  supported_windows <- all_windows[vapply(all_windows, `[[`, logical(1),
    "supported")]
  if (want_local) {
    .trajectory_assert(length(supported_windows) > 0L,
      "trajectory_preranked_gsea: no local window passes support thresholds.")
    observed_local <- .trajectory_gsea_local_effects(adjusted$values,
      supported_windows)
    null_sum <- null_sumsq <- matrix(0, nrow(matrix),
      length(supported_windows))
    local_seed <- .trajectory_gsea_seed(seed, "local_grid", "ranking")
    null <- .trajectory_regions_with_seed(local_seed, {
      for (permutation in seq_len(permutations)) {
        permuted <- .trajectory_regions_permute(pt_normalized, blocks)
        perm_windows <- .trajectory_gsea_windows(permuted, sample_values,
          width, step, min_cells, min_samples)
        wanted <- vapply(supported_windows, `[[`, character(1),
          "ranking_id")
        available <- vapply(perm_windows, `[[`, character(1), "ranking_id")
        perm_windows <- perm_windows[match(wanted, available)]
        value <- .trajectory_gsea_local_effects(adjusted$values,
          perm_windows)
        value[!is.finite(value)] <- 0
        null_sum <- null_sum + value
        null_sumsq <- null_sumsq + value^2
      }
      list(sum = null_sum, sumsq = null_sumsq)
    })
    for (ii in seq_along(supported_windows)) {
      window <- supported_windows[[ii]]
      ranking_id <- window$ranking_id
      rankings[[ranking_id]] <- .trajectory_gsea_finalize_ranking(
        observed_local[, ii], null$sum[, ii], null$sumsq[, ii],
        permutations, base_reason, feature_ids, "local", namespace,
        pseudotime_key, ranking_id)
      window$start <- pt_min + window$start_fraction * pt_range
      window$end <- pt_min + window$end_fraction * pt_range
      ranking_windows[[ranking_id]] <- window
      endpoint_support[[ranking_id]] <- list(early_fraction = NA_real_,
        late_fraction = NA_real_, early_n_cells = NA_integer_,
        late_n_cells = NA_integer_, early_n_samples = NA_integer_,
        late_n_samples = NA_integer_)
      ranking_seeds[ranking_id] <- local_seed
      fgsea_seeds[ranking_id] <- .trajectory_gsea_seed(seed, ranking_id,
        "fgsea")
    }
  }

  enrichment_results <- leading_results <- list()
  for (ranking_id in names(rankings)) {
    scope <- unique(rankings[[ranking_id]]$table$analysis_scope)
    value <- .trajectory_gsea_fgsea(rankings[[ranking_id]], sets, scope,
      namespace, pseudotime_key, ranking_windows[[ranking_id]], score_type,
      gsea_param, sample_size, n_perm_simple, eps, min_size, max_size,
      p_adjust_method, fgsea_seeds[[ranking_id]])
    enrichment_results[[ranking_id]] <- value$enrichment
    leading_results[[ranking_id]] <- value$leading
  }
  enrichment <- if (length(enrichment_results)) do.call(rbind,
    enrichment_results) else .trajectory_gsea_empty_enrichment()
  rownames(enrichment) <- NULL
  leading <- if (length(leading_results)) do.call(rbind,
    leading_results) else .trajectory_gsea_empty_leading()
  rownames(leading) <- NULL
  global_enrichment <- enrichment[enrichment$analysis_scope == "global",
    , drop = FALSE]
  window_enrichment <- enrichment[enrichment$analysis_scope == "local",
    , drop = FALSE]
  if (nrow(global_enrichment)) {
    valid <- global_enrichment$status == "ok" &
      is.finite(global_enrichment$p_value)
    global_enrichment$q_global_family[valid] <- stats::p.adjust(
      global_enrichment$p_value[valid], method = p_adjust_method)
    global_enrichment$supported[valid] <-
      global_enrichment$q_global_family[valid] <= global_q_cutoff
  }
  if (nrow(window_enrichment)) {
    valid <- window_enrichment$status == "ok" &
      is.finite(window_enrichment$p_value)
    window_enrichment$q_global_family[valid] <- stats::p.adjust(
      window_enrichment$p_value[valid], method = p_adjust_method)
    window_enrichment$supported[valid] <-
      window_enrichment$q_global_family[valid] <= local_q_cutoff
  }
  merged <- .trajectory_gsea_merge_intervals(window_enrichment, leading,
    supported_windows, pt_min, pt_range, sample_values,
    merge_gap_fraction, consensus_leading_edge_fraction)
  leading <- merged$leading
  ranked_features <- if (length(rankings)) do.call(rbind,
    lapply(rankings, `[[`, "table")) else .trajectory_gsea_empty_ranked()
  rownames(ranked_features) <- NULL

  summary_rows <- list()
  summary_n <- 0L
  for (scope in c("global", "local")) {
    requested <- if (identical(scope, "global")) want_global else want_local
    scope_ids <- names(rankings)[vapply(rankings, function(value)
      identical(unique(value$table$analysis_scope), scope), logical(1))]
    if (!requested) {
      summary_n <- summary_n + 1L
      row <- .trajectory_gsea_empty_summary()[NA_integer_, , drop = FALSE]
      row[1L, ] <- NA
      row$analysis_scope <- scope; row$namespace <- namespace
      row$lineage_id <- pseudotime_key; row$ranking_id <- ""
      row$status <- "not_requested"; row$rank_method <- ""
      summary_rows[[summary_n]] <- row
      next
    }
    for (ranking_id in scope_ids) {
      ranking <- rankings[[ranking_id]]
      window <- ranking_windows[[ranking_id]]
      endpoints <- endpoint_support[[ranking_id]]
      rows <- if (identical(scope, "global")) global_enrichment else
        window_enrichment[window_enrichment$ranking_id == ranking_id,
          , drop = FALSE]
      summary_n <- summary_n + 1L
      summary_rows[[summary_n]] <- data.frame(analysis_scope = scope,
        namespace = namespace, lineage_id = pseudotime_key,
        ranking_id = ranking_id, status = "ok",
        rank_method = if (identical(scope, "global")) global_rank_method else
          "local_vs_trajectory_background",
        window_index = as.integer(window$window_index %||% NA_integer_),
        window_start_fraction = as.numeric(
          window$start_fraction %||% NA_real_),
        window_end_fraction = as.numeric(window$end_fraction %||% NA_real_),
        window_start = as.numeric(window$start %||% NA_real_),
        window_end = as.numeric(window$end %||% NA_real_),
        n_cells = as.integer(window$n_cells),
        n_samples = as.integer(window$n_samples),
        early_fraction = as.numeric(endpoints$early_fraction),
        late_fraction = as.numeric(endpoints$late_fraction),
        early_n_cells = as.integer(endpoints$early_n_cells),
        late_n_cells = as.integer(endpoints$late_n_cells),
        early_n_samples = as.integer(endpoints$early_n_samples),
        late_n_samples = as.integer(endpoints$late_n_samples),
        matrix_id = matrix_id,
        pseudotime_artifact_id = pseudotime_artifact_id,
        membership_artifact_id = membership_artifact_id,
        sample_id_artifact_id = sample_id_artifact_id,
        covariates_artifact_id = covariates_artifact_id,
        pseudotime_direction = pseudotime_direction,
        permutation_mode = actual_permutation_mode,
        n_permutations = permutations,
        ranking_seed = as.integer(ranking_seeds[[ranking_id]]),
        fgsea_seed = as.integer(fgsea_seeds[[ranking_id]]),
        covariate_adjustment = covariate_adjustment,
        covariate_spline_df = spline_df,
        feature_universe_size = length(feature_ids),
        eligible_feature_count = ranking$diagnostics$eligible,
        excluded_feature_count = ranking$diagnostics$excluded,
        tied_feature_count = ranking$diagnostics$tied_features,
        tied_group_count = ranking$diagnostics$tied_groups,
        maximum_tie_size = ranking$diagnostics$maximum_tie_size,
        tested_set_count = sum(rows$status == "ok"),
        excluded_set_count = sum(rows$status != "ok"),
        fgsea_version = as.character(utils::packageVersion("fgsea")),
        score_type = score_type, gsea_param = gsea_param,
        sample_size = sample_size, n_perm_simple = n_perm_simple, eps = eps,
        min_set_size = min_size, max_set_size = max_size,
        p_adjust_method = p_adjust_method, stringsAsFactors = FALSE)
    }
  }
  fit_summary <- do.call(rbind, summary_rows)
  rownames(fit_summary) <- NULL
  global_status <- if (want_global) "ok" else "not_requested"
  local_status <- if (want_local) "ok" else "not_requested"
  outputs <- list(
    global_enrichment = scworkbench::wb_prim_obs_table(global_enrichment,
      table_type = "trajectory_preranked_gsea_global_enrichment",
      entity_kind = "feature_set", meta = list(status = global_status)),
    window_enrichment = scworkbench::wb_prim_obs_table(window_enrichment,
      table_type = "trajectory_preranked_gsea_window_enrichment",
      entity_kind = "feature_set_window", meta = list(status = local_status)),
    enrichment_intervals = scworkbench::wb_prim_obs_table(merged$intervals,
      table_type = "trajectory_preranked_gsea_enrichment_intervals",
      entity_kind = "feature_set_interval", meta = list(
        status = local_status)),
    leading_edge = scworkbench::wb_prim_obs_table(leading,
      table_type = "trajectory_preranked_gsea_leading_edge",
      entity_kind = "feature_set_feature", meta = list(status = "ok")),
    ranked_features = scworkbench::wb_prim_obs_table(ranked_features,
      table_type = "trajectory_preranked_gsea_ranked_features",
      entity_kind = "ranking_feature", meta = list(status = "ok")),
    fit_summary = scworkbench::wb_prim_obs_table(fit_summary,
      table_type = "trajectory_preranked_gsea_summary",
      entity_kind = "ranking"))
  list(outputs = outputs, params = list(namespace = namespace,
    pseudotime_key = pseudotime_key, analysis_scope = analysis_scope,
    pseudotime_direction = pseudotime_direction,
    global_rank_method = global_rank_method,
    permutation_mode = actual_permutation_mode,
    n_permutations = permutations, seed = seed,
    p_adjust_method = p_adjust_method))
}

.run_trajectory_preranked_gsea <- function(readers, writers, params,
                                           context) {
  registry <- readers$feature_sets()
  all_features <- as.character(
    context$typed_slice_reads$matrix_chunks$features)
  .trajectory_assert(length(all_features) > 1L &&
      !anyDuplicated(all_features),
    "trajectory_preranked_gsea: authenticated matrix feature axis is invalid.")
  cells <- as.character(context$effective_cells)
  matrix <- as.matrix(.trajectory_chunked_matrix(readers$matrix_chunks,
    all_features, all_features, cells, cells, params$chunk_size %||% 128L))
  call_params <- params
  call_params$chunk_size <- NULL
  result <- .trajectory_s2_call(.trajectory_preranked_gsea,
    list(matrix = matrix, registry = registry, effective_cells = cells,
      pseudotime = readers$pseudotime(),
      membership = .trajectory_regions_optional_reader(readers,
        "membership"),
      sample_id = .trajectory_regions_optional_reader(readers,
        "sample_id"),
      covariates = .trajectory_regions_optional_reader(readers,
        "covariates"),
      matrix_id = as.character(registry$matrix_id)[1L],
      pseudotime_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "pseudotime_ancestry"),
      membership_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "membership_ancestry"),
      sample_id_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "sample_id_ancestry"),
      covariates_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "covariates_ancestry")), call_params)
  .trajectory_s2_publish(result, writers, c("global_enrichment",
    "window_enrichment", "enrichment_intervals", "leading_edge",
    "ranked_features", "fit_summary"))
}
