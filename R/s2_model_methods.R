# Numerical S2 model consumers. Each receives a fitted lineage projection
# extracted from a G1 reader and ordinary typed inputs; no bundle is exposed.
.trajectory_s2_model_features <- function(snapshot, universe,
                                          feature_set = NULL, features = NULL,
                                          required = FALSE, prefix) {
  .trajectory_assert(is.list(snapshot) && is.list(snapshot$sets) &&
      identical(as.character(unlist(snapshot$feature_universe,
        use.names = FALSE)), as.character(universe)),
    paste0(prefix, ": authenticated feature-set axis is inconsistent."))
  set_features <- NULL
  if (!is.null(feature_set)) {
    name <- as.character(feature_set)[1L]
    .trajectory_assert(!is.na(name) && nzchar(name), paste0(prefix,
      ": feature_set must be nonempty."))
    matches <- Filter(function(item) identical(
      as.character(item$feature_set %||% "")[1L], name), snapshot$sets)
    .trajectory_assert(length(matches) == 1L, paste0(prefix,
      ": feature_set does not resolve uniquely."))
    set_features <- unique(vapply(matches[[1L]]$memberships %||% list(),
      function(item) as.character(item$feature_id)[1L], character(1)))
    set_features <- set_features[set_features %in% universe]
    .trajectory_assert(length(set_features) > 0L,
      paste0(prefix, ": feature_set has no matrix members."))
  }
  explicit <- NULL
  if (!is.null(features)) {
    explicit <- unique(as.character(unlist(features, use.names = FALSE)))
    explicit <- explicit[!is.na(explicit) & nzchar(explicit)]
    .trajectory_assert(length(explicit) > 0L &&
        !length(setdiff(explicit, universe)),
      paste0(prefix, ": features is empty or outside the matrix axis."))
  }
  selected <- if (is.null(set_features)) explicit else if (is.null(explicit))
    set_features else intersect(set_features, explicit)
  if (is.null(selected)) selected <- if (required) character() else universe
  .trajectory_assert(length(selected) > 0L, paste0(prefix,
    ": no features remain after selection."))
  selected
}

.trajectory_s2_model_matrix <- function(matrix, cells, features) {
  .trajectory_assert(is.matrix(matrix) && !is.null(rownames(matrix)) &&
      !is.null(colnames(matrix)) && !anyDuplicated(rownames(matrix)) &&
      !anyDuplicated(colnames(matrix)) &&
      !anyNA(match(features, rownames(matrix))) &&
      !anyNA(match(cells, colnames(matrix))),
    "S2 model method received an incomplete typed matrix slice.")
  matrix[features, cells, drop = FALSE]
}

.trajectory_s2_de_model <- function(matrix, source, effective_cells,
                                    matrix_id, feature_set = NULL,
                                    features = NULL, test = "spearman",
                                    lineage = NULL, min_pct_nonzero = 0,
                                    min_sd = 0, p_adjust_method = "BH") {
  .trajectory_assert(identical(as.character(test)[1L], "spearman"),
    "trajectory_de supports test='spearman' only.")
  .trajectory_assert(as.character(p_adjust_method)[1L] %in%
      stats::p.adjust.methods, "trajectory_de: invalid p_adjust_method.")
  all_axis <- source$axis[source$axis$cell_id %in% effective_cells &
    is.finite(source$axis$pseudotime), , drop = FALSE]
  n_lineages_input <- length(unique(as.character(all_axis$lineage_id)))
  lineages <- .trajectory_s2_model_lineages(source, lineage,
    effective_cells)
  if (!is.null(lineage)) {
    requested <- unique(as.character(unlist(lineage, use.names = FALSE)))
    lineages <- lineages[requested]
  }
  rows <- list(); summaries <- list(); ri <- 0L; si <- 0L
  feats <- rownames(matrix)
  for (lin in names(lineages)) {
    item <- lineages[[lin]]
    pt <- as.numeric(item$pseudotime)
    names(pt) <- as.character(item$cell_id)
    pt <- pt[is.finite(pt)]
    pt <- pt[!duplicated(names(pt))]
    if (length(pt) < 10L || !is.finite(stats::sd(pt)) ||
        stats::sd(pt) <= 0) next
    m <- .trajectory_s2_model_matrix(matrix, names(pt), feats)
    one_rows <- list(); oi <- 0L
    filtered_cells <- filtered_nonzero <- filtered_sd <- 0L
    for (k in seq_len(nrow(m))) {
      y0 <- as.numeric(m[k, ])
      ok <- is.finite(y0) & is.finite(pt)
      n_used <- sum(ok)
      n_dropped <- length(ok) - n_used
      if (n_used < 10L) { filtered_cells <- filtered_cells + 1L; next }
      y <- y0[ok]; pt_ok <- pt[ok]
      nz <- mean(y > 0)
      if (is.finite(min_pct_nonzero) && min_pct_nonzero > 0 &&
          nz < as.numeric(min_pct_nonzero)) {
        filtered_nonzero <- filtered_nonzero + 1L; next
      }
      sd_y <- stats::sd(y)
      if (!is.finite(sd_y) || sd_y < as.numeric(min_sd)) {
        filtered_sd <- filtered_sd + 1L; next
      }
      rho <- suppressWarnings(stats::cor(y, pt_ok, method = "spearman"))
      if (!is.finite(rho)) next
      pval <- if (abs(rho) >= 1) 0 else {
        tval <- rho * sqrt((n_used - 2) /
          pmax(1 - rho^2, .Machine$double.eps))
        2 * stats::pt(abs(tval), df = n_used - 2, lower.tail = FALSE)
      }
      slope <- tryCatch(as.numeric(stats::coef(
        stats::lm(y ~ pt_ok))[2]), error = function(e) NA_real_)
      oi <- oi + 1L
      one_rows[[oi]] <- data.frame(lineage_id = as.character(lin),
        feature_id = rownames(m)[k], rho = rho, abs_rho = abs(rho),
        direction = ifelse(rho > 0, "up", ifelse(rho < 0, "down", "flat")),
        slope = slope, p_val = pval, pct_nonzero = nz, sd_expr = sd_y,
        mean_expr = mean(y), n_used = as.integer(n_used),
        n_dropped_na = as.integer(n_dropped), stringsAsFactors = FALSE)
    }
    tbl <- if (oi) do.call(rbind, one_rows) else data.frame(
      lineage_id = character(), feature_id = character(), rho = numeric(),
      abs_rho = numeric(), direction = character(), slope = numeric(),
      p_val = numeric(), p_val_adj = numeric(), pct_nonzero = numeric(),
      sd_expr = numeric(), mean_expr = numeric(), n_used = integer(),
      n_dropped_na = integer(), stringsAsFactors = FALSE)
    ri <- ri + 1L; rows[[ri]] <- tbl
    si <- si + 1L
    summaries[[si]] <- data.frame(lineage_id = as.character(lin),
      n_cells_input = as.integer(length(item$cell_id)),
      n_cells_used = as.integer(length(pt)),
      n_features_input = as.integer(length(feats)),
      n_features_tested = as.integer(nrow(tbl)),
      n_filtered_min_cells = as.integer(filtered_cells),
      n_filtered_min_pct_nonzero = as.integer(filtered_nonzero),
      n_filtered_min_sd = as.integer(filtered_sd), stringsAsFactors = FALSE)
  }
  .trajectory_assert(ri > 0L,
    "trajectory_de: no lineage produced testable cells/features.")
  de <- do.call(rbind, rows); summary <- do.call(rbind, summaries)
  de$p_val_adj <- stats::ave(de$p_val, de$lineage_id,
    FUN = function(x) stats::p.adjust(x,
      method = as.character(p_adjust_method)[1L]))
  de <- de[order(de$lineage_id, de$p_val_adj, -de$abs_rho), , drop = FALSE]
  rownames(de) <- rownames(summary) <- NULL
  list(outputs = list(
    de_table = scworkbench::wb_prim_feature_table(de,
      table_type = "trajectory_de", meta = list(matrix_id = matrix_id,
        model_engine = source$engine,
        n_lineages_input = as.integer(n_lineages_input),
        n_lineages_tested = as.integer(length(unique(de$lineage_id))),
        n_features_tested = as.integer(nrow(de)),
        feature_set = as.character(feature_set %||% "")[1L],
        min_pct_nonzero = as.numeric(min_pct_nonzero),
        min_sd = as.numeric(min_sd),
        p_adjust_method = as.character(p_adjust_method)[1L],
        test = "spearman")),
    fit_summary = scworkbench::wb_prim_obs_table(summary,
      table_type = "trajectory_de_summary", entity_kind = "lineage",
      meta = list(matrix_id = matrix_id, model_engine = source$engine,
        test = "spearman"))))
}

.trajectory_s2_model_pt <- function(item, lineage_id, prefix) {
  pt <- suppressWarnings(as.numeric(item$pseudotime))
  names(pt) <- as.character(item$cell_id)
  pt <- pt[is.finite(pt)]
  .trajectory_assert(length(pt) >= 10L,
    paste0(prefix, ": lineage '", lineage_id,
      "' has <10 finite pseudotime cells."))
  .trajectory_assert(stats::sd(pt) > 0,
    paste0(prefix, ": lineage '", lineage_id,
      "' has zero pseudotime variance."))
  pt
}

.trajectory_s2_smooth_model <- function(matrix, source, effective_cells,
                                        matrix_id, lineage = NULL,
                                        n_bins = 50L, summary_fun = "mean",
                                        min_cells_per_bin = 5L) {
  summary_fun <- .trajectory_s2_summary_fun(summary_fun)
  min_cells_per_bin <- suppressWarnings(as.integer(min_cells_per_bin))
  .trajectory_assert(is.finite(min_cells_per_bin) &&
      min_cells_per_bin >= 1L,
    "trajectory_smooth: min_cells_per_bin must be >= 1.")
  lineages <- .trajectory_s2_model_lineages(source, lineage,
    effective_cells)
  feats <- rownames(matrix)
  rows <- list(); fits <- list(); ri <- 0L; si <- 0L
  for (lin in names(lineages)) {
    pt <- .trajectory_s2_model_pt(lineages[[lin]], lin,
      "trajectory_smooth")
    bb <- .trajectory_s2_bins(pt, n_bins)
    m <- .trajectory_s2_model_matrix(matrix, names(pt), feats)
    n_rows <- n_entities <- 0L
    for (k in seq_len(nrow(m))) {
      agg <- .trajectory_s2_aggregate_one(rownames(m)[k], lin, pt,
        as.numeric(m[k, ]), bb$bins, min_cells_per_bin,
        value_col = "expr", summary_fun = summary_fun)
      if (is.null(agg) || !nrow(agg)) next
      agg$feature_id <- as.character(agg$entity_id)
      agg$label <- as.character(agg$feature_id)
      agg <- agg[, c("feature_id", "lineage_id", "bin", "pseudotime",
        "pseudotime_min", "pseudotime_max", "expr", "n_cells", "label"),
        drop = FALSE]
      ri <- ri + 1L; rows[[ri]] <- agg
      n_rows <- n_rows + nrow(agg); n_entities <- n_entities + 1L
    }
    si <- si + 1L
    fits[[si]] <- data.frame(lineage_id = as.character(lin),
      n_cells_used = as.integer(length(pt)),
      n_features_input = as.integer(length(feats)),
      n_features_with_rows = as.integer(n_entities),
      n_rows = as.integer(n_rows), n_bins = as.integer(bb$n_bins),
      summary_fun = as.character(summary_fun),
      min_cells_per_bin = as.integer(min_cells_per_bin),
      matrix_id = as.character(matrix_id), stringsAsFactors = FALSE)
  }
  tbl <- if (ri) do.call(rbind, rows) else data.frame(
    feature_id = character(), lineage_id = character(), bin = integer(),
    pseudotime = numeric(), pseudotime_min = numeric(),
    pseudotime_max = numeric(), expr = numeric(), n_cells = integer(),
    label = character(), stringsAsFactors = FALSE)
  fit <- do.call(rbind, fits)
  rownames(tbl) <- rownames(fit) <- NULL
  list(outputs = list(
    trend_table = scworkbench::wb_prim_obs_table(tbl,
      table_type = "trajectory_smooth", entity_kind = "trajectory_bin",
      meta = list(matrix_id = matrix_id, summary_fun = summary_fun,
        min_cells_per_bin = as.integer(min_cells_per_bin),
        n_bins = as.integer(unique(fit$n_bins)[1L]))),
    fit_summary = scworkbench::wb_prim_obs_table(fit,
      table_type = "trajectory_smooth_summary", entity_kind = "lineage")))
}

.trajectory_s2_signature_model <- function(matrix, source,
                                           effective_cells, matrix_id,
                                           feature_set, lineage = NULL,
                                           score_method = "mean",
                                           n_bins = 50L,
                                           summary_fun = "mean",
                                           min_cells_per_bin = 5L) {
  .trajectory_assert(is.character(feature_set) &&
      length(feature_set) == 1L && !is.na(feature_set) &&
      nzchar(feature_set),
    "trajectory_signature_smooth: feature_set is required.")
  score_method <- as.character(score_method %||% "mean")[1L]
  .trajectory_assert(score_method %in% c("mean", "sum"),
    "trajectory_signature_smooth: score_method must be mean or sum.")
  summary_fun <- .trajectory_s2_summary_fun(summary_fun)
  min_cells_per_bin <- suppressWarnings(as.integer(min_cells_per_bin))
  .trajectory_assert(is.finite(min_cells_per_bin) &&
      min_cells_per_bin >= 1L,
    "trajectory_signature_smooth: min_cells_per_bin must be >= 1.")
  lineages <- .trajectory_s2_model_lineages(source, lineage,
    effective_cells)
  feats <- rownames(matrix)
  rows <- list(); fits <- list(); ri <- si <- 0L
  for (lin in names(lineages)) {
    pt <- .trajectory_s2_model_pt(lineages[[lin]], lin,
      "trajectory_signature_smooth")
    bb <- .trajectory_s2_bins(pt, n_bins)
    m <- .trajectory_s2_model_matrix(matrix, names(pt), feats)
    signal <- if (identical(score_method, "sum"))
      colSums(m, na.rm = TRUE) else colMeans(m, na.rm = TRUE)
    agg <- .trajectory_s2_aggregate_one(feature_set, lin, pt,
      as.numeric(signal), bb$bins, min_cells_per_bin,
      value_col = "value", summary_fun = summary_fun)
    .trajectory_assert(!is.null(agg) && nrow(agg) > 0L,
      paste0("trajectory_signature_smooth: no output rows for lineage '",
        lin, "'."))
    agg$feature_set <- as.character(agg$entity_id)
    agg$label <- as.character(agg$feature_set)
    agg <- agg[, c("feature_set", "lineage_id", "bin", "pseudotime",
      "pseudotime_min", "pseudotime_max", "value", "n_cells", "label"),
      drop = FALSE]
    ri <- ri + 1L; rows[[ri]] <- agg
    si <- si + 1L
    fits[[si]] <- data.frame(lineage_id = as.character(lin),
      feature_set = as.character(feature_set),
      n_cells_used = as.integer(length(pt)),
      n_features = as.integer(length(feats)),
      n_rows = as.integer(nrow(agg)), n_bins = as.integer(bb$n_bins),
      score_method = as.character(score_method),
      summary_fun = as.character(summary_fun),
      min_cells_per_bin = as.integer(min_cells_per_bin),
      matrix_id = as.character(matrix_id), stringsAsFactors = FALSE)
  }
  tbl <- do.call(rbind, rows); fit <- do.call(rbind, fits)
  rownames(tbl) <- rownames(fit) <- NULL
  list(outputs = list(
    signature_table = scworkbench::wb_prim_obs_table(tbl,
      table_type = "trajectory_signature_smooth",
      entity_kind = "trajectory_bin", meta = list(matrix_id = matrix_id,
        feature_set = feature_set, n_features = as.integer(length(feats)),
        score_method = score_method, summary_fun = summary_fun,
        min_cells_per_bin = as.integer(min_cells_per_bin),
        n_bins = as.integer(unique(fit$n_bins)[1L]))),
    fit_summary = scworkbench::wb_prim_obs_table(fit,
      table_type = "trajectory_signature_smooth_summary",
      entity_kind = "lineage")))
}

.trajectory_s2_vector_model <- function(signal, source, effective_cells,
                                        lineage = NULL, signal_keys = NULL,
                                        n_bins = 50L,
                                        summary_fun = "mean",
                                        min_cells_per_bin = 5L) {
  .trajectory_kind(signal, "cell_vector_family")
  summary_fun <- .trajectory_s2_summary_fun(summary_fun)
  min_cells_per_bin <- suppressWarnings(as.integer(min_cells_per_bin))
  .trajectory_assert(is.finite(min_cells_per_bin) &&
      min_cells_per_bin >= 1L,
    "trajectory_vector_smooth: min_cells_per_bin must be >= 1.")
  tbl <- scworkbench::wb_primitive_to_table(signal,
    include_context = FALSE)
  tbl$key <- as.character(tbl$key)
  tbl$cell_id <- as.character(tbl$cell_id)
  tbl$value <- suppressWarnings(as.numeric(tbl$value))
  tbl <- tbl[nzchar(tbl$key) & nzchar(tbl$cell_id) &
    is.finite(tbl$value), , drop = FALSE]
  .trajectory_assert(nrow(tbl) > 0L,
    "trajectory_vector_smooth: signal has no finite rows.")
  available <- unique(as.character(tbl$key))
  if (!is.null(signal_keys)) {
    requested <- unique(as.character(unlist(signal_keys,
      use.names = FALSE)))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    .trajectory_assert(length(requested) > 0L &&
        !length(setdiff(requested, available)),
      "trajectory_vector_smooth: signal_keys is empty or unknown.")
    tbl <- tbl[tbl$key %in% requested, , drop = FALSE]
  }
  keys <- unique(as.character(tbl$key))
  lineages <- .trajectory_s2_model_lineages(source, lineage,
    effective_cells)
  rows <- list(); fits <- list(); ri <- si <- 0L
  for (lin in names(lineages)) {
    pt <- .trajectory_s2_model_pt(lineages[[lin]], lin,
      "trajectory_vector_smooth")
    bb <- .trajectory_s2_bins(pt, n_bins)
    line_signal <- tbl[tbl$cell_id %in% names(pt),
      c("key", "cell_id", "value"), drop = FALSE]
    .trajectory_assert(nrow(line_signal) > 0L,
      paste0("trajectory_vector_smooth: no signal rows overlap lineage '",
        lin, "'."))
    n_rows <- n_signals <- 0L
    for (key in keys) {
      item <- line_signal[line_signal$key == key, , drop = FALSE]
      if (!nrow(item)) next
      values <- item$value[match(names(pt), item$cell_id)]
      keep <- is.finite(values) & is.finite(pt)
      if (sum(keep) < min_cells_per_bin) next
      agg <- .trajectory_s2_aggregate_one(key, lin, pt[keep],
        values[keep], bb$bins[keep], min_cells_per_bin,
        value_col = "value", summary_fun = summary_fun)
      if (is.null(agg) || !nrow(agg)) next
      agg$signal_key <- as.character(agg$entity_id)
      agg$label <- as.character(agg$signal_key)
      agg <- agg[, c("signal_key", "lineage_id", "bin", "pseudotime",
        "pseudotime_min", "pseudotime_max", "value", "n_cells", "label"),
        drop = FALSE]
      ri <- ri + 1L; rows[[ri]] <- agg
      n_rows <- n_rows + nrow(agg); n_signals <- n_signals + 1L
    }
    si <- si + 1L
    fits[[si]] <- data.frame(lineage_id = as.character(lin),
      n_cells_used = as.integer(length(pt)),
      n_signal_keys_input = as.integer(length(keys)),
      n_signal_keys_with_rows = as.integer(n_signals),
      n_rows = as.integer(n_rows), n_bins = as.integer(bb$n_bins),
      summary_fun = as.character(summary_fun),
      min_cells_per_bin = as.integer(min_cells_per_bin),
      stringsAsFactors = FALSE)
  }
  signal_table <- if (ri) do.call(rbind, rows) else data.frame(
    signal_key = character(), lineage_id = character(), bin = integer(),
    pseudotime = numeric(), pseudotime_min = numeric(),
    pseudotime_max = numeric(), value = numeric(), n_cells = integer(),
    label = character(), stringsAsFactors = FALSE)
  fit <- do.call(rbind, fits)
  rownames(signal_table) <- rownames(fit) <- NULL
  list(outputs = list(
    signal_table = scworkbench::wb_prim_obs_table(signal_table,
      table_type = "trajectory_vector_smooth",
      entity_kind = "trajectory_bin", meta = list(
        signal_kind = as.character(signal$kind),
        n_signal_keys = as.integer(length(keys)),
        summary_fun = summary_fun,
        min_cells_per_bin = as.integer(min_cells_per_bin),
        n_bins = as.integer(unique(fit$n_bins)[1L]))),
    fit_summary = scworkbench::wb_prim_obs_table(fit,
      table_type = "trajectory_vector_smooth_summary",
      entity_kind = "lineage")))
}
