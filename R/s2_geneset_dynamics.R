# Pack-side S2 method, with source selection supplied by authenticated readers.
.trajectory_s2_geneset_dynamics <- function(matrix, registry, effective_cells,
                                              matrix_id,
                                              depth = NULL,
                                              depth_artifact_id = "",
                                              pseudotime = NULL,
                                              model_source = NULL,
                                              pseudotime_key = NULL,
                                              membership = NULL,
                                              namespace = NULL,
                                              pseudotime_artifact_id = NULL,
                                              membership_artifact_id = NULL,
                                              lineage = NULL,
                                              feature_sets = NULL,
                                              score_method = "mean",
                                              adjust_depth = TRUE,
                                              depth_field = "nCount_RNA",
                                              n_bins = 50L,
                                              summary_fun = "mean",
                                              min_cells_per_bin = 5L,
                                              min_set_size = 5L,
                                              max_set_size = 5000L,
                                              top_n = 50L,
                                              on_off_min_delta = 0,
                                              monotonic_min_abs_rho = 0.5,
                                              top_driver_genes = 5L,
                                              pseudotime_direction = "as_is",
                                              start_label = NULL,
                                              end_label = NULL,
                                              min_amplitude = 0,
                                              min_bins_observed = 5L,
                                              min_features_detected = 1L,
                                              min_detected_fraction = 0,
                                              min_mean_score = 0,
                                              min_max_score = 0,
                                              min_signal_to_noise = 0,
                                              low_abundance_max_score = 0.1,
                                              smooth_method = "none",
                                              loess_span = 0.5,
                                              bootstrap_n = 0L,
                                              bootstrap_seed = 1L,
                                              ci_level = 0.95,
                                              class_min_amplitude = 0,
                                              n_permutations = 0L,
                                              permutation_seed = 1L,
                                              p_adjust_method = "BH") {
  .trajectory_assert(is.matrix(matrix) && !is.null(rownames(matrix)) &&
    !is.null(colnames(matrix)), "trajectory_geneset_dynamics: matrix is malformed.")
  .trajectory_assert(identical(colnames(matrix), as.character(effective_cells)),
    "trajectory_geneset_dynamics: matrix cell axis differs from the authenticated effective axis.")
  .trajectory_assert(xor(is.null(pseudotime), is.null(model_source)),
    "trajectory_geneset_dynamics requires exactly one pseudotime source.")
  if (is.null(model_source)) .trajectory_kind(pseudotime,
    "cell_vector_family")
  .trajectory_assert(is.character(matrix_id) && length(matrix_id) == 1 && nzchar(matrix_id), "matrix_id invalid")
  .trajectory_assert(is.character(namespace) && length(namespace) == 1 && nzchar(namespace), "namespace invalid")

  adjust_depth <- isTRUE(adjust_depth)
  depth_field <- trimws(as.character(depth_field %||% "nCount_RNA")[1L])
  .trajectory_assert(!is.na(depth_field) && nzchar(depth_field),
    "trajectory_geneset_dynamics: depth_field is required.")
  if (adjust_depth) {
    depth <- as.numeric(depth[match(as.character(effective_cells),
      names(depth))])
    names(depth) <- as.character(effective_cells)
    .trajectory_assert(length(depth) == length(effective_cells) &&
        all(is.finite(depth)) && all(depth >= 0),
      "trajectory_geneset_dynamics: depth must be complete, finite, non-negative, and cell-aligned.")
  } else depth <- stats::setNames(rep(0, length(effective_cells)),
    as.character(effective_cells))

  score_method <- as.character(score_method %||% "mean")[1]
  if (!nzchar(score_method)) score_method <- "mean"
  .trajectory_assert(score_method %in% c("mean", "sum", "rank_mean"),
             "trajectory_geneset_dynamics: score_method must be 'mean', 'sum', or 'rank_mean'.")
  .trajectory_assert(!adjust_depth || identical(score_method, "mean"),
    paste0("trajectory_geneset_dynamics: adjust_depth=TRUE requires ",
      "score_method='mean'; use adjust_depth=FALSE for sum or rank_mean."))
  summary_fun <- .trajectory_s2_summary_fun(summary_fun)
  min_cells_per_bin <- suppressWarnings(as.integer(min_cells_per_bin))
  min_set_size <- suppressWarnings(as.integer(min_set_size))
  max_set_size <- suppressWarnings(as.integer(max_set_size))
  top_n <- suppressWarnings(as.integer(top_n))
  on_off_min_delta <- suppressWarnings(as.numeric(on_off_min_delta))
  monotonic_min_abs_rho <- suppressWarnings(as.numeric(monotonic_min_abs_rho))
  top_driver_genes <- suppressWarnings(as.integer(top_driver_genes))
  pseudotime_direction <- as.character(pseudotime_direction %||% "as_is")[1]
  if (!nzchar(pseudotime_direction)) pseudotime_direction <- "as_is"
  .trajectory_assert(pseudotime_direction %in% c("as_is", "invert"),
             "trajectory_geneset_dynamics: pseudotime_direction must be 'as_is' or 'invert'.")
  start_label <- as.character(start_label %||% NA_character_)[1]
  end_label <- as.character(end_label %||% NA_character_)[1]
  min_amplitude <- suppressWarnings(as.numeric(min_amplitude)[1])
  min_bins_observed <- suppressWarnings(as.integer(min_bins_observed)[1])
  min_features_detected <- suppressWarnings(as.integer(min_features_detected)[1])
  min_detected_fraction <- suppressWarnings(as.numeric(min_detected_fraction)[1])
  min_mean_score <- suppressWarnings(as.numeric(min_mean_score)[1])
  min_max_score <- suppressWarnings(as.numeric(min_max_score)[1])
  min_signal_to_noise <- suppressWarnings(as.numeric(min_signal_to_noise)[1])
  low_abundance_max_score <- suppressWarnings(as.numeric(low_abundance_max_score)[1])
  smooth_method <- .trajectory_s2_smooth_fun(smooth_method)
  loess_span <- suppressWarnings(as.numeric(loess_span)[1])
  bootstrap_n <- suppressWarnings(as.integer(bootstrap_n)[1])
  bootstrap_seed <- suppressWarnings(as.integer(bootstrap_seed)[1])
  ci_level <- suppressWarnings(as.numeric(ci_level)[1])
  class_min_amplitude <- suppressWarnings(as.numeric(class_min_amplitude)[1])
  n_permutations <- suppressWarnings(as.integer(n_permutations)[1])
  permutation_seed <- suppressWarnings(as.integer(permutation_seed)[1])
  p_adjust_method <- as.character(p_adjust_method %||% "BH")[1]
  if (!nzchar(p_adjust_method)) p_adjust_method <- "BH"
  .trajectory_assert(is.finite(min_cells_per_bin) && min_cells_per_bin >= 1L,
             "trajectory_geneset_dynamics: min_cells_per_bin must be >= 1.")
  .trajectory_assert(is.finite(min_set_size) && min_set_size >= 1L,
             "trajectory_geneset_dynamics: min_set_size must be >= 1.")
  .trajectory_assert(is.finite(max_set_size) && max_set_size >= min_set_size,
             "trajectory_geneset_dynamics: max_set_size must be >= min_set_size.")
  if (!is.finite(top_n) || top_n < 1L) top_n <- 50L
  if (!is.finite(on_off_min_delta) || on_off_min_delta < 0) on_off_min_delta <- 0
  if (!is.finite(monotonic_min_abs_rho) || monotonic_min_abs_rho < 0) monotonic_min_abs_rho <- 0.5
  if (!is.finite(top_driver_genes) || top_driver_genes < 0L) top_driver_genes <- 5L
  if (!is.finite(min_amplitude) || min_amplitude < 0) min_amplitude <- 0
  if (!is.finite(min_bins_observed) || min_bins_observed < 1L) min_bins_observed <- 5L
  if (!is.finite(min_features_detected) || min_features_detected < 0L) min_features_detected <- 1L
  if (!is.finite(min_detected_fraction) || min_detected_fraction < 0) min_detected_fraction <- 0
  if (min_detected_fraction > 1) min_detected_fraction <- 1
  if (!is.finite(min_mean_score) || min_mean_score < 0) min_mean_score <- 0
  if (!is.finite(min_max_score) || min_max_score < 0) min_max_score <- 0
  if (!is.finite(min_signal_to_noise) || min_signal_to_noise < 0) min_signal_to_noise <- 0
  if (!is.finite(low_abundance_max_score) || low_abundance_max_score < 0) low_abundance_max_score <- 0.1
  if (!is.finite(loess_span) || loess_span <= 0 || loess_span > 1) loess_span <- 0.5
  if (!is.finite(bootstrap_n) || bootstrap_n < 0L) bootstrap_n <- 0L
  if (!is.finite(bootstrap_seed)) bootstrap_seed <- 1L
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) ci_level <- 0.95
  if (!is.finite(class_min_amplitude) || class_min_amplitude < 0) class_min_amplitude <- 0
  if (!is.finite(n_permutations) || n_permutations < 0L) n_permutations <- 0L
  if (!is.finite(permutation_seed)) permutation_seed <- 1L
  .trajectory_assert(p_adjust_method %in% stats::p.adjust.methods,
             "trajectory_geneset_dynamics: p_adjust_method must be one of stats::p.adjust.methods.")

  universe <- as.character(rownames(matrix))
  .trajectory_assert(length(universe) > 0L && !anyDuplicated(universe),
    "trajectory_geneset_dynamics: no feature universe for matrix_id.")
  selected_sets <- .trajectory_s2_registry_sets(registry, namespace,
    feature_sets, universe, min_set_size, max_set_size)
  set_names <- selected_sets$names
  set_members <- selected_sets$members

  pt_res <- if (is.null(model_source)) .trajectory_s2_explicit_source(
    pseudotime = pseudotime, pseudotime_key = pseudotime_key,
    membership = membership, lineage = lineage,
    cells = effective_cells, pseudotime_artifact_id = pseudotime_artifact_id,
    membership_artifact_id = membership_artifact_id) else
    .trajectory_s2_model_source(model_source, lineage, effective_cells,
      membership)
  lineages <- pt_res$lineages
  pt_meta <- pt_res$meta

  curve_rows <- list()
  dyn_rows <- list()
  driver_rows <- list()
  driver_curve_rows <- list()
  score_rows <- list()
  adjustment_rows <- list()
  fit_rows <- list()
  ci <- 0L
  di <- 0L
  gi <- 0L
  dci <- 0L
  sci <- 0L
  ai <- 0L
  si <- 0L
  if (n_permutations > 0L) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv,
      inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(permutation_seed)
  }

  for (lin in names(lineages)) {
    df_lin <- lineages[[lin]]
    pt <- suppressWarnings(as.numeric(df_lin$pseudotime))
    names(pt) <- as.character(df_lin$cell_id)
    keep <- is.finite(pt)
    pt <- pt[keep]
    if (identical(pseudotime_direction, "invert")) {
      pt <- max(pt, na.rm = TRUE) - pt
    }
    cells_lin <- as.character(names(pt))
    .trajectory_assert(length(cells_lin) >= 10L,
               paste0("trajectory_geneset_dynamics: lineage '", lin, "' has <10 finite pseudotime cells."))
    .trajectory_assert(stats::sd(pt) > 0,
               paste0("trajectory_geneset_dynamics: lineage '", lin, "' has zero pseudotime variance."))

    bb <- .trajectory_s2_bins(pt, n_bins = n_bins)
    bins <- bb$bins
    n_bins_used <- bb$n_bins
    n_sets_with_rows <- 0L
    n_curve_rows <- 0L
    scored <- .trajectory_geneset_score_lineage(matrix, set_members,
      cells_lin, pt, depth, adjust_depth, depth_field, score_method,
      lin, label = "trajectory_geneset_dynamics")
    adjusted_matrix <- scored$matrix
    sci <- sci + 1L
    score_rows[[sci]] <- scored$score_table
    if (nrow(scored$diagnostics)) {
      ai <- ai + 1L
      adjustment_rows[[ai]] <- scored$diagnostics
    }

    for (sn in names(set_members)) {
      feats <- set_members[[sn]]
      raw_m <- as.matrix(matrix[feats, cells_lin, drop = FALSE])
      m <- adjusted_matrix[feats, , drop = FALSE]
      n_features_detected <- scored$n_features_detected[[sn]]
      detected_fraction <- scored$detected_fraction[[sn]]
      sig <- scored$scores[[sn]]

      agg <- .trajectory_s2_curve(
        entity_id = sn,
        lineage_id = lin,
        pt = pt,
        values = as.numeric(sig),
        bins = bins,
        n_bins_used = n_bins_used,
        min_cells_per_bin = min_cells_per_bin,
        summary_fun = summary_fun,
        smooth_method = smooth_method,
        loess_span = loess_span,
        bootstrap_n = bootstrap_n,
        ci_level = ci_level,
        seed = bootstrap_seed + di + 1L
      )
      if (is.null(agg) || !nrow(agg)) next

      agg$feature_set <- as.character(agg$entity_id)
      agg$namespace <- as.character(namespace)
      agg$n_features <- as.integer(length(feats))
      agg$n_features_detected <- as.integer(n_features_detected)
      agg$detected_fraction <- as.numeric(detected_fraction)
      agg$label <- as.character(agg$feature_set)
      agg <- agg[, c("feature_set", "namespace", "lineage_id", "bin", "pseudotime", "pseudotime_min",
                     "pseudotime_max", "value", "value_raw", "value_smooth", "value_se", "value_ci_low",
                     "value_ci_high", "n_cells", "n_features", "n_features_detected",
                     "detected_fraction", "label"), drop = FALSE]

      mm <- .trajectory_s2_dynamic_metrics(agg)
      null_scores <- .trajectory_s2_null_dynamic_scores(
        entity_id = sn,
        lineage_id = lin,
        pt = pt,
        values = as.numeric(sig),
        bins = bins,
        min_cells_per_bin = min_cells_per_bin,
        summary_fun = summary_fun,
        smooth_method = smooth_method,
        loess_span = loess_span,
        n_permutations = n_permutations
      )
      null_scores <- null_scores[is.finite(null_scores)]
      p_dynamic <- if (length(null_scores) > 0L && is.finite(mm$dynamic_score)) {
        (1 + sum(null_scores >= mm$dynamic_score, na.rm = TRUE)) / (1 + length(null_scores))
      } else {
        NA_real_
      }
      null_mean_dynamic_score <- if (length(null_scores)) mean(null_scores, na.rm = TRUE) else NA_real_
      null_sd_dynamic_score <- if (length(null_scores) >= 2L) stats::sd(null_scores, na.rm = TRUE) else NA_real_
      z_dynamic_score <- if (is.finite(null_sd_dynamic_score) && null_sd_dynamic_score > 0) {
        (mm$dynamic_score - null_mean_dynamic_score) / null_sd_dynamic_score
      } else {
        NA_real_
      }
      curve_class <- .trajectory_s2_curve_class(
        late_minus_early = mm$late_minus_early,
        spearman_rho = mm$spearman_rho,
        amplitude = mm$amplitude,
        peak_bin = mm$peak_bin,
        trough_bin = mm$trough_bin,
        n_bins_observed = mm$n_bins_observed,
        on_off_min_delta = on_off_min_delta,
        monotonic_min_abs_rho = monotonic_min_abs_rho,
        class_min_amplitude = class_min_amplitude
      )
      filter_reasons <- character(0)
      if (!is.finite(mm$amplitude) || mm$amplitude < min_amplitude) filter_reasons <- c(filter_reasons, "low_amplitude")
      if (!is.finite(mm$n_bins_observed) || mm$n_bins_observed < min_bins_observed) filter_reasons <- c(filter_reasons, "too_few_bins")
      if (n_features_detected < min_features_detected) filter_reasons <- c(filter_reasons, "too_few_detected_features")
      if (detected_fraction < min_detected_fraction) filter_reasons <- c(filter_reasons, "low_detected_fraction")
      if (min_mean_score > 0 && (!is.finite(mm$mean_score) || mm$mean_score < min_mean_score)) filter_reasons <- c(filter_reasons, "low_mean_score")
      if (min_max_score > 0 && (!is.finite(mm$max_score) || mm$max_score < min_max_score)) filter_reasons <- c(filter_reasons, "low_max_score")
      if (min_signal_to_noise > 0 && (!is.finite(mm$signal_to_noise) || mm$signal_to_noise < min_signal_to_noise)) filter_reasons <- c(filter_reasons, "low_signal_to_noise")
      passes_filters <- length(filter_reasons) == 0L
      low_abundance_flag <- is.finite(mm$max_score) && mm$max_score < low_abundance_max_score
      di <- di + 1L
      dyn_rows[[di]] <- data.frame(
        feature_set = as.character(sn),
        namespace = as.character(namespace),
        lineage_id = as.character(lin),
        curve_class = as.character(curve_class),
        n_features = as.integer(length(feats)),
        n_features_detected = as.integer(n_features_detected),
        detected_fraction = as.numeric(detected_fraction),
        n_bins_observed = as.integer(mm$n_bins_observed),
        passes_filters = as.logical(passes_filters),
        filter_reason = if (passes_filters) "" else paste(filter_reasons, collapse = ";"),
        amplitude = as.numeric(mm$amplitude),
        min_score = as.numeric(mm$min_score),
        max_score = as.numeric(mm$max_score),
        mean_score = as.numeric(mm$mean_score),
        auc_score = as.numeric(mm$auc_score),
        relative_amplitude = as.numeric(mm$relative_amplitude),
        sd_value = as.numeric(mm$sd_value),
        mean_se = as.numeric(mm$mean_se),
        signal_to_noise = as.numeric(mm$signal_to_noise),
        low_abundance_flag = as.logical(low_abundance_flag),
        early_mean = as.numeric(mm$early_mean),
        late_mean = as.numeric(mm$late_mean),
        late_minus_early = as.numeric(mm$late_minus_early),
        spearman_rho = as.numeric(mm$spearman_rho),
        abs_spearman_rho = as.numeric(abs(mm$spearman_rho)),
        dynamic_score = as.numeric(mm$dynamic_score),
        p_dynamic = as.numeric(p_dynamic),
        q_dynamic = NA_real_,
        null_mean_dynamic_score = as.numeric(null_mean_dynamic_score),
        null_sd_dynamic_score = as.numeric(null_sd_dynamic_score),
        z_dynamic_score = as.numeric(z_dynamic_score),
        n_permutations = as.integer(length(null_scores)),
        peak_bin = as.integer(mm$peak_bin),
        trough_bin = as.integer(mm$trough_bin),
        peak_pseudotime = as.numeric(mm$peak_pseudotime),
        trough_pseudotime = as.numeric(mm$trough_pseudotime),
        stringsAsFactors = FALSE
      )

      if (top_driver_genes > 0L && length(feats) > 0L) {
        gene_rows <- list()
        gene_curve_by_feature <- list()
        for (g in seq_len(nrow(m))) {
          gagg <- .trajectory_s2_curve(
            entity_id = rownames(m)[g],
            lineage_id = lin,
            pt = pt,
            values = as.numeric(m[g, ]),
            bins = bins,
            n_bins_used = n_bins_used,
            min_cells_per_bin = min_cells_per_bin,
            summary_fun = summary_fun,
            smooth_method = smooth_method,
            loess_span = loess_span,
            bootstrap_n = bootstrap_n,
            ci_level = ci_level,
            seed = bootstrap_seed + di + g + 1L
          )
          if (is.null(gagg) || !nrow(gagg)) next
          gm <- .trajectory_s2_dynamic_metrics(gagg)
          gene_curve_by_feature[[as.character(rownames(m)[g])]] <- gagg
          gene_rows[[length(gene_rows) + 1L]] <- data.frame(
            feature_set = as.character(sn),
            namespace = as.character(namespace),
            lineage_id = as.character(lin),
            feature_id = as.character(rownames(m)[g]),
            n_bins_observed = as.integer(gm$n_bins_observed),
            amplitude = as.numeric(gm$amplitude),
            late_minus_early = as.numeric(gm$late_minus_early),
            spearman_rho = as.numeric(gm$spearman_rho),
            dynamic_score = as.numeric(gm$dynamic_score),
            stringsAsFactors = FALSE
          )
        }
        if (length(gene_rows)) {
          gt <- do.call(rbind, gene_rows)
          gt <- gt[order(-gt$dynamic_score, -gt$amplitude), , drop = FALSE]
          gt$driver_rank <- seq_len(nrow(gt))
          gt <- utils::head(gt, min(top_driver_genes, nrow(gt)))
          gi <- gi + 1L
          driver_rows[[gi]] <- gt
          for (ii in seq_len(nrow(gt))) {
            fid <- as.character(gt$feature_id[ii])
            gcurve <- gene_curve_by_feature[[fid]]
            if (is.null(gcurve) || !nrow(gcurve)) next
            gcurve$feature_set <- as.character(sn)
            gcurve$namespace <- as.character(namespace)
            gcurve$feature_id <- fid
            gcurve$driver_rank <- as.integer(gt$driver_rank[ii])
            gcurve$label <- fid
            gcurve <- gcurve[, c("feature_set", "namespace", "feature_id", "driver_rank", "lineage_id",
                                 "bin", "pseudotime", "pseudotime_min", "pseudotime_max", "value",
                                 "value_raw", "value_smooth", "value_se", "value_ci_low", "value_ci_high",
                                 "n_cells", "label"), drop = FALSE]
            dci <- dci + 1L
            driver_curve_rows[[dci]] <- gcurve
          }
        }
      }

      ci <- ci + 1L
      curve_rows[[ci]] <- agg
      n_sets_with_rows <- n_sets_with_rows + 1L
      n_curve_rows <- n_curve_rows + nrow(agg)
    }

    si <- si + 1L
    fit_rows[[si]] <- data.frame(
      lineage_id = as.character(lin),
      namespace = as.character(namespace),
      pseudotime_source_kind = as.character(pt_meta$pseudotime_source_kind),
      pseudotime_artifact_id = as.character(pt_meta$pseudotime_artifact_id),
      pseudotime_key = as.character(pt_meta$pseudotime_key),
      membership_artifact_id = as.character(pt_meta$membership_artifact_id),
      n_pseudotime_cells = as.integer(pt_meta$n_pseudotime_cells),
      n_membership_cells = as.integer(pt_meta$n_membership_cells),
      n_cells_used = as.integer(length(cells_lin)),
      n_feature_sets_input = as.integer(length(set_names)),
      n_feature_sets_screened = as.integer(length(set_members)),
      n_feature_sets_with_rows = as.integer(n_sets_with_rows),
      n_rows = as.integer(n_curve_rows),
      n_bins = as.integer(n_bins_used),
      score_method = as.character(score_method),
      summary_fun = as.character(summary_fun),
      min_cells_per_bin = as.integer(min_cells_per_bin),
      min_set_size = as.integer(min_set_size),
      max_set_size = as.integer(max_set_size),
      on_off_min_delta = as.numeric(on_off_min_delta),
      monotonic_min_abs_rho = as.numeric(monotonic_min_abs_rho),
      pseudotime_direction = as.character(pseudotime_direction),
      start_label = as.character(start_label),
      end_label = as.character(end_label),
      min_amplitude = as.numeric(min_amplitude),
      min_bins_observed = as.integer(min_bins_observed),
      min_features_detected = as.integer(min_features_detected),
      min_detected_fraction = as.numeric(min_detected_fraction),
      min_mean_score = as.numeric(min_mean_score),
      min_max_score = as.numeric(min_max_score),
      min_signal_to_noise = as.numeric(min_signal_to_noise),
      low_abundance_max_score = as.numeric(low_abundance_max_score),
      smooth_method = as.character(smooth_method),
      loess_span = as.numeric(loess_span),
      bootstrap_n = as.integer(bootstrap_n),
      ci_level = as.numeric(ci_level),
      class_min_amplitude = as.numeric(class_min_amplitude),
      n_permutations = as.integer(n_permutations),
      permutation_seed = as.integer(permutation_seed),
      p_adjust_method = as.character(p_adjust_method),
      top_driver_genes = as.integer(top_driver_genes),
      matrix_id = as.character(matrix_id),
      stringsAsFactors = FALSE
    )
  }

  curve_tbl <- if (ci == 0L) {
    data.frame(
      feature_set = character(0), namespace = character(0), lineage_id = character(0),
      bin = integer(0), pseudotime = numeric(0), pseudotime_min = numeric(0),
      pseudotime_max = numeric(0), value = numeric(0), value_raw = numeric(0),
      value_smooth = numeric(0), value_se = numeric(0), value_ci_low = numeric(0),
      value_ci_high = numeric(0), n_cells = integer(0), n_features = integer(0),
      n_features_detected = integer(0), detected_fraction = numeric(0),
      label = character(0), stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, curve_rows[seq_len(ci)])
    rownames(out) <- NULL
    out
  }
  dyn_tbl <- if (di == 0L) {
    data.frame(
      feature_set = character(0), namespace = character(0), lineage_id = character(0),
      curve_class = character(0), n_features = integer(0), n_features_detected = integer(0),
      detected_fraction = numeric(0), n_bins_observed = integer(0), passes_filters = logical(0),
      filter_reason = character(0), amplitude = numeric(0),
      min_score = numeric(0), max_score = numeric(0), mean_score = numeric(0),
      auc_score = numeric(0), relative_amplitude = numeric(0),
      sd_value = numeric(0), mean_se = numeric(0), signal_to_noise = numeric(0),
      low_abundance_flag = logical(0), early_mean = numeric(0), late_mean = numeric(0),
      late_minus_early = numeric(0), spearman_rho = numeric(0), abs_spearman_rho = numeric(0),
      dynamic_score = numeric(0), peak_bin = integer(0), trough_bin = integer(0),
      p_dynamic = numeric(0), q_dynamic = numeric(0), null_mean_dynamic_score = numeric(0),
      null_sd_dynamic_score = numeric(0), z_dynamic_score = numeric(0), n_permutations = integer(0),
      peak_pseudotime = numeric(0), trough_pseudotime = numeric(0), rank = integer(0),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, dyn_rows[seq_len(di)])
    out$q_dynamic <- NA_real_
    finite_p <- is.finite(out$p_dynamic)
    if (any(finite_p)) out$q_dynamic[finite_p] <- stats::p.adjust(out$p_dynamic[finite_p], method = p_adjust_method)
    out <- out[order(as.character(out$lineage_id), -out$dynamic_score, -out$amplitude, -out$sd_value), , drop = FALSE]
    out$rank <- stats::ave(seq_len(nrow(out)), as.character(out$lineage_id), FUN = seq_along)
    rownames(out) <- NULL
    out
  }
  eligible_dyn <- dyn_tbl[dyn_tbl$passes_filters %in% TRUE, , drop = FALSE]
  if (!nrow(eligible_dyn)) eligible_dyn <- dyn_tbl
  driver_tbl <- if (gi == 0L) {
    data.frame(
      feature_set = character(0), namespace = character(0), lineage_id = character(0),
      feature_id = character(0), n_bins_observed = integer(0), amplitude = numeric(0),
      late_minus_early = numeric(0), spearman_rho = numeric(0), dynamic_score = numeric(0),
      driver_rank = integer(0), stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, driver_rows[seq_len(gi)])
    rownames(out) <- NULL
    out
  }
  driver_curve_tbl <- if (dci == 0L) {
    data.frame(
      feature_set = character(0), namespace = character(0), feature_id = character(0),
      driver_rank = integer(0), lineage_id = character(0), bin = integer(0),
      pseudotime = numeric(0), pseudotime_min = numeric(0), pseudotime_max = numeric(0),
      value = numeric(0), value_raw = numeric(0), value_smooth = numeric(0),
      value_se = numeric(0), value_ci_low = numeric(0), value_ci_high = numeric(0),
      n_cells = integer(0), label = character(0), stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, driver_curve_rows[seq_len(dci)])
    rownames(out) <- NULL
    out
  }
  top_tbl <- utils::head(eligible_dyn[order(as.character(eligible_dyn$lineage_id), eligible_dyn$rank), , drop = FALSE], top_n)
  fit_summary <- do.call(rbind, fit_rows[seq_len(si)])
  rownames(fit_summary) <- NULL
  score_table <- do.call(rbind, score_rows[seq_len(sci)])
  rownames(score_table) <- NULL
  adjustment_table <- if (ai == 0L)
    .trajectory_geneset_empty_adjustment_table() else
      do.call(rbind, adjustment_rows[seq_len(ai)])
  rownames(adjustment_table) <- NULL

  list(
    outputs = list(
      signature_table = scworkbench::wb_prim_obs_table(
        table = curve_tbl,
        table_type = "trajectory_geneset_dynamics_curve",
        entity_kind = "trajectory_bin",
        meta = list(
          matrix_id = as.character(matrix_id),
          namespace = as.character(namespace),
          pseudotime_source_kind = as.character(pt_meta$pseudotime_source_kind),
          pseudotime_artifact_id = as.character(pt_meta$pseudotime_artifact_id),
          pseudotime_key = as.character(pt_meta$pseudotime_key),
          membership_artifact_id = as.character(pt_meta$membership_artifact_id),
          score_method = as.character(score_method),
          summary_fun = as.character(summary_fun),
          min_cells_per_bin = as.integer(min_cells_per_bin),
          min_set_size = as.integer(min_set_size),
          max_set_size = as.integer(max_set_size),
          on_off_min_delta = as.numeric(on_off_min_delta),
          monotonic_min_abs_rho = as.numeric(monotonic_min_abs_rho),
          pseudotime_direction = as.character(pseudotime_direction),
          start_label = as.character(start_label),
          end_label = as.character(end_label),
          min_amplitude = as.numeric(min_amplitude),
          min_bins_observed = as.integer(min_bins_observed),
          min_features_detected = as.integer(min_features_detected),
          min_detected_fraction = as.numeric(min_detected_fraction),
          min_mean_score = as.numeric(min_mean_score),
          min_max_score = as.numeric(min_max_score),
          min_signal_to_noise = as.numeric(min_signal_to_noise),
          low_abundance_max_score = as.numeric(low_abundance_max_score),
          smooth_method = as.character(smooth_method),
          loess_span = as.numeric(loess_span),
          bootstrap_n = as.integer(bootstrap_n),
          ci_level = as.numeric(ci_level),
          n_permutations = as.integer(n_permutations),
          permutation_seed = as.integer(permutation_seed),
          p_adjust_method = as.character(p_adjust_method),
          n_bins = as.integer(unique(fit_summary$n_bins %||% n_bins)[1] %||% n_bins)
        )
      ),
      dynamics_table = scworkbench::wb_prim_obs_table(
        table = dyn_tbl,
        table_type = "trajectory_geneset_dynamics",
        entity_kind = "feature_set",
        meta = list(
          matrix_id = as.character(matrix_id),
          namespace = as.character(namespace),
          pseudotime_source_kind = as.character(pt_meta$pseudotime_source_kind),
          pseudotime_artifact_id = as.character(pt_meta$pseudotime_artifact_id),
          pseudotime_key = as.character(pt_meta$pseudotime_key),
          membership_artifact_id = as.character(pt_meta$membership_artifact_id),
          dynamic_score = "amplitude + sd_value + abs(late_minus_early) + abs(spearman_rho) * sd_value",
          abundance_metrics = "mean_score, max_score, min_score, auc_score, relative_amplitude, mean_se, signal_to_noise, low_abundance_flag",
          empirical_significance = "p_dynamic and q_dynamic from permuting cell-level scores relative to fixed pseudotime/bin structure; disabled when n_permutations=0",
          curve_class = "turns_on/turns_off/transient_peak/transient_dip/weak_on/weak_off/flat_or_noisy",
          passes_filters = "min_amplitude, min_bins_observed, min_features_detected, min_detected_fraction, min_mean_score, min_max_score, min_signal_to_noise"
        )
      ),
      top_sets = scworkbench::wb_prim_obs_table(
        table = top_tbl,
        table_type = "trajectory_geneset_dynamics_top",
        entity_kind = "feature_set",
        meta = list(
          matrix_id = as.character(matrix_id),
          namespace = as.character(namespace),
          top_n = as.integer(top_n)
        )
      ),
      driver_table = scworkbench::wb_prim_obs_table(
        table = driver_tbl,
        table_type = "trajectory_geneset_dynamics_drivers",
        entity_kind = "feature",
        meta = list(
          matrix_id = as.character(matrix_id),
          namespace = as.character(namespace),
          top_driver_genes = as.integer(top_driver_genes)
        )
      ),
      driver_curve_table = scworkbench::wb_prim_obs_table(
        table = driver_curve_tbl,
        table_type = "trajectory_geneset_dynamics_driver_curve",
        entity_kind = "trajectory_bin",
        meta = list(
          matrix_id = as.character(matrix_id),
          namespace = as.character(namespace),
          top_driver_genes = as.integer(top_driver_genes),
          smooth_method = as.character(smooth_method),
          bootstrap_n = as.integer(bootstrap_n)
        )
      ),
      adjusted_scores = scworkbench::wb_prim_cell_vector_family(
        key = score_table$key, cell_id = score_table$cell_id,
        value = score_table$value, value_type = "continuous",
        meta = list(
          depth_adjusted = adjust_depth,
          depth_field = depth_field,
          depth_artifact_id = as.character(depth_artifact_id),
          adjustment_spline_df = if (adjust_depth) 4L else NA_integer_,
          score_key_format = "lineage_id::feature_set",
          effective_score_method = if (adjust_depth) "adjusted_z_mean" else
            as.character(score_method))),
      adjustment_diagnostics = scworkbench::wb_prim_obs_table(
        table = adjustment_table,
        table_type = "trajectory_geneset_depth_adjustment",
        entity_kind = "feature",
        meta = list(matrix_id = as.character(matrix_id),
          depth_adjusted = adjust_depth, depth_field = depth_field,
          depth_artifact_id = as.character(depth_artifact_id),
          adjustment_spline_df = if (adjust_depth) 4L else NA_integer_)),
      fit_summary = scworkbench::wb_prim_obs_table(
        table = fit_summary,
        table_type = "trajectory_geneset_dynamics_summary",
        entity_kind = "lineage"
      )
    ),
    params = list(
      matrix_id = as.character(matrix_id),
      namespace = as.character(namespace),
      pseudotime_source_kind = as.character(pt_meta$pseudotime_source_kind),
      pseudotime_artifact_id = as.character(pt_meta$pseudotime_artifact_id),
      pseudotime_key = as.character(pt_meta$pseudotime_key),
      membership_artifact_id = as.character(pt_meta$membership_artifact_id),
      n_feature_sets_input = as.integer(length(set_names)),
      n_feature_sets_screened = as.integer(length(set_members)),
      n_lineages = as.integer(length(lineages)),
      n_bins = as.integer(unique(fit_summary$n_bins %||% n_bins)[1] %||% n_bins),
      score_method = as.character(score_method),
      effective_score_method = if (adjust_depth) "adjusted_z_mean" else
        as.character(score_method),
      adjust_depth = adjust_depth,
      depth_field = depth_field,
      depth_artifact_id = as.character(depth_artifact_id),
      adjustment_spline_df = if (adjust_depth) 4L else NA_integer_,
      summary_fun = as.character(summary_fun),
      min_cells_per_bin = as.integer(min_cells_per_bin),
      min_set_size = as.integer(min_set_size),
      max_set_size = as.integer(max_set_size),
      top_n = as.integer(top_n),
      on_off_min_delta = as.numeric(on_off_min_delta),
      monotonic_min_abs_rho = as.numeric(monotonic_min_abs_rho),
      top_driver_genes = as.integer(top_driver_genes),
      pseudotime_direction = as.character(pseudotime_direction),
      start_label = as.character(start_label),
      end_label = as.character(end_label),
      min_amplitude = as.numeric(min_amplitude),
      min_bins_observed = as.integer(min_bins_observed),
      min_features_detected = as.integer(min_features_detected),
      min_detected_fraction = as.numeric(min_detected_fraction),
      min_mean_score = as.numeric(min_mean_score),
      min_max_score = as.numeric(min_max_score),
      min_signal_to_noise = as.numeric(min_signal_to_noise),
      low_abundance_max_score = as.numeric(low_abundance_max_score),
      smooth_method = as.character(smooth_method),
      loess_span = as.numeric(loess_span),
      bootstrap_n = as.integer(bootstrap_n),
      ci_level = as.numeric(ci_level),
      class_min_amplitude = as.numeric(class_min_amplitude),
      n_permutations = as.integer(n_permutations),
      permutation_seed = as.integer(permutation_seed),
      p_adjust_method = as.character(p_adjust_method)
    )
  )
}
