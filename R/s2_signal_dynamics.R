# Pack-side S2 method, with source selection supplied by authenticated readers.
.trajectory_s2_signal_dynamics <- function(effective_cells,
                                             signal,
                                             pseudotime = NULL,
                                             model_source = NULL,
                                             pseudotime_key = NULL,
                                             membership = NULL,
                                             pseudotime_artifact_id = NULL,
                                             membership_artifact_id = NULL,
                                             lineage = NULL,
                                             signal_keys = NULL,
                                             n_bins = 50L,
                                             summary_fun = "mean",
                                             min_cells_per_bin = 5L,
                                             top_n = 50L,
                                             on_off_min_delta = 0,
                                             monotonic_min_abs_rho = 0.5,
                                             pseudotime_direction = "as_is",
                                             start_label = NULL,
                                             end_label = NULL,
                                             min_amplitude = 0,
                                             min_bins_observed = 5L,
                                             smooth_method = "none",
                                             loess_span = 0.5,
                                             bootstrap_n = 0L,
                                             bootstrap_seed = 1L,
                                             ci_level = 0.95,
                                             class_min_amplitude = 0) {
  .trajectory_assert(xor(is.null(pseudotime), is.null(model_source)),
    "trajectory_signal_dynamics requires exactly one pseudotime source.")
  if (is.null(model_source)) .trajectory_kind(pseudotime,
    "cell_vector_family")
  .trajectory_assert(inherits(signal, "wb_primitive"), "trajectory_signal_dynamics: signal must be a primitive.")

  if (identical(as.character(signal$kind %||% "")[1], "cell_vector")) {
    label <- as.character(signal$label %||% "signal")[1]
    if (!nzchar(label)) label <- "signal"
    signal <- scworkbench::wb_prim_cell_vector_family(
      key = rep(label, length(signal$cell_id)),
      cell_id = as.character(signal$cell_id),
      value = as.numeric(signal$value),
      value_type = as.character(signal$value_type %||% "continuous")[1]
    )
  }
  .trajectory_kind(signal, "cell_vector_family")

  summary_fun <- .trajectory_s2_summary_fun(summary_fun)
  min_cells_per_bin <- suppressWarnings(as.integer(min_cells_per_bin))
  top_n <- suppressWarnings(as.integer(top_n))
  on_off_min_delta <- suppressWarnings(as.numeric(on_off_min_delta))
  monotonic_min_abs_rho <- suppressWarnings(as.numeric(monotonic_min_abs_rho))
  pseudotime_direction <- as.character(pseudotime_direction %||% "as_is")[1]
  if (!nzchar(pseudotime_direction)) pseudotime_direction <- "as_is"
  .trajectory_assert(pseudotime_direction %in% c("as_is", "invert"),
             "trajectory_signal_dynamics: pseudotime_direction must be 'as_is' or 'invert'.")
  start_label <- as.character(start_label %||% NA_character_)[1]
  end_label <- as.character(end_label %||% NA_character_)[1]
  min_amplitude <- suppressWarnings(as.numeric(min_amplitude)[1])
  min_bins_observed <- suppressWarnings(as.integer(min_bins_observed)[1])
  smooth_method <- .trajectory_s2_smooth_fun(smooth_method)
  loess_span <- suppressWarnings(as.numeric(loess_span)[1])
  bootstrap_n <- suppressWarnings(as.integer(bootstrap_n)[1])
  bootstrap_seed <- suppressWarnings(as.integer(bootstrap_seed)[1])
  ci_level <- suppressWarnings(as.numeric(ci_level)[1])
  class_min_amplitude <- suppressWarnings(as.numeric(class_min_amplitude)[1])
  .trajectory_assert(is.finite(min_cells_per_bin) && min_cells_per_bin >= 1L,
             "trajectory_signal_dynamics: min_cells_per_bin must be >= 1.")
  if (!is.finite(top_n) || top_n < 1L) top_n <- 50L
  if (!is.finite(on_off_min_delta) || on_off_min_delta < 0) on_off_min_delta <- 0
  if (!is.finite(monotonic_min_abs_rho) || monotonic_min_abs_rho < 0) monotonic_min_abs_rho <- 0.5
  if (!is.finite(min_amplitude) || min_amplitude < 0) min_amplitude <- 0
  if (!is.finite(min_bins_observed) || min_bins_observed < 1L) min_bins_observed <- 5L
  if (!is.finite(loess_span) || loess_span <= 0 || loess_span > 1) loess_span <- 0.5
  if (!is.finite(bootstrap_n) || bootstrap_n < 0L) bootstrap_n <- 0L
  if (!is.finite(bootstrap_seed)) bootstrap_seed <- 1L
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) ci_level <- 0.95
  if (!is.finite(class_min_amplitude) || class_min_amplitude < 0) class_min_amplitude <- 0

  sig_tbl <- scworkbench::wb_primitive_to_table(signal, include_context = FALSE)
  .trajectory_assert(all(c("key", "cell_id", "value") %in% colnames(sig_tbl)),
             "trajectory_signal_dynamics: signal family table missing key/cell_id/value.")
  sig_tbl$key <- as.character(sig_tbl$key)
  sig_tbl$cell_id <- as.character(sig_tbl$cell_id)
  sig_tbl$value <- suppressWarnings(as.numeric(sig_tbl$value))
  sig_tbl <- sig_tbl[nzchar(sig_tbl$key) & nzchar(sig_tbl$cell_id) & is.finite(sig_tbl$value), , drop = FALSE]
  .trajectory_assert(nrow(sig_tbl) > 0L, "trajectory_signal_dynamics: signal has no finite rows.")

  if (!is.null(signal_keys)) {
    signal_keys <- unique(as.character(signal_keys))
    signal_keys <- signal_keys[!is.na(signal_keys) & nzchar(signal_keys)]
    .trajectory_assert(length(signal_keys) > 0L, "trajectory_signal_dynamics: signal_keys is empty after cleaning.")
    .trajectory_assert(!length(setdiff(signal_keys, unique(sig_tbl$key))),
      "trajectory_signal_dynamics: requested signal_key not found.")
    sig_tbl <- sig_tbl[sig_tbl$key %in% signal_keys, , drop = FALSE]
    .trajectory_assert(nrow(sig_tbl) > 0L, "trajectory_signal_dynamics: no signal rows remain after signal_keys filtering.")
  }
  signal_keys <- unique(as.character(sig_tbl$key))

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
  fit_rows <- list()
  ci <- 0L
  di <- 0L
  si <- 0L

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
               paste0("trajectory_signal_dynamics: lineage '", lin, "' has <10 finite pseudotime cells."))
    .trajectory_assert(stats::sd(pt) > 0,
               paste0("trajectory_signal_dynamics: lineage '", lin, "' has zero pseudotime variance."))

    bb <- .trajectory_s2_bins(pt, n_bins = n_bins)
    bins <- bb$bins
    n_bins_used <- bb$n_bins
    n_rows_lineage <- 0L
    n_signals_with_rows <- 0L
    sig_lin <- sig_tbl[sig_tbl$cell_id %in% cells_lin, c("key", "cell_id", "value"), drop = FALSE]
    .trajectory_assert(nrow(sig_lin) > 0L,
               paste0("trajectory_signal_dynamics: no signal rows overlap lineage '", lin, "'."))

    for (sig_key in signal_keys) {
      sub <- sig_lin[sig_lin$key == sig_key, , drop = FALSE]
      if (!nrow(sub)) next
      vals <- sub$value[match(cells_lin, sub$cell_id)]
      keep2 <- is.finite(vals) & is.finite(pt) & is.finite(bins)
      if (sum(keep2) < min_cells_per_bin) next

      agg <- .trajectory_s2_curve(
        entity_id = sig_key,
        lineage_id = lin,
        pt = pt[keep2],
        values = vals[keep2],
        bins = bins[keep2],
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

      agg$signal_key <- as.character(agg$entity_id)
      agg$feature_set <- as.character(agg$signal_key)
      agg$label <- as.character(agg$signal_key)
      agg <- agg[, c("signal_key", "feature_set", "lineage_id", "bin", "pseudotime",
                     "pseudotime_min", "pseudotime_max", "value", "value_raw",
                     "value_smooth", "value_se", "value_ci_low", "value_ci_high",
                     "n_cells", "label"), drop = FALSE]

      mm <- .trajectory_s2_dynamic_metrics(agg)
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
      passes_filters <- length(filter_reasons) == 0L

      di <- di + 1L
      dyn_rows[[di]] <- data.frame(
        signal_key = as.character(sig_key),
        feature_set = as.character(sig_key),
        lineage_id = as.character(lin),
        curve_class = as.character(curve_class),
        n_bins_observed = as.integer(mm$n_bins_observed),
        passes_filters = as.logical(passes_filters),
        filter_reason = if (passes_filters) "" else paste(filter_reasons, collapse = ";"),
        amplitude = as.numeric(mm$amplitude),
        sd_value = as.numeric(mm$sd_value),
        early_mean = as.numeric(mm$early_mean),
        late_mean = as.numeric(mm$late_mean),
        late_minus_early = as.numeric(mm$late_minus_early),
        spearman_rho = as.numeric(mm$spearman_rho),
        abs_spearman_rho = as.numeric(abs(mm$spearman_rho)),
        dynamic_score = as.numeric(mm$dynamic_score),
        peak_bin = as.integer(mm$peak_bin),
        trough_bin = as.integer(mm$trough_bin),
        peak_pseudotime = as.numeric(mm$peak_pseudotime),
        trough_pseudotime = as.numeric(mm$trough_pseudotime),
        stringsAsFactors = FALSE
      )

      ci <- ci + 1L
      curve_rows[[ci]] <- agg
      n_signals_with_rows <- n_signals_with_rows + 1L
      n_rows_lineage <- n_rows_lineage + nrow(agg)
    }

    si <- si + 1L
    fit_rows[[si]] <- data.frame(
      lineage_id = as.character(lin),
      pseudotime_source_kind = as.character(pt_meta$pseudotime_source_kind),
      pseudotime_artifact_id = as.character(pt_meta$pseudotime_artifact_id),
      pseudotime_key = as.character(pt_meta$pseudotime_key),
      membership_artifact_id = as.character(pt_meta$membership_artifact_id),
      n_pseudotime_cells = as.integer(pt_meta$n_pseudotime_cells),
      n_membership_cells = as.integer(pt_meta$n_membership_cells),
      n_cells_used = as.integer(length(cells_lin)),
      n_signal_keys_input = as.integer(length(signal_keys)),
      n_signal_keys_with_rows = as.integer(n_signals_with_rows),
      n_rows = as.integer(n_rows_lineage),
      n_bins = as.integer(n_bins_used),
      summary_fun = as.character(summary_fun),
      min_cells_per_bin = as.integer(min_cells_per_bin),
      on_off_min_delta = as.numeric(on_off_min_delta),
      monotonic_min_abs_rho = as.numeric(monotonic_min_abs_rho),
      pseudotime_direction = as.character(pseudotime_direction),
      start_label = as.character(start_label),
      end_label = as.character(end_label),
      min_amplitude = as.numeric(min_amplitude),
      min_bins_observed = as.integer(min_bins_observed),
      smooth_method = as.character(smooth_method),
      loess_span = as.numeric(loess_span),
      bootstrap_n = as.integer(bootstrap_n),
      ci_level = as.numeric(ci_level),
      class_min_amplitude = as.numeric(class_min_amplitude),
      stringsAsFactors = FALSE
    )
  }

  curve_tbl <- if (ci == 0L) {
    data.frame(
      signal_key = character(0), feature_set = character(0), lineage_id = character(0),
      bin = integer(0), pseudotime = numeric(0), pseudotime_min = numeric(0),
      pseudotime_max = numeric(0), value = numeric(0), value_raw = numeric(0),
      value_smooth = numeric(0), value_se = numeric(0), value_ci_low = numeric(0),
      value_ci_high = numeric(0), n_cells = integer(0), label = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, curve_rows[seq_len(ci)])
    rownames(out) <- NULL
    out
  }
  dyn_tbl <- if (di == 0L) {
    data.frame(
      signal_key = character(0), feature_set = character(0), lineage_id = character(0),
      curve_class = character(0), n_bins_observed = integer(0), passes_filters = logical(0),
      filter_reason = character(0), amplitude = numeric(0), sd_value = numeric(0),
      early_mean = numeric(0), late_mean = numeric(0), late_minus_early = numeric(0),
      spearman_rho = numeric(0), abs_spearman_rho = numeric(0), dynamic_score = numeric(0),
      peak_bin = integer(0), trough_bin = integer(0), peak_pseudotime = numeric(0),
      trough_pseudotime = numeric(0), rank = integer(0), stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, dyn_rows[seq_len(di)])
    out <- out[order(as.character(out$lineage_id), -out$dynamic_score, -out$amplitude, -out$sd_value), , drop = FALSE]
    out$rank <- ave(seq_len(nrow(out)), as.character(out$lineage_id), FUN = seq_along)
    rownames(out) <- NULL
    out
  }
  eligible_dyn <- dyn_tbl[dyn_tbl$passes_filters %in% TRUE, , drop = FALSE]
  if (!nrow(eligible_dyn)) eligible_dyn <- dyn_tbl
  top_tbl <- head(eligible_dyn[order(as.character(eligible_dyn$lineage_id), eligible_dyn$rank), , drop = FALSE], top_n)
  fit_summary <- do.call(rbind, fit_rows[seq_len(si)])
  rownames(fit_summary) <- NULL

  list(
    outputs = list(
      signal_table = scworkbench::wb_prim_obs_table(
        table = curve_tbl,
        table_type = "trajectory_signal_dynamics_curve",
        entity_kind = "trajectory_bin",
        meta = list(
          signal_kind = as.character(signal$kind),
          pseudotime_source_kind = as.character(pt_meta$pseudotime_source_kind),
          pseudotime_artifact_id = as.character(pt_meta$pseudotime_artifact_id),
          pseudotime_key = as.character(pt_meta$pseudotime_key),
          membership_artifact_id = as.character(pt_meta$membership_artifact_id),
          n_signal_keys = as.integer(length(signal_keys)),
          summary_fun = as.character(summary_fun),
          min_cells_per_bin = as.integer(min_cells_per_bin),
          on_off_min_delta = as.numeric(on_off_min_delta),
          monotonic_min_abs_rho = as.numeric(monotonic_min_abs_rho),
          pseudotime_direction = as.character(pseudotime_direction),
          start_label = as.character(start_label),
          end_label = as.character(end_label),
          min_amplitude = as.numeric(min_amplitude),
          min_bins_observed = as.integer(min_bins_observed),
          smooth_method = as.character(smooth_method),
          loess_span = as.numeric(loess_span),
          bootstrap_n = as.integer(bootstrap_n),
          ci_level = as.numeric(ci_level),
          n_bins = as.integer(unique(fit_summary$n_bins %||% n_bins)[1] %||% n_bins)
        )
      ),
      dynamics_table = scworkbench::wb_prim_obs_table(
        table = dyn_tbl,
        table_type = "trajectory_signal_dynamics",
        entity_kind = "signal",
        meta = list(
          dynamic_score = "amplitude + sd_value + abs(late_minus_early) + abs(spearman_rho) * sd_value",
          pseudotime_source_kind = as.character(pt_meta$pseudotime_source_kind),
          pseudotime_artifact_id = as.character(pt_meta$pseudotime_artifact_id),
          pseudotime_key = as.character(pt_meta$pseudotime_key),
          membership_artifact_id = as.character(pt_meta$membership_artifact_id),
          curve_class = "turns_on/turns_off/transient_peak/transient_dip/weak_on/weak_off/flat_or_noisy",
          passes_filters = "min_amplitude, min_bins_observed"
        )
      ),
      top_signals = scworkbench::wb_prim_obs_table(
        table = top_tbl,
        table_type = "trajectory_signal_dynamics_top",
        entity_kind = "signal",
        meta = list(top_n = as.integer(top_n))
      ),
      fit_summary = scworkbench::wb_prim_obs_table(
        table = fit_summary,
        table_type = "trajectory_signal_dynamics_summary",
        entity_kind = "lineage"
      )
    ),
    params = list(
      signal_kind = as.character(signal$kind),
      pseudotime_source_kind = as.character(pt_meta$pseudotime_source_kind),
      pseudotime_artifact_id = as.character(pt_meta$pseudotime_artifact_id),
      pseudotime_key = as.character(pt_meta$pseudotime_key),
      membership_artifact_id = as.character(pt_meta$membership_artifact_id),
      signal_keys = as.character(signal_keys),
      n_signal_keys_input = as.integer(length(signal_keys)),
      n_lineages = as.integer(length(lineages)),
      n_bins = as.integer(unique(fit_summary$n_bins %||% n_bins)[1] %||% n_bins),
      summary_fun = as.character(summary_fun),
      min_cells_per_bin = as.integer(min_cells_per_bin),
      top_n = as.integer(top_n),
      on_off_min_delta = as.numeric(on_off_min_delta),
      monotonic_min_abs_rho = as.numeric(monotonic_min_abs_rho),
      pseudotime_direction = as.character(pseudotime_direction),
      start_label = as.character(start_label),
      end_label = as.character(end_label),
      min_amplitude = as.numeric(min_amplitude),
      min_bins_observed = as.integer(min_bins_observed),
      smooth_method = as.character(smooth_method),
      loess_span = as.numeric(loess_span),
      bootstrap_n = as.integer(bootstrap_n),
      ci_level = as.numeric(ci_level),
      class_min_amplitude = as.numeric(class_min_amplitude)
    )
  )
}
