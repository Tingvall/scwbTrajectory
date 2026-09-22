# S2 method equations copied from the frozen core implementation, with only
# pack-local method calls. No bundle or state authority crosses this boundary.

.trajectory_s2_dynamic_metrics <- function(curve) {
  yy <- suppressWarnings(as.numeric(curve$value))
  se <- if ("value_se" %in% colnames(curve)) {
    suppressWarnings(as.numeric(curve$value_se))
  } else {
    rep(NA_real_, length(yy))
  }
  xx <- suppressWarnings(as.numeric(curve$pseudotime))
  bb <- suppressWarnings(as.integer(curve$bin))
  keep <- is.finite(yy) & is.finite(xx) & is.finite(bb)
  yy <- yy[keep]
  se <- se[keep]
  xx <- xx[keep]
  bb <- bb[keep]

  if (!length(yy)) {
    return(list(
      n_bins_observed = 0L,
      amplitude = NA_real_,
      min_score = NA_real_,
      max_score = NA_real_,
      mean_score = NA_real_,
      auc_score = NA_real_,
      relative_amplitude = NA_real_,
      sd_value = NA_real_,
      mean_se = NA_real_,
      signal_to_noise = NA_real_,
      early_mean = NA_real_,
      late_mean = NA_real_,
      late_minus_early = NA_real_,
      spearman_rho = NA_real_,
      dynamic_score = NA_real_,
      peak_bin = NA_integer_,
      trough_bin = NA_integer_,
      peak_pseudotime = NA_real_,
      trough_pseudotime = NA_real_
    ))
  }

  oo <- order(bb, xx)
  yy <- yy[oo]
  xx <- xx[oo]
  bb <- bb[oo]

  n <- length(yy)
  qn <- max(1L, ceiling(n * 0.2))
  early <- yy[seq_len(qn)]
  late <- yy[seq.int(max(1L, n - qn + 1L), n)]
  pk <- which.max(yy)
  tr <- which.min(yy)
  amp <- suppressWarnings(max(yy, na.rm = TRUE) - min(yy, na.rm = TRUE))
  min_score <- suppressWarnings(min(yy, na.rm = TRUE))
  max_score <- suppressWarnings(max(yy, na.rm = TRUE))
  mean_score <- suppressWarnings(mean(yy, na.rm = TRUE))
  auc_score <- mean_score
  relative_amplitude <- if (is.finite(mean_score)) amp / (abs(mean_score) + 1e-8) else NA_real_
  sdv <- if (n >= 2L) suppressWarnings(stats::sd(yy, na.rm = TRUE)) else 0
  mean_se <- suppressWarnings(mean(se[is.finite(se)], na.rm = TRUE))
  if (!is.finite(mean_se)) mean_se <- NA_real_
  signal_to_noise <- if (is.finite(mean_se) && mean_se > 0) amp / mean_se else NA_real_
  rho <- if (n >= 3L && stats::sd(xx) > 0 && stats::sd(yy) > 0) {
    suppressWarnings(stats::cor(xx, yy, method = "spearman", use = "complete.obs"))
  } else {
    NA_real_
  }
  delta <- suppressWarnings(mean(late, na.rm = TRUE) - mean(early, na.rm = TRUE))
  dyn <- amp + sdv + abs(delta) + if (is.finite(rho)) abs(rho) * sdv else 0

  list(
    n_bins_observed = as.integer(n),
    amplitude = as.numeric(amp),
    min_score = as.numeric(min_score),
    max_score = as.numeric(max_score),
    mean_score = as.numeric(mean_score),
    auc_score = as.numeric(auc_score),
    relative_amplitude = as.numeric(relative_amplitude),
    sd_value = as.numeric(sdv),
    mean_se = as.numeric(mean_se),
    signal_to_noise = as.numeric(signal_to_noise),
    early_mean = as.numeric(mean(early, na.rm = TRUE)),
    late_mean = as.numeric(mean(late, na.rm = TRUE)),
    late_minus_early = as.numeric(delta),
    spearman_rho = as.numeric(rho),
    dynamic_score = as.numeric(dyn),
    peak_bin = as.integer(bb[pk]),
    trough_bin = as.integer(bb[tr]),
    peak_pseudotime = as.numeric(xx[pk]),
    trough_pseudotime = as.numeric(xx[tr])
  )
}

.trajectory_s2_curve_class <- function(late_minus_early,
                                               spearman_rho,
                                               amplitude = NA_real_,
                                               peak_bin,
                                               trough_bin,
                                               n_bins_observed,
                                               on_off_min_delta = 0,
                                               monotonic_min_abs_rho = 0.5,
                                               class_min_amplitude = 0) {
  delta <- suppressWarnings(as.numeric(late_minus_early)[1])
  rho <- suppressWarnings(as.numeric(spearman_rho)[1])
  amp <- suppressWarnings(as.numeric(amplitude)[1])
  peak_bin <- suppressWarnings(as.integer(peak_bin)[1])
  trough_bin <- suppressWarnings(as.integer(trough_bin)[1])
  n_bins <- suppressWarnings(as.integer(n_bins_observed)[1])
  if (!is.finite(delta)) delta <- 0
  if (!is.finite(rho)) rho <- 0
  if (!is.finite(amp)) amp <- 0
  if (!is.finite(on_off_min_delta) || on_off_min_delta < 0) on_off_min_delta <- 0
  if (!is.finite(monotonic_min_abs_rho) || monotonic_min_abs_rho < 0) monotonic_min_abs_rho <- 0.5
  if (!is.finite(class_min_amplitude) || class_min_amplitude < 0) class_min_amplitude <- 0

  peak_internal <- is.finite(peak_bin) && is.finite(n_bins) && peak_bin > 1L && peak_bin < n_bins
  trough_internal <- is.finite(trough_bin) && is.finite(n_bins) && trough_bin > 1L && trough_bin < n_bins

  if (amp < class_min_amplitude) return("flat_or_noisy")
  if (abs(delta) >= on_off_min_delta && rho >= monotonic_min_abs_rho) return("turns_on")
  if (abs(delta) >= on_off_min_delta && rho <= -monotonic_min_abs_rho) return("turns_off")
  if (peak_internal && amp >= class_min_amplitude) return("transient_peak")
  if (trough_internal && amp >= class_min_amplitude) return("transient_dip")
  if (delta > 0) return("weak_on")
  if (delta < 0) return("weak_off")
  "flat_or_noisy"
}

.trajectory_s2_smooth_fun <- function(smooth_method) {
  sm <- as.character(smooth_method %||% "none")[1]
  if (!nzchar(sm)) sm <- "none"
  .trajectory_assert(sm %in% c("none", "loess", "spline"),
             "trajectory_geneset_dynamics: smooth_method must be 'none', 'loess', or 'spline'.")
  sm
}

.trajectory_s2_score_values <- function(m, score_method, rank_bg = NULL) {
  score_method <- as.character(score_method %||% "mean")[1]
  if (!nzchar(score_method)) score_method <- "mean"
  .trajectory_assert(score_method %in% c("mean", "sum", "rank_mean"),
             "trajectory_geneset_dynamics: score_method must be 'mean', 'sum', or 'rank_mean'.")
  if (identical(score_method, "sum")) return(colSums(m, na.rm = TRUE))
  if (identical(score_method, "rank_mean")) {
    bg <- rank_bg
    if (is.null(bg)) bg <- m
    .trajectory_assert(is.matrix(bg) || is.data.frame(bg), "trajectory_geneset_dynamics: rank_mean background must be matrix-like.")
    bg <- as.matrix(bg)
    rr <- apply(bg, 2L, function(z) rank(z, ties.method = "average", na.last = "keep"))
    if (is.null(dim(rr))) rr <- matrix(rr, ncol = 1L, dimnames = list(rownames(bg), colnames(bg)))
    idx <- match(rownames(m), rownames(rr))
    idx <- idx[is.finite(idx)]
    .trajectory_assert(length(idx) > 0L, "trajectory_geneset_dynamics: rank_mean set has no overlap with rank background.")
    return(colMeans(rr[idx, , drop = FALSE], na.rm = TRUE))
  }
  colMeans(m, na.rm = TRUE)
}

.trajectory_s2_apply_smoothing <- function(agg, smooth_method = "none", loess_span = 0.5) {
  if (is.null(agg) || !nrow(agg)) return(agg)
  smooth_method <- .trajectory_s2_smooth_fun(smooth_method)
  agg$value_raw <- suppressWarnings(as.numeric(agg$value))
  agg$value_smooth <- agg$value_raw
  if (!identical(smooth_method, "none") && nrow(agg) >= 5L) {
    xx <- suppressWarnings(as.numeric(agg$pseudotime))
    yy <- suppressWarnings(as.numeric(agg$value_raw))
    keep <- is.finite(xx) & is.finite(yy)
    if (sum(keep) >= 5L && stats::sd(xx[keep]) > 0 && stats::sd(yy[keep]) > 0) {
      pred <- rep(NA_real_, length(yy))
      if (identical(smooth_method, "loess")) {
        sp <- suppressWarnings(as.numeric(loess_span)[1])
        if (!is.finite(sp) || sp <= 0 || sp > 1) sp <- 0.5
        fit <- try(stats::loess(yy[keep] ~ xx[keep], span = sp, degree = 2L, surface = "direct"), silent = TRUE)
        if (!inherits(fit, "try-error")) pred[keep] <- suppressWarnings(stats::predict(fit, newdata = xx[keep]))
      } else if (identical(smooth_method, "spline")) {
        fit <- try(stats::smooth.spline(x = xx[keep], y = yy[keep]), silent = TRUE)
        if (!inherits(fit, "try-error")) pred[keep] <- suppressWarnings(stats::predict(fit, x = xx[keep])$y)
      }
      ok <- is.finite(pred)
      agg$value_smooth[ok] <- pred[ok]
    }
  }
  agg$value <- if (identical(smooth_method, "none")) agg$value_raw else agg$value_smooth
  agg
}

.trajectory_s2_bootstrap_ci <- function(pt,
                                                values,
                                                bins,
                                                n_bins_used,
                                                summary_fun = "mean",
                                                bootstrap_n = 0L,
                                                ci_level = 0.95,
                                                seed = 1L) {
  bootstrap_n <- suppressWarnings(as.integer(bootstrap_n))
  if (!is.finite(bootstrap_n) || bootstrap_n <= 0L) {
    return(data.frame(bin = seq_len(n_bins_used), value_se = NA_real_, value_ci_low = NA_real_,
                      value_ci_high = NA_real_, stringsAsFactors = FALSE))
  }
  ci_level <- suppressWarnings(as.numeric(ci_level)[1])
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) ci_level <- 0.95
  ss_fun <- if (identical(summary_fun, "median")) stats::median else base::mean
  alpha <- (1 - ci_level) / 2
  had_seed <- exists(".Random.seed", envir = .GlobalEnv,
    inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  seed <- suppressWarnings(as.integer(seed)[1])
  if (!is.finite(seed)) seed <- 1L
  set.seed(seed)
  rows <- lapply(seq_len(n_bins_used), function(bb) {
    idx <- which(is.finite(pt) & is.finite(values) & is.finite(bins) & bins == bb)
    if (!length(idx)) {
      return(data.frame(bin = bb, value_se = NA_real_, value_ci_low = NA_real_, value_ci_high = NA_real_))
    }
    vals <- as.numeric(values[idx])
    boot <- replicate(bootstrap_n, ss_fun(sample(vals, length(vals), replace = TRUE), na.rm = TRUE))
    data.frame(
      bin = as.integer(bb),
      value_se = as.numeric(stats::sd(boot, na.rm = TRUE)),
      value_ci_low = as.numeric(stats::quantile(boot, probs = alpha, na.rm = TRUE, names = FALSE)),
      value_ci_high = as.numeric(stats::quantile(boot, probs = 1 - alpha, na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.trajectory_s2_curve <- function(entity_id,
                                         lineage_id,
                                         pt,
                                         values,
                                         bins,
                                         n_bins_used,
                                         min_cells_per_bin,
                                         summary_fun,
                                         smooth_method,
                                         loess_span,
                                         bootstrap_n,
                                         ci_level,
                                         seed) {
  agg <- .trajectory_s2_aggregate_one(
    entity_id = entity_id,
    lineage_id = lineage_id,
    pt = pt,
    values = values,
    bins = bins,
    min_cells_per_bin = min_cells_per_bin,
    value_col = "value",
    summary_fun = summary_fun
  )
  if (is.null(agg) || !nrow(agg)) return(NULL)
  agg <- .trajectory_s2_apply_smoothing(agg, smooth_method = smooth_method, loess_span = loess_span)
  ci <- .trajectory_s2_bootstrap_ci(
    pt = pt,
    values = values,
    bins = bins,
    n_bins_used = n_bins_used,
    summary_fun = summary_fun,
    bootstrap_n = bootstrap_n,
    ci_level = ci_level,
    seed = seed
  )
  merge(agg, ci, by = "bin", all.x = TRUE, sort = FALSE)
}

.trajectory_s2_null_dynamic_scores <- function(entity_id,
                                                       lineage_id,
                                                       pt,
                                                       values,
                                                       bins,
                                                       min_cells_per_bin,
                                                       summary_fun,
                                                       smooth_method,
                                                       loess_span,
                                                       n_permutations) {
  n_permutations <- suppressWarnings(as.integer(n_permutations)[1])
  if (!is.finite(n_permutations) || n_permutations <= 0L) return(numeric(0))
  values <- suppressWarnings(as.numeric(values))
  ok <- is.finite(values) & is.finite(pt) & is.finite(bins)
  if (sum(ok) < min_cells_per_bin) return(rep(NA_real_, n_permutations))

  out <- rep(NA_real_, n_permutations)
  for (ii in seq_len(n_permutations)) {
    perm_values <- values
    perm_values[ok] <- sample(values[ok], length(values[ok]), replace = FALSE)
    agg <- .trajectory_s2_aggregate_one(
      entity_id = entity_id,
      lineage_id = lineage_id,
      pt = pt,
      values = perm_values,
      bins = bins,
      min_cells_per_bin = min_cells_per_bin,
      value_col = "value",
      summary_fun = summary_fun
    )
    if (is.null(agg) || !nrow(agg)) next
    agg <- .trajectory_s2_apply_smoothing(agg, smooth_method = smooth_method, loess_span = loess_span)
    mm <- .trajectory_s2_dynamic_metrics(agg)
    out[ii] <- suppressWarnings(as.numeric(mm$dynamic_score)[1])
  }
  out
}
