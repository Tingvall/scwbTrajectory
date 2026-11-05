.trajectory_adjustment_spline_df <- function(value = 4L,
                                             label = "depth adjustment") {
  raw <- suppressWarnings(as.numeric(value)[1L])
  .trajectory_assert(length(raw) == 1L && is.finite(raw) &&
      raw == floor(raw) && raw >= 2L && raw <= 10L,
    paste0(label, ": adjustment_spline_df must be an integer from 2 to 10."))
  as.integer(raw)
}

.trajectory_depth_vector <- function(value, cells, label = "depth") {
  .trajectory_kind(value, "cell_vector")
  ids <- as.character(value$cell_id)
  values <- suppressWarnings(as.numeric(value$value))
  .trajectory_assert(length(ids) == length(values) && !anyNA(ids) &&
      all(nzchar(ids)) && !anyDuplicated(ids),
    paste0(label, ": depth vector has an invalid cell axis."))
  index <- match(as.character(cells), ids)
  .trajectory_assert(!anyNA(index),
    paste0(label, ": depth vector does not cover the effective cell axis."))
  values <- values[index]
  .trajectory_assert(all(is.finite(values)) && all(values >= 0),
    paste0(label, ": depth values must be finite and non-negative."))
  stats::setNames(values, as.character(cells))
}

.trajectory_depth_design <- function(pseudotime, depth, spline_df = 4L,
                                     label = "depth adjustment") {
  pseudotime <- as.numeric(pseudotime)
  depth <- as.numeric(depth)
  spline_df <- as.integer(spline_df)[1L]
  .trajectory_assert(length(pseudotime) == length(depth) &&
      length(depth) >= spline_df + 2L && all(is.finite(pseudotime)) &&
      all(is.finite(depth)) && all(depth >= 0) &&
      length(unique(pseudotime)) > spline_df,
    paste0(label, ": insufficient finite pseudotime/depth support."))
  log_depth <- log1p(depth)
  centered_depth <- log_depth - mean(log_depth)
  basis <- splines::ns(pseudotime, df = spline_df)
  design <- cbind(`(Intercept)` = 1, basis,
    centered_log1p_depth = centered_depth)
  fit <- qr(design)
  if (stats::sd(centered_depth) == 0) {
    design <- design[, -ncol(design), drop = FALSE]
    fit <- qr(design)
    .trajectory_assert(fit$rank == ncol(design),
      paste0(label, ": pseudotime spline design is rank deficient."))
    return(list(design = design, fit = fit, reduced = fit,
      centered_depth = centered_depth, constant_depth = TRUE,
      spline_df = spline_df))
  }
  .trajectory_assert(fit$rank == ncol(design),
    paste0(label, ": pseudotime and depth design is rank deficient."))
  reduced_design <- design[, -ncol(design), drop = FALSE]
  reduced <- qr(reduced_design)
  .trajectory_assert(reduced$rank == ncol(reduced_design),
    paste0(label, ": pseudotime spline design is rank deficient."))
  list(design = design, fit = fit, reduced = reduced,
    centered_depth = centered_depth, constant_depth = FALSE,
    spline_df = spline_df)
}

.trajectory_depth_adjust_matrix <- function(matrix, pseudotime, depth,
                                            spline_df = 4L,
                                            standardize = FALSE,
                                            label = "depth adjustment") {
  matrix <- as.matrix(matrix)
  storage.mode(matrix) <- "double"
  .trajectory_assert(ncol(matrix) == length(pseudotime) &&
      ncol(matrix) == length(depth) && all(is.finite(matrix)),
    paste0(label, ": expression values must be finite and cell-aligned."))
  design <- .trajectory_depth_design(pseudotime, depth, spline_df, label)
  if (design$constant_depth) {
    beta <- rep(0, nrow(matrix))
    adjusted <- matrix
    full_residual <- t(qr.resid(design$fit, t(matrix)))
    reduced_residual <- full_residual
  } else {
    coefficients <- qr.coef(design$fit, t(matrix))
    beta <- as.numeric(coefficients[nrow(coefficients), ])
    adjusted <- matrix - tcrossprod(beta, design$centered_depth)
    full_residual <- t(qr.resid(design$fit, t(matrix)))
    reduced_residual <- t(qr.resid(design$reduced, t(matrix)))
  }
  raw_mean <- rowMeans(matrix)
  raw_sd <- apply(matrix, 1L, stats::sd)
  adjusted_mean <- rowMeans(adjusted)
  adjusted_sd <- apply(adjusted, 1L, stats::sd)
  zero_variance <- !is.finite(adjusted_sd) | adjusted_sd <=
    sqrt(.Machine$double.eps)
  values <- adjusted
  if (isTRUE(standardize)) {
    values <- sweep(values, 1L, adjusted_mean, "-")
    good <- !zero_variance
    values[good, ] <- sweep(values[good, , drop = FALSE], 1L,
      adjusted_sd[good], "/")
    values[zero_variance, ] <- 0
  }
  full_sse <- rowSums(full_residual^2)
  reduced_sse <- rowSums(reduced_residual^2)
  partial_r_squared <- ifelse(reduced_sse > 0,
    pmax(0, pmin(1, (reduced_sse - full_sse) / reduced_sse)), 0)
  diagnostics <- data.frame(
    feature_id = as.character(rownames(matrix)),
    status = ifelse(zero_variance, "zero_adjusted_variance",
      if (design$constant_depth) "constant_depth" else "ok"),
    depth_coefficient = as.numeric(beta),
    depth_partial_r_squared = as.numeric(partial_r_squared),
    raw_mean = as.numeric(raw_mean), raw_sd = as.numeric(raw_sd),
    adjusted_mean = as.numeric(adjusted_mean),
    adjusted_sd = as.numeric(adjusted_sd),
    standardized = isTRUE(standardize),
    n_cells = as.integer(ncol(matrix)),
    design_rank = as.integer(design$fit$rank),
    design_columns = as.integer(ncol(design$design)),
    spline_df = as.integer(design$spline_df),
    stringsAsFactors = FALSE)
  list(values = values, adjusted = adjusted, diagnostics = diagnostics,
    constant_depth = design$constant_depth)
}

.trajectory_depth_adjust_vector <- function(values, pseudotime, depth,
                                            spline_df = 4L,
                                            label = "depth adjustment") {
  matrix <- matrix(as.numeric(values), nrow = 1L,
    dimnames = list("signal", NULL))
  result <- .trajectory_depth_adjust_matrix(matrix, pseudotime, depth,
    spline_df = spline_df, standardize = FALSE, label = label)
  list(values = as.numeric(result$values[1L, ]),
    diagnostics = result$diagnostics[1L, , drop = FALSE])
}
