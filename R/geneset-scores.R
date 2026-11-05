.trajectory_geneset_empty_adjustment_table <- function() {
  data.frame(
    lineage_id = character(), feature_id = character(),
    depth_field = character(), status = character(),
    adjustment_method = character(), depth_coefficient = numeric(),
    depth_partial_r_squared = numeric(), raw_mean = numeric(),
    raw_sd = numeric(), adjusted_mean = numeric(), adjusted_sd = numeric(),
    standardized = logical(), n_cells = integer(), design_rank = integer(),
    design_columns = integer(), spline_df = integer(),
    stringsAsFactors = FALSE)
}

.trajectory_geneset_score_lineage <- function(matrix, set_members, cells,
                                               pseudotime, depth,
                                               adjust_depth, depth_field,
                                               score_method, lineage_id,
                                               label =
                                                 "trajectory_geneset_scores") {
  raw_matrix <- as.matrix(matrix[, cells, drop = FALSE])
  adjusted_matrix <- raw_matrix
  diagnostics <- .trajectory_geneset_empty_adjustment_table()
  if (adjust_depth) {
    adjusted <- .trajectory_depth_adjust_matrix(raw_matrix,
      pseudotime = pseudotime, depth = depth[cells], spline_df = 4L,
      standardize = TRUE,
      label = paste0(label, " lineage '", lineage_id, "'"))
    adjusted_matrix <- adjusted$values
    diagnostics <- adjusted$diagnostics
    diagnostics$lineage_id <- as.character(lineage_id)
    diagnostics$depth_field <- depth_field
    diagnostics$adjustment_method <- paste0(
      "expression~ns(pseudotime,df=4)+centered_log1p(depth);",
      "remove_depth_component;zscore")
    diagnostics <- diagnostics[, c("lineage_id", "feature_id",
      "depth_field", "status", "adjustment_method", "depth_coefficient",
      "depth_partial_r_squared", "raw_mean", "raw_sd", "adjusted_mean",
      "adjusted_sd", "standardized", "n_cells", "design_rank",
      "design_columns", "spline_df"), drop = FALSE]
  }

  rank_background <- NULL
  if (identical(score_method, "rank_mean")) {
    features <- unique(unlist(set_members, use.names = FALSE))
    rank_background <- raw_matrix[features[nzchar(features)], , drop = FALSE]
  }
  scores <- list()
  detected <- integer(length(set_members))
  fractions <- numeric(length(set_members))
  names(detected) <- names(fractions) <- names(set_members)
  rows <- vector("list", length(set_members))
  index <- 0L
  for (name in names(set_members)) {
    features <- set_members[[name]]
    raw <- raw_matrix[features, , drop = FALSE]
    effective <- if (adjust_depth)
      adjusted_matrix[features, , drop = FALSE] else raw
    detected[[name]] <- as.integer(sum(rowSums(raw > 0,
      na.rm = TRUE) > 0))
    fractions[[name]] <- if (nrow(raw)) detected[[name]] / nrow(raw) else 0
    score <- if (adjust_depth) colMeans(effective) else
      .trajectory_s2_score_values(effective, score_method = score_method,
        rank_bg = rank_background)
    scores[[name]] <- as.numeric(score)
    index <- index + 1L
    rows[[index]] <- data.frame(
      key = rep(paste(lineage_id, name, sep = "::"), length(cells)),
      cell_id = as.character(cells), value = as.numeric(score),
      lineage_id = rep(as.character(lineage_id), length(cells)),
      feature_set = rep(as.character(name), length(cells)),
      stringsAsFactors = FALSE)
  }
  score_table <- do.call(rbind, rows)
  rownames(score_table) <- NULL
  list(matrix = adjusted_matrix, scores = scores,
    n_features_detected = detected, detected_fraction = fractions,
    score_table = score_table, diagnostics = diagnostics)
}

.trajectory_s2_geneset_scores <- function(matrix, registry, effective_cells,
                                           matrix_id, pseudotime,
                                           depth = NULL,
                                           depth_artifact_id = "",
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
                                           min_set_size = 5L,
                                           max_set_size = 5000L,
                                           pseudotime_direction = "as_is") {
  .trajectory_assert(is.matrix(matrix) && !is.null(rownames(matrix)) &&
      !is.null(colnames(matrix)),
    "trajectory_geneset_scores: matrix is malformed.")
  .trajectory_assert(identical(colnames(matrix),
      as.character(effective_cells)), paste0(
        "trajectory_geneset_scores: matrix cell axis differs from the ",
        "authenticated effective axis."))
  .trajectory_kind(pseudotime, "cell_vector_family")
  .trajectory_assert(is.character(matrix_id) && length(matrix_id) == 1L &&
      nzchar(matrix_id), "trajectory_geneset_scores: matrix_id is invalid.")
  .trajectory_assert(is.character(namespace) && length(namespace) == 1L &&
      nzchar(namespace), "trajectory_geneset_scores: namespace is required.")

  adjust_depth <- isTRUE(adjust_depth)
  depth_field <- trimws(as.character(depth_field %||% "nCount_RNA")[1L])
  .trajectory_assert(!is.na(depth_field) && nzchar(depth_field),
    "trajectory_geneset_scores: depth_field is required.")
  if (adjust_depth) {
    depth <- as.numeric(depth[match(as.character(effective_cells),
      names(depth))])
    names(depth) <- as.character(effective_cells)
    .trajectory_assert(length(depth) == length(effective_cells) &&
        all(is.finite(depth)) && all(depth >= 0), paste0(
          "trajectory_geneset_scores: depth must be complete, finite, ",
          "non-negative, and cell-aligned."))
  } else {
    depth <- stats::setNames(rep(0, length(effective_cells)),
      as.character(effective_cells))
  }
  score_method <- as.character(score_method %||% "mean")[1L]
  .trajectory_assert(score_method %in% c("mean", "sum", "rank_mean"),
    paste0("trajectory_geneset_scores: score_method must be 'mean', 'sum', ",
      "or 'rank_mean'."))
  .trajectory_assert(!adjust_depth || identical(score_method, "mean"),
    paste0("trajectory_geneset_scores: adjust_depth=TRUE requires ",
      "score_method='mean'."))
  min_set_size <- suppressWarnings(as.integer(min_set_size)[1L])
  max_set_size <- suppressWarnings(as.integer(max_set_size)[1L])
  .trajectory_assert(is.finite(min_set_size) && min_set_size >= 1L &&
      is.finite(max_set_size) && max_set_size >= min_set_size,
    "trajectory_geneset_scores: invalid feature-set size bounds.")
  pseudotime_direction <- as.character(
    pseudotime_direction %||% "as_is")[1L]
  .trajectory_assert(pseudotime_direction %in% c("as_is", "invert"),
    "trajectory_geneset_scores: invalid pseudotime_direction.")

  universe <- as.character(rownames(matrix))
  .trajectory_assert(length(universe) && !anyDuplicated(universe),
    "trajectory_geneset_scores: matrix feature universe is invalid.")
  selected <- .trajectory_s2_registry_sets(registry, namespace,
    feature_sets, universe, min_set_size, max_set_size)
  source <- .trajectory_s2_explicit_source(pseudotime = pseudotime,
    pseudotime_key = pseudotime_key, membership = membership,
    lineage = lineage, cells = effective_cells,
    pseudotime_artifact_id = pseudotime_artifact_id,
    membership_artifact_id = membership_artifact_id)

  score_rows <- list()
  adjustment_rows <- list()
  summary_rows <- list()
  score_index <- adjustment_index <- summary_index <- 0L
  for (lineage_id in names(source$lineages)) {
    lineage_table <- source$lineages[[lineage_id]]
    pseudotime_value <- suppressWarnings(
      as.numeric(lineage_table$pseudotime))
    names(pseudotime_value) <- as.character(lineage_table$cell_id)
    pseudotime_value <- pseudotime_value[is.finite(pseudotime_value)]
    if (identical(pseudotime_direction, "invert"))
      pseudotime_value <- max(pseudotime_value) - pseudotime_value
    cells <- names(pseudotime_value)
    .trajectory_assert(length(cells) >= 10L,
      paste0("trajectory_geneset_scores: lineage '", lineage_id,
        "' has fewer than 10 finite pseudotime cells."))
    .trajectory_assert(stats::sd(pseudotime_value) > 0,
      paste0("trajectory_geneset_scores: lineage '", lineage_id,
        "' has zero pseudotime variance."))
    scored <- .trajectory_geneset_score_lineage(matrix,
      selected$members, cells, pseudotime_value, depth, adjust_depth,
      depth_field, score_method, lineage_id)
    score_index <- score_index + 1L
    score_rows[[score_index]] <- scored$score_table
    if (nrow(scored$diagnostics)) {
      adjustment_index <- adjustment_index + 1L
      adjustment_rows[[adjustment_index]] <- scored$diagnostics
    }
    summary_index <- summary_index + 1L
    summary_rows[[summary_index]] <- data.frame(
      lineage_id = as.character(lineage_id), namespace = namespace,
      pseudotime_source_kind = as.character(
        source$meta$pseudotime_source_kind),
      pseudotime_artifact_id = as.character(
        source$meta$pseudotime_artifact_id),
      pseudotime_key = as.character(source$meta$pseudotime_key),
      membership_artifact_id = as.character(
        source$meta$membership_artifact_id),
      depth_artifact_id = as.character(depth_artifact_id),
      n_pseudotime_cells = as.integer(source$meta$n_pseudotime_cells),
      n_membership_cells = as.integer(source$meta$n_membership_cells),
      n_cells_used = as.integer(length(cells)),
      n_feature_sets_input = as.integer(length(selected$names)),
      n_feature_sets_scored = as.integer(length(selected$members)),
      score_method = score_method,
      effective_score_method = if (adjust_depth) "adjusted_z_mean" else
        score_method,
      adjust_depth = adjust_depth, depth_field = depth_field,
      adjustment_spline_df = if (adjust_depth) 4L else NA_integer_,
      pseudotime_direction = pseudotime_direction,
      min_set_size = min_set_size, max_set_size = max_set_size,
      matrix_id = matrix_id, stringsAsFactors = FALSE)
  }
  score_table <- do.call(rbind, score_rows)
  rownames(score_table) <- NULL
  adjustment_table <- if (adjustment_index)
    do.call(rbind, adjustment_rows) else
      .trajectory_geneset_empty_adjustment_table()
  rownames(adjustment_table) <- NULL
  fit_summary <- do.call(rbind, summary_rows)
  rownames(fit_summary) <- NULL
  score_meta <- list(matrix_id = matrix_id, namespace = namespace,
    depth_adjusted = adjust_depth, depth_field = depth_field,
    depth_artifact_id = as.character(depth_artifact_id),
    adjustment_spline_df = if (adjust_depth) 4L else NA_integer_,
    pseudotime_artifact_id = as.character(
      source$meta$pseudotime_artifact_id),
    pseudotime_key = as.character(source$meta$pseudotime_key),
    membership_artifact_id = as.character(
      source$meta$membership_artifact_id),
    score_key_format = "lineage_id::feature_set",
    effective_score_method = if (adjust_depth) "adjusted_z_mean" else
      score_method)
  list(outputs = list(
    adjusted_scores = scworkbench::wb_prim_cell_vector_family(
      key = score_table$key, cell_id = score_table$cell_id,
      value = score_table$value, value_type = "continuous",
      meta = score_meta),
    adjustment_diagnostics = scworkbench::wb_prim_obs_table(
      table = adjustment_table,
      table_type = "trajectory_geneset_depth_adjustment",
      entity_kind = "feature", meta = score_meta),
    fit_summary = scworkbench::wb_prim_obs_table(table = fit_summary,
      table_type = "trajectory_geneset_scores_summary",
      entity_kind = "lineage", meta = score_meta)),
    params = c(score_meta, list(
      n_feature_sets_input = as.integer(length(selected$names)),
      n_feature_sets_scored = as.integer(length(selected$members)),
      n_lineages = as.integer(length(source$lineages)),
      score_method = score_method, adjust_depth = adjust_depth,
      pseudotime_direction = pseudotime_direction,
      min_set_size = min_set_size, max_set_size = max_set_size)))
}
