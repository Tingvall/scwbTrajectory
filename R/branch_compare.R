.trajectory_branch_cells <- function(membership, branch_a, branch_b,
                                     effective_cells) {
  .trajectory_kind(membership, "relation_family")
  table <- scworkbench::wb_primitive_to_table(membership,
    include_context = FALSE)
  .trajectory_assert(nrow(table) > 0L &&
      all(c("key", "cell_id") %in% names(table)),
    "trajectory_branch_compare: branch_membership is empty or malformed.")
  table <- table[as.character(table$key) %in% c(branch_a, branch_b) &
    as.character(table$cell_id) %in% effective_cells, , drop = FALSE]
  .trajectory_assert(nrow(table) > 0L,
    "trajectory_branch_compare: no requested branch rows remain.")
  table$key <- as.character(table$key)
  table$cell_id <- as.character(table$cell_id)
  table$weight_use <- if ("weight" %in% names(table))
    suppressWarnings(as.numeric(table$weight)) else 1
  table$weight_use[!is.finite(table$weight_use)] <- 1
  aggregate <- stats::aggregate(weight_use ~ cell_id + key, table, max)
  names(aggregate)[names(aggregate) == "weight_use"] <- "weight"
  a <- aggregate[aggregate$key == branch_a, c("cell_id", "weight"),
    drop = FALSE]
  b <- aggregate[aggregate$key == branch_b, c("cell_id", "weight"),
    drop = FALSE]
  names(a)[2L] <- "weight_a"; names(b)[2L] <- "weight_b"
  merged <- merge(a, b, by = "cell_id", all = TRUE)
  in_a <- !is.na(merged$weight_a); in_b <- !is.na(merged$weight_b)
  assignment <- rep(NA_character_, nrow(merged))
  assignment[in_a & !in_b] <- "a"; assignment[!in_a & in_b] <- "b"
  both <- in_a & in_b
  assignment[both & merged$weight_a > merged$weight_b] <- "a"
  assignment[both & merged$weight_b > merged$weight_a] <- "b"
  list(cells_a = unique(merged$cell_id[assignment == "a" &
      !is.na(assignment)]),
    cells_b = unique(merged$cell_id[assignment == "b" &
      !is.na(assignment)]),
    n_overlap_dropped = as.integer(sum(both & is.na(assignment))))
}

.trajectory_branch_pt_summary <- function(pseudotime, branch_a, branch_b,
                                          cells_a, cells_b) {
  .trajectory_kind(pseudotime, "cell_vector_family")
  table <- scworkbench::wb_primitive_to_table(pseudotime,
    include_context = FALSE)
  .trajectory_assert(all(c("key", "cell_id", "value") %in% names(table)),
    "trajectory_branch_compare: pseudotime family is malformed.")
  value <- suppressWarnings(as.numeric(table$value))
  a <- value[as.character(table$key) == branch_a &
    as.character(table$cell_id) %in% cells_a]
  b <- value[as.character(table$key) == branch_b &
    as.character(table$cell_id) %in% cells_b]
  list(mean_pt_a = mean(a[is.finite(a)]),
    mean_pt_b = mean(b[is.finite(b)]),
    median_pt_a = stats::median(a[is.finite(a)]),
    median_pt_b = stats::median(b[is.finite(b)]))
}

.trajectory_de <- function(matrix, cells_a, cells_b) {
  output <- vector("list", nrow(matrix))
  for (i in seq_len(nrow(matrix))) {
    a <- as.numeric(matrix[i, cells_a, drop = TRUE])
    b <- as.numeric(matrix[i, cells_b, drop = TRUE])
    p <- tryCatch(suppressWarnings(stats::wilcox.test(a, b,
      alternative = "two.sided", exact = FALSE)$p.value),
      error = function(error) 1)
    output[[i]] <- data.frame(feature_id = rownames(matrix)[i],
      pct_a = mean(a > 0), pct_b = mean(b > 0), avg_expr_a = mean(a),
      avg_expr_b = mean(b),
      avg_log1p_fc = log1p(mean(a)) - log1p(mean(b)),
      p_val = as.numeric(p), stringsAsFactors = FALSE)
  }
  do.call(rbind, output)
}

.run_trajectory_branch_compare <- function(readers, writers, params, context) {
  pseudotime <- readers$pseudotime()
  membership <- readers$branch_membership()
  branch_a <- as.character(params$branch_a %||% "")[1L]
  branch_b <- as.character(params$branch_b %||% "")[1L]
  .trajectory_assert(nzchar(branch_a) && nzchar(branch_b),
    "trajectory_branch_compare: branch_a and branch_b are required.")
  .trajectory_assert(!identical(branch_a, branch_b),
    "trajectory_branch_compare: branch_a and branch_b must differ.")
  all_features <- as.character(
    context$typed_slice_reads$matrix_chunks$features)
  features <- .trajectory_feature_sets(readers, params, all_features,
    required = FALSE, prefix = "trajectory_branch_compare")
  cells <- as.character(context$effective_cells)
  resolved <- .trajectory_branch_cells(membership, branch_a, branch_b, cells)
  .trajectory_assert(length(resolved$cells_a) >= 2L,
    "trajectory_branch_compare: branch_a has fewer than 2 assigned cells.")
  .trajectory_assert(length(resolved$cells_b) >= 2L,
    "trajectory_branch_compare: branch_b has fewer than 2 assigned cells.")
  matrix <- .trajectory_chunked_matrix(readers$matrix_chunks, all_features,
    features, cells, unique(c(resolved$cells_a, resolved$cells_b)),
    params$chunk_size %||% 128L)
  de <- .trajectory_de(matrix, resolved$cells_a, resolved$cells_b)
  min_pct <- as.numeric(params$min_pct %||% 0)
  if (is.finite(min_pct) && min_pct > 0)
    de <- de[pmax(de$pct_a, de$pct_b) >= min_pct, , drop = FALSE]
  de$p_val_adj <- if (nrow(de)) stats::p.adjust(de$p_val,
    method = as.character(params$p_adjust_method %||% "BH")[1L]) else
    numeric()
  de$branch_a <- branch_a; de$branch_b <- branch_b
  de$n_cells_a <- as.integer(length(resolved$cells_a))
  de$n_cells_b <- as.integer(length(resolved$cells_b))
  de <- de[order(de$p_val_adj, -de$avg_log1p_fc), , drop = FALSE]
  rownames(de) <- NULL
  pt <- .trajectory_branch_pt_summary(pseudotime, branch_a, branch_b,
    resolved$cells_a, resolved$cells_b)
  summary <- data.frame(branch_a = branch_a, branch_b = branch_b,
    n_cells_a = as.integer(length(resolved$cells_a)),
    n_cells_b = as.integer(length(resolved$cells_b)),
    n_overlap_dropped = resolved$n_overlap_dropped,
    n_features_input = as.integer(length(features)),
    n_features_tested = as.integer(nrow(de)),
    mean_pt_a = as.numeric(pt$mean_pt_a), mean_pt_b = as.numeric(pt$mean_pt_b),
    median_pt_a = as.numeric(pt$median_pt_a),
    median_pt_b = as.numeric(pt$median_pt_b), stringsAsFactors = FALSE)
  matrix_id <- as.character(context$input_ids$matrix)[1L]
  writers$write_terminal_primitive("branch_de_table",
    scworkbench::wb_prim_feature_table(de,
      table_type = "trajectory_branch_compare",
      meta = list(matrix_id = matrix_id, branch_a = branch_a,
        branch_b = branch_b)))
  writers$write_terminal_primitive("branch_summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "trajectory_branch_compare_summary",
      entity_kind = "branch_pair", meta = list(matrix_id = matrix_id,
        branch_a = branch_a, branch_b = branch_b)))
  invisible(list(status = "ok"))
}
