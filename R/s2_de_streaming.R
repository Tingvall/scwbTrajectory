# Consume the authenticated feature axis once, retaining only final DE rows.
# The reader's chunk discipline requires every axis partition to be requested,
# including partitions with no selected features.
.trajectory_s2_de_stream <- function(reader, all_features, selected_features,
                                     source, effective_cells, matrix_id,
                                     feature_set = NULL, features = NULL,
                                     test = "spearman", lineage = NULL,
                                     min_pct_nonzero = 0, min_sd = 0,
                                     p_adjust_method = "BH",
                                     chunk_size = 128L,
                                     observe_block = NULL) {
  .trajectory_assert(identical(as.character(test)[1L], "spearman"),
    "trajectory_de supports test='spearman' only.")
  .trajectory_assert(as.character(p_adjust_method)[1L] %in%
    stats::p.adjust.methods, "trajectory_de: invalid p_adjust_method.")
  .trajectory_assert(length(chunk_size) == 1L &&
    is.finite(as.numeric(chunk_size)) && as.integer(chunk_size) >= 1L,
    "trajectory_de: chunk_size must be positive.")
  all_features <- as.character(all_features)
  selected_features <- as.character(selected_features)
  effective_cells <- as.character(effective_cells)
  .trajectory_assert(length(all_features) > 0L &&
      !anyDuplicated(all_features) &&
      length(selected_features) > 0L &&
      !anyDuplicated(selected_features) &&
      all(selected_features %in% all_features) &&
      length(effective_cells) > 0L && !anyDuplicated(effective_cells),
    "trajectory_de: authenticated matrix axes are invalid.")
  all_axis <- source$axis[source$axis$cell_id %in% effective_cells &
    is.finite(source$axis$pseudotime), , drop = FALSE]
  n_lineages_input <- length(unique(as.character(all_axis$lineage_id)))
  lineages <- .trajectory_s2_model_lineages(source, lineage,
    effective_cells)
  if (!is.null(lineage)) {
    requested <- unique(as.character(unlist(lineage, use.names = FALSE)))
    lineages <- lineages[requested]
  }
  states <- lapply(lineages, function(item) {
    pt <- as.numeric(item$pseudotime)
    names(pt) <- as.character(item$cell_id)
    pt <- pt[is.finite(pt)]
    pt <- pt[!duplicated(names(pt))]
    if (length(pt) < 10L || !is.finite(stats::sd(pt)) ||
        stats::sd(pt) <= 0) return(NULL)
    index <- match(names(pt), effective_cells)
    .trajectory_assert(!anyNA(index),
      "trajectory_de: fitted cells exceed the authenticated matrix axis.")
    list(item = item, pt = pt, cell_index = index,
      rows = vector("list", length(selected_features)),
      filtered_cells = 0L, filtered_nonzero = 0L, filtered_sd = 0L)
  })
  active <- names(states)[!vapply(states, is.null, logical(1))]
  .trajectory_assert(length(active) > 0L,
    "trajectory_de: no lineage produced testable cells/features.")
  chunks <- scworkbench::wb_partition_axis(all_features,
    as.integer(chunk_size))
  for (chunk in chunks) {
    block <- reader(features = chunk, cells = effective_cells)$values
    .trajectory_assert(identical(as.character(rownames(block)),
        as.character(chunk)) &&
        identical(as.character(colnames(block)), effective_cells),
      "trajectory_de: typed feature chunk has an unexpected axis.")
    keep <- which(chunk %in% selected_features)
    if (!length(keep)) next
    dense <- as.matrix(block[keep, , drop = FALSE])
    .trajectory_assert(identical(dim(dense),
        c(as.integer(length(keep)), as.integer(length(effective_cells)))),
      "trajectory_de: typed feature chunk has an unexpected shape.")
    if (is.function(observe_block))
      observe_block(nrow(dense), ncol(dense))
    positions <- match(as.character(chunk[keep]), selected_features)
    for (lin in active) {
      state <- states[[lin]]
      for (k in seq_along(keep)) {
        one <- .trajectory_s2_de_one_feature(
          as.numeric(dense[k, state$cell_index]), state$pt,
          lin, as.character(chunk[keep[[k]]]),
          min_pct_nonzero, min_sd)
        if (identical(one$status, "filtered_cells"))
          state$filtered_cells <- state$filtered_cells + 1L
        else if (identical(one$status, "filtered_nonzero"))
          state$filtered_nonzero <- state$filtered_nonzero + 1L
        else if (identical(one$status, "filtered_sd"))
          state$filtered_sd <- state$filtered_sd + 1L
        else if (identical(one$status, "tested"))
          state$rows[[positions[[k]]]] <- one$row
      }
      states[[lin]] <- state
    }
  }
  rows <- list(); summaries <- list(); ri <- 0L
  for (lin in active) {
    state <- states[[lin]]
    one_rows <- Filter(Negate(is.null), state$rows)
    tbl <- if (length(one_rows)) do.call(rbind, one_rows) else data.frame(
      lineage_id = character(), feature_id = character(), rho = numeric(),
      abs_rho = numeric(), direction = character(), slope = numeric(),
      p_val = numeric(), p_val_adj = numeric(), pct_nonzero = numeric(),
      sd_expr = numeric(), mean_expr = numeric(), n_used = integer(),
      n_dropped_na = integer(), stringsAsFactors = FALSE)
    ri <- ri + 1L; rows[[ri]] <- tbl
    summaries[[ri]] <- data.frame(lineage_id = as.character(lin),
      n_cells_input = as.integer(length(state$item$cell_id)),
      n_cells_used = as.integer(length(state$pt)),
      n_features_input = as.integer(length(selected_features)),
      n_features_tested = as.integer(nrow(tbl)),
      n_filtered_min_cells = as.integer(state$filtered_cells),
      n_filtered_min_pct_nonzero = as.integer(state$filtered_nonzero),
      n_filtered_min_sd = as.integer(state$filtered_sd),
      stringsAsFactors = FALSE)
  }
  .trajectory_s2_de_finish(rows, summaries, matrix_id, source,
    feature_set, min_pct_nonzero, min_sd, p_adjust_method,
    n_lineages_input)
}
