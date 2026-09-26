# Monocle3 method implementation owned by scwbTrajectory. All bundle reads and
# writes are supplied by the authenticated operation-pack boundary.

.trajectory_monocle3_require_pkgs <- function() {
  .trajectory_assert(requireNamespace("monocle3", quietly = TRUE), "Monocle3 is required.")
  .trajectory_assert(requireNamespace("Matrix", quietly = TRUE), "Matrix is required.")
  .trajectory_assert(requireNamespace("SingleCellExperiment", quietly = TRUE), "SingleCellExperiment is required.")
  .trajectory_assert(requireNamespace("igraph", quietly = TRUE), "igraph is required.")
}

.trajectory_monocle3_as_int_vec <- function(x, allow_null = TRUE, name = "dims") {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop(name, " is required.", call. = FALSE)
  }
  if (!length(x)) return(integer(0))
  xi <- suppressWarnings(as.integer(x))
  .trajectory_assert(!anyNA(xi), paste0(name, " must be coercible to integer(s)."))
  xi
}

.trajectory_monocle3_dummy_cds <- function(cell_ids) {
  n <- length(cell_ids)
  .trajectory_assert(n > 0, "monocle3_fit requires at least one cell.")

  # Minimal non-zero assay to keep monocle3 internals happy while we inject a
  # scworkbench latent space into reducedDims(cds)$UMAP.
  mat <- Matrix::Matrix(1, nrow = 1L, ncol = n, sparse = TRUE)
  rownames(mat) <- "dummy_feature"
  colnames(mat) <- as.character(cell_ids)

  cell_md <- data.frame(row.names = as.character(cell_ids))
  gene_md <- data.frame(gene_short_name = "dummy_feature", row.names = "dummy_feature")
  monocle3::new_cell_data_set(
    expression_data = mat,
    cell_metadata = cell_md,
    gene_metadata = gene_md
  )
}

.trajectory_monocle3_lineage_paths <- function(cds) {
  pg <- monocle3::principal_graph(cds)[["UMAP"]] %||% NULL
  .trajectory_assert(!is.null(pg), "monocle3_fit: principal_graph(cds)[['UMAP']] is missing after learn_graph().")

  aux <- monocle3::principal_graph_aux(cds)[["UMAP"]] %||% NULL
  .trajectory_assert(!is.null(aux), "monocle3_fit: principal_graph_aux(cds)[['UMAP']] is missing after order_cells().")

  root_nodes <- as.character(aux$root_pr_nodes %||% character(0))
  root_nodes <- root_nodes[nzchar(root_nodes)]
  root_nodes <- unique(root_nodes)
  .trajectory_assert(length(root_nodes) >= 1L,
             "monocle3_fit requires at least one root principal node after order_cells().")

  vnames <- as.character(igraph::V(pg)$name %||% character(0))
  .trajectory_assert(length(vnames) > 0, "monocle3_fit: principal graph vertices are unnamed.")
  .trajectory_assert(all(root_nodes %in% vnames),
             "monocle3_fit: one or more root principal nodes not found in principal graph.")

  # V1 semantics assume a tree (no reconverging branches).
  is_conn <- isTRUE(igraph::is_connected(pg))
  is_tree <- isTRUE(igraph::is_tree(pg, mode = "all"))
  .trajectory_assert(is_conn && is_tree,
             "monocle3_fit v1 requires a connected tree principal graph (use close_loop=FALSE).")

  deg <- igraph::degree(pg)
  leaves_all <- as.character(names(deg)[deg == 1L])
  leaves <- setdiff(leaves_all, root_nodes)
  used_root_leaves <- FALSE
  if (!length(leaves)) {
    # Multi-root cohorts can mark all terminal leaves as roots. In that case,
    # keep leaf keys and allow root==leaf lineages so branch-membership remains
    # representable without changing the primitive contract.
    leaves <- leaves_all
    used_root_leaves <- TRUE
  }
  .trajectory_assert(length(leaves) >= 1L,
             "monocle3_fit: no terminal leaf nodes available for lineage derivation.")

  lineage_paths <- vector("list", length(leaves))
  lineage_roots <- character(length(leaves))
  lineage_leaves <- as.character(leaves)
  for (i in seq_along(leaves)) {
    leaf <- leaves[[i]]
    # Multi-root support: assign each leaf to the closest root in graph steps
    # (deterministic tie-break by lexicographic root node id).
    dd <- as.numeric(igraph::distances(pg, v = root_nodes, to = leaf, mode = "all"))
    .trajectory_assert(length(dd) == length(root_nodes), "monocle3_fit: root distance computation failed.")
    .trajectory_assert(any(is.finite(dd)), paste0("monocle3_fit: leaf '", leaf, "' is unreachable from all roots."))
    best <- which(is.finite(dd) & dd == min(dd[is.finite(dd)]))
    root_choice <- sort(as.character(root_nodes[best]))[[1]]
    lineage_roots[[i]] <- root_choice
    sp <- igraph::shortest_paths(pg, from = root_choice, to = leaf, output = "vpath")$vpath[[1]]
    nodes <- as.character(igraph::V(pg)[sp]$name)
    min_nodes <- if (isTRUE(used_root_leaves)) 1L else 2L
    .trajectory_assert(length(nodes) >= min_nodes,
               paste0("monocle3_fit: invalid root-to-leaf path for leaf '", leaf, "'."))
    lineage_paths[[i]] <- nodes
  }
  names(lineage_paths) <- paste0("root:", lineage_roots, "::lineage:", lineage_leaves)

  list(
    pg = pg,
    root_nodes = root_nodes,
    leaves = leaves,
    used_root_leaves = used_root_leaves,
    lineage_roots = lineage_roots,
    lineage_leaves = lineage_leaves,
    lineage_paths = lineage_paths
  )
}

.trajectory_monocle3_closest_vertex_names <- function(cds, pg) {
  aux <- monocle3::principal_graph_aux(cds)[["UMAP"]] %||% NULL
  cv <- aux$pr_graph_cell_proj_closest_vertex %||% NULL
  .trajectory_assert(!is.null(cv), "monocle3_fit: missing pr_graph_cell_proj_closest_vertex in principal_graph_aux.")

  if (is.matrix(cv) || is.array(cv)) {
    cell_ids <- rownames(cv) %||% character(0)
    vals <- cv[, 1, drop = TRUE]
  } else {
    vals <- cv
    cell_ids <- names(vals) %||% character(0)
  }
  cell_ids <- as.character(cell_ids)
  .trajectory_assert(length(cell_ids) == length(vals), "monocle3_fit: closest-vertex mapping malformed.")

  vnames <- as.character(igraph::V(pg)$name %||% character(0))
  .trajectory_assert(length(vnames) > 0, "monocle3_fit: principal graph vertex names missing.")

  if (is.numeric(vals) || is.integer(vals)) {
    idx <- suppressWarnings(as.integer(vals))
    ok <- !is.na(idx) & idx >= 1L & idx <= length(vnames)
    out <- rep(NA_character_, length(idx))
    out[ok] <- vnames[idx[ok]]
  } else {
    out <- as.character(vals)
    out[!(out %in% vnames)] <- NA_character_
  }

  names(out) <- cell_ids
  out
}

.trajectory_monocle3_branch_membership_prim <- function(cds, lineage_info) {
  pg <- lineage_info$pg
  lineage_paths <- lineage_info$lineage_paths
  closest_vertex <- .trajectory_monocle3_closest_vertex_names(cds, pg)
  cell_ids <- names(closest_vertex)
  .trajectory_assert(length(cell_ids) > 0, "monocle3_fit: no cells found in closest-vertex mapping.")

  rows <- vector("list", length(lineage_paths))
  keep_n <- 0L
  for (i in seq_along(lineage_paths)) {
    key <- names(lineage_paths)[i]
    path_nodes <- lineage_paths[[i]]
    hit <- !is.na(closest_vertex) & (closest_vertex %in% path_nodes)
    if (!any(hit)) next
    keep_n <- keep_n + 1L
    rows[[keep_n]] <- data.frame(
      key = rep(key, sum(hit)),
      cell_id = as.character(cell_ids[hit]),
      weight = rep(1, sum(hit)),
      stringsAsFactors = FALSE
    )
  }
  .trajectory_assert(keep_n > 0L, "monocle3_fit: derived zero branch-membership relations.")

  rel <- do.call(rbind, rows[seq_len(keep_n)])
  rownames(rel) <- NULL
  scworkbench::wb_prim_relation_family(
    key = rel$key,
    cell_id = rel$cell_id,
    weight = rel$weight,
    meta = list(
      key_type = "lineage",
      relation_semantics = "cell_to_lineage_membership",
      weight_semantics = "closest_vertex_path_membership_binary",
      normalized_per_cell = FALSE
    )
  )
}

.trajectory_monocle3_lineage_curves_2d_prim <- function(cds,
                                                lineage_info,
                                                plot_xy,
                                                approx_points = 200L) {
  .trajectory_assert(is.matrix(plot_xy) && ncol(plot_xy) == 2L, "plot_xy must be a 2-column matrix.")
  .trajectory_assert(!is.null(rownames(plot_xy)), "plot_xy must have rownames (cell IDs).")

  closest_vertex <- .trajectory_monocle3_closest_vertex_names(cds, lineage_info$pg)
  cc <- closest_vertex[rownames(plot_xy)]
  ok_cells <- !is.na(cc)
  .trajectory_assert(any(ok_cells), "monocle3_project: no cells could be mapped to principal graph vertices.")

  plot_xy <- plot_xy[ok_cells, , drop = FALSE]
  cc <- as.character(cc[ok_cells])

  rows <- list()
  row_i <- 0L
  for (k in names(lineage_info$lineage_paths)) {
    nodes <- lineage_info$lineage_paths[[k]]
    xs <- numeric(0)
    ys <- numeric(0)
    pts <- integer(0)
    for (j in seq_along(nodes)) {
      hit <- which(cc == nodes[[j]])
      if (!length(hit)) next
      row_i <- row_i + 1L
      rows[[row_i]] <- data.frame(
        key = as.character(k),
        point = as.integer(length(pts) + 1L),
        x = mean(plot_xy[hit, 1], na.rm = TRUE),
        y = mean(plot_xy[hit, 2], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      pts <- c(pts, j)
    }
  }
  .trajectory_assert(length(rows) > 0L, "monocle3_project: failed to derive any projected lineage curves.")
  df <- do.call(rbind, rows)
  rownames(df) <- NULL

  prim <- scworkbench::wb_prim_path_2d_family(key = df$key, point = df$point, x = df$x, y = df$y)
  prim <- .trajectory_monocle3_compact(prim, max_points_per_lineage = approx_points)
  prim
}

.trajectory_monocle3_compact <- function(curves, max_points_per_lineage) {
  nmax <- suppressWarnings(as.integer(max_points_per_lineage)[1L])
  if (is.na(nmax) || nmax <= 0L) return(curves)
  key <- as.character(curves$key)
  rows <- lapply(unique(key), function(k) {
    index <- which(key == k)
    index <- index[order(as.integer(curves$point[index]), seq_along(index))]
    if (length(index) > nmax)
      index <- index[unique(round(seq.int(1L, length(index),
        length.out = nmax)))]
    data.frame(key = rep(k, length(index)), point = seq_along(index),
      x = as.numeric(curves$x[index]), y = as.numeric(curves$y[index]))
  })
  rows <- do.call(rbind, rows)
  rownames(rows) <- NULL
  scworkbench::wb_prim_path_2d_family(rows$key, rows$point,
    rows$x, rows$y)
}

.trajectory_monocle3_fit_input <- function(space, effective_cells, params) {
  coordinates <- .trajectory_space_matrix(space)
  cells <- as.character(effective_cells)
  .trajectory_assert(length(cells) >= 10L && !anyDuplicated(cells) &&
      all(cells %in% rownames(coordinates)),
    "monocle3_fit requires at least ten unique authenticated fit cells.")
  dims <- .trajectory_monocle3_as_int_vec(params$dims %||% NULL)
  if (is.null(dims)) dims <- seq_len(ncol(coordinates))
  .trajectory_assert(length(dims) >= 2L && !anyDuplicated(dims) &&
      all(dims >= 1L & dims <= ncol(coordinates)),
    "monocle3_fit requires at least two unique available dimensions.")
  has_cells <- length(params$root_cells %||% character()) > 0L
  has_nodes <- length(params$root_pr_nodes %||% character()) > 0L
  .trajectory_assert(xor(has_cells, has_nodes),
    "monocle3_fit requires exactly one of root_cells or root_pr_nodes.")
  roots <- if (has_cells) intersect(as.character(params$root_cells), cells)
    else NULL
  if (has_cells) .trajectory_assert(length(roots) > 0L,
    "monocle3_fit: root_cells has no overlap with fitted cells.")
  values <- coordinates[match(cells, rownames(coordinates)), dims,
    drop = FALSE]
  list(cells = cells, values = values, root_cells = roots,
    root_pr_nodes = if (has_nodes) as.character(params$root_pr_nodes) else
      NULL)
}

.trajectory_monocle3_fit_cds <- function(fit, params) {
  .trajectory_monocle3_require_pkgs()
  cds <- .trajectory_monocle3_dummy_cds(fit$cells)
  SingleCellExperiment::reducedDims(cds)$UMAP <- fit$values
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP",
    k = as.integer(params$k %||% 20L),
    cluster_method = as.character(params$cluster_method %||% "leiden"),
    num_iter = as.integer(params$num_iter %||% 2L),
    partition_qval = as.numeric(params$partition_qval %||% 0.05),
    weight = isTRUE(params$weight %||% FALSE),
    resolution = params$resolution %||% NULL,
    random_seed = as.integer(params$random_seed %||% 42L),
    verbose = isTRUE(params$verbose %||% FALSE), nn_control = list())
  cds <- monocle3::learn_graph(cds,
    use_partition = isTRUE(params$use_partition %||% FALSE),
    close_loop = isTRUE(params$close_loop %||% FALSE),
    learn_graph_control = NULL,
    verbose = isTRUE(params$verbose %||% FALSE))
  monocle3::order_cells(cds, reduction_method = "UMAP",
    root_pr_nodes = fit$root_pr_nodes, root_cells = fit$root_cells,
    verbose = isTRUE(params$verbose %||% FALSE))
}

.trajectory_monocle3_relation_chunks <- function(relation,
    maximum_chunk_pairs = 8192L) {
  row <- seq_along(relation$key)
  chunks <- split(row, ceiling(row / maximum_chunk_pairs))
  lapply(chunks, function(index) unname(Map(function(key, cell_id)
    list(key = key, cell_id = cell_id),
    as.character(relation$key[index]),
    as.character(relation$cell_id[index]))))
}

.run_monocle3_fit <- function(readers, writers, params, context) {
  fit <- .trajectory_monocle3_fit_input(readers$fit_space(),
    context$effective_cells, params)
  writers$declare_derived_fitted_axis(fit$cells)
  cds <- .trajectory_monocle3_fit_cds(fit, params)
  .trajectory_assert(identical(as.character(colnames(cds)), fit$cells),
    "monocle3_fit changed the declared fitted-cell axis.")
  values <- monocle3::pseudotime(cds)
  .trajectory_assert(length(values) == ncol(cds),
    "monocle3_fit: pseudotime length mismatch.")
  values <- values[colnames(cds)]
  pseudotime <- scworkbench::wb_prim_cell_vector(values,
    value_type = "continuous")
  lineage <- .trajectory_monocle3_lineage_paths(cds)
  relation <- .trajectory_monocle3_branch_membership_prim(cds, lineage)
  writers$declare_relation_family_rows("branch_membership",
    ordered_keys = unique(as.character(relation$key)),
    row_chunks = .trajectory_monocle3_relation_chunks(relation),
    numeric_fields = c(weight = TRUE, value = FALSE,
      distance = FALSE))
  writers$write_terminal_primitive("pseudotime", pseudotime)
  writers$write_terminal_primitive("branch_membership", relation)
  fit_space_id <- as.character(context$input_ids$fit_space)[1L]
  tree <- writers$open_model_tree("model")
  writeBin(serialize(cds, NULL, version = 3L),
    tree$path(tree$primary_path))
  reference <- scworkbench::wb_prim_model_ref(fit_space_id,
    fit_space_id, engine = "monocle3")
  reference$model_type <- "monocle3_cds"
  reference$fit_space_id <- fit_space_id
  reference$cells <- fit$cells
  reference$root_pr_nodes <- as.character(lineage$root_nodes)
  reference$lineage_keys <- names(lineage$lineage_paths)
  reference$lineage_paths <- lineage$lineage_paths
  tree$finish(reference)
  invisible(list(status = "ok"))
}

.trajectory_monocle3_read_cds <- function(verified) {
  .trajectory_assert(identical(verified$primary_path, "model.rds"),
    "Authenticated Monocle3 model primary entry changed.")
  cds <- unserialize(verified$read_bytes(verified$primary_path))
  .trajectory_assert(inherits(cds, "cell_data_set"),
    "Authenticated Monocle3 model has an unsupported R object type.")
  fitted <- verified$fitted_source_space
  .trajectory_assert(is.matrix(fitted) && nrow(fitted) > 0L &&
      !is.null(rownames(fitted)) && !anyDuplicated(rownames(fitted)) &&
      all(is.finite(fitted)),
    "Authenticated Monocle3 fitted source matrix is malformed.")
  .trajectory_assert(identical(as.character(colnames(cds)),
      as.character(rownames(fitted))),
    "Authenticated Monocle3 fitted axis differs from the contained CDS.")
  cds
}

.run_monocle3_project <- function(readers, writers, params, context) {
  verified <- readers$fitted_model()
  cds <- .trajectory_monocle3_read_cds(verified)
  fitted <- as.character(colnames(cds))
  display <- .trajectory_space_matrix(readers$plot_space())
  index <- match(fitted, rownames(display))
  .trajectory_assert(!anyNA(index) && ncol(display) == 2L,
    "Display space is not two-dimensional or lacks fitted cells.")
  display <- display[index, , drop = FALSE]
  rownames(display) <- fitted
  lineage <- .trajectory_monocle3_lineage_paths(cds)
  curves <- .trajectory_monocle3_lineage_curves_2d_prim(cds, lineage,
    display, approx_points = as.integer(params$approx_points %||% 200L))
  writers$write_path_2d_family("curves", curves)
  invisible(list(status = "ok"))
}
