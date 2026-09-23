# Project a verified CDS onto the old-core root-to-leaf lineage relation.
# G1 supplies the fitted axis independently of this serialized model.
.trajectory_s2_monocle3_axis <- function(cds, fitted_cells) {
  .trajectory_assert(requireNamespace("monocle3", quietly = TRUE) &&
      requireNamespace("igraph", quietly = TRUE),
    "S2 Monocle3 projection requires monocle3 and igraph.")
  pt <- monocle3::pseudotime(cds)
  cells <- names(pt)
  .trajectory_assert(is.numeric(pt) && length(pt) > 0L &&
      !is.null(cells) && !anyNA(cells) && all(nzchar(cells)) &&
      !anyDuplicated(cells),
    "S2 G1 Monocle3 model has no valid fitted pseudotime axis.")
  .trajectory_assert(identical(as.character(cells), as.character(fitted_cells)),
    "S2 G1 Monocle3 pseudotime cells differ from the authenticated fitted axis.")
  pg <- monocle3::principal_graph(cds)[["UMAP"]] %||% NULL
  aux <- monocle3::principal_graph_aux(cds)[["UMAP"]] %||% NULL
  .trajectory_assert(inherits(pg, "igraph") && is.list(aux),
    "S2 G1 Monocle3 principal graph or auxiliary graph is missing.")
  vertices <- as.character(igraph::V(pg)$name %||% character())
  roots <- unique(as.character(aux$root_pr_nodes %||% character()))
  .trajectory_assert(length(vertices) > 0L && !anyNA(vertices) &&
      all(nzchar(vertices)) && !anyDuplicated(vertices) &&
      length(roots) > 0L && !anyNA(roots) && all(nzchar(roots)) &&
      all(roots %in% vertices) && isTRUE(igraph::is_connected(pg)) &&
      isTRUE(igraph::is_tree(pg, mode = "all")),
    "S2 G1 Monocle3 principal graph must be a named connected tree with valid roots.")
  degrees <- igraph::degree(pg)
  leaves_all <- as.character(names(degrees)[degrees == 1L])
  leaves <- setdiff(leaves_all, roots)
  root_leaves <- !length(leaves)
  if (root_leaves) leaves <- leaves_all
  .trajectory_assert(length(leaves) > 0L,
    "S2 G1 Monocle3 principal graph has no terminal leaves.")
  paths <- vector("list", length(leaves))
  keys <- character(length(leaves))
  for (i in seq_along(leaves)) {
    distance <- as.numeric(igraph::distances(pg, v = roots,
      to = leaves[[i]], mode = "all"))
    .trajectory_assert(length(distance) == length(roots) &&
        any(is.finite(distance)),
      "S2 G1 Monocle3 root-to-leaf path is unreachable.")
    root <- sort(roots[which(is.finite(distance) &
      distance == min(distance[is.finite(distance)]))])[[1L]]
    walk <- igraph::shortest_paths(pg, from = root, to = leaves[[i]],
      output = "vpath")$vpath[[1L]]
    paths[[i]] <- as.character(igraph::V(pg)[walk]$name)
    .trajectory_assert(length(paths[[i]]) >= if (root_leaves) 1L else 2L,
      "S2 G1 Monocle3 root-to-leaf path is invalid.")
    keys[[i]] <- paste0("root:", root, "::lineage:", leaves[[i]])
  }
  closest <- aux$pr_graph_cell_proj_closest_vertex %||% NULL
  .trajectory_assert(!is.null(closest),
    "S2 G1 Monocle3 closest-vertex membership is missing.")
  if (is.matrix(closest) || is.array(closest)) {
    .trajectory_assert(ncol(closest) >= 1L,
      "S2 G1 Monocle3 closest-vertex membership is malformed.")
    membership_cells <- rownames(closest)
    values <- closest[, 1L, drop = TRUE]
  } else {
    membership_cells <- names(closest)
    values <- closest
  }
  .trajectory_assert(!is.null(membership_cells) &&
      length(membership_cells) == length(values) &&
      !anyNA(membership_cells) && all(nzchar(membership_cells)) &&
      !anyDuplicated(membership_cells) &&
      setequal(membership_cells, fitted_cells),
    "S2 G1 Monocle3 closest-vertex cells differ from the authenticated fitted axis.")
  if (is.numeric(values) || is.integer(values)) {
    indices <- suppressWarnings(as.integer(values))
    .trajectory_assert(!anyNA(indices) && all(indices == values) &&
        all(indices >= 1L & indices <= length(vertices)),
      "S2 G1 Monocle3 closest-vertex membership is invalid.")
    closest_names <- vertices[indices]
  } else closest_names <- as.character(values)
  .trajectory_assert(!anyNA(closest_names) &&
      all(closest_names %in% vertices),
    "S2 G1 Monocle3 closest-vertex membership is invalid.")
  rows <- lapply(seq_along(paths), function(i) {
    hit <- closest_names %in% paths[[i]]
    if (!any(hit)) return(NULL)
    data.frame(lineage_id = keys[[i]],
      cell_id = as.character(membership_cells[hit]),
      pseudotime = as.numeric(pt[membership_cells[hit]]),
      stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  .trajectory_assert(length(rows) > 0L,
    "S2 G1 Monocle3 graph produced no branch membership.")
  axis <- do.call(rbind, rows)
  rownames(axis) <- NULL
  list(engine = "monocle3", axis = axis)
}
