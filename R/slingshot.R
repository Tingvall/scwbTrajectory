.trajectory_slingshot_packages <- function() {
  for (package in c("SingleCellExperiment", "S4Vectors", "slingshot"))
    .trajectory_assert(requireNamespace(package, quietly = TRUE),
      paste0("Slingshot operation requires package '", package, "'."))
}

.trajectory_slingshot_fit_input <- function(space, clusters, cells, params) {
  .trajectory_kind(clusters, "cell_vector")
  coordinates <- .trajectory_space_matrix(space)
  cells <- as.character(cells)
  .trajectory_assert(length(cells) > 0L && !anyDuplicated(cells) &&
      all(cells %in% rownames(coordinates)),
    "Fitted cells must be a unique ordered subset of the authenticated space.")
  labels <- stats::setNames(as.character(clusters$value),
    as.character(clusters$cell_id))
  .trajectory_assert(!anyNA(labels[cells]),
    "Clusters are missing for selected fit cells.")
  exclude_cells <- as.character(params$exclude_cells %||% character())
  include_clusters <- as.character(params$include_clusters %||% character())
  exclude_clusters <- as.character(params$exclude_clusters %||% character())
  if (length(exclude_cells)) cells <- cells[!cells %in% exclude_cells]
  if (length(include_clusters)) cells <- cells[labels[cells] %in%
    include_clusters]
  if (length(exclude_clusters)) cells <- cells[!labels[cells] %in%
    exclude_clusters]
  .trajectory_assert(length(cells) >= 10L,
    "Too few cells after Slingshot filtering (need >= 10).")
  selected <- coordinates[match(cells, rownames(coordinates)), ,
    drop = FALSE]
  cluster <- factor(unname(labels[cells]))
  present <- sort(unique(as.character(cluster)))
  start <- as.character(params[["start.clus"]] %||% character())
  end <- as.character(params[["end.clus"]] %||% character())
  .trajectory_assert(length(start) <= 1L && all(start %in% present),
    "start.clus is absent after filtering.")
  .trajectory_assert(all(end %in% present),
    "end.clus is absent after filtering.")
  list(values = selected, cells = cells, clusters = cluster,
    source_count = length(as.character(space$cell_id %||%
      rownames(coordinates))))
}

.trajectory_slingshot_fit <- function(prepared, params) {
  .trajectory_slingshot_packages()
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(dummy =
    matrix(0, nrow = 1L, ncol = length(prepared$cells))))
  colnames(sce) <- prepared$cells
  SingleCellExperiment::reducedDims(sce) <- S4Vectors::SimpleList(
    space = prepared$values)
  slingshot::slingshot(sce, clusterLabels = prepared$clusters,
    reducedDim = "space", start.clus = params[["start.clus"]] %||% NULL,
    end.clus = params[["end.clus"]] %||% NULL,
    stretch = as.numeric(params$stretch %||% 2),
    reweight = isTRUE(params$reweight %||% TRUE),
    shrink = as.numeric(params$shrink %||% 1),
    maxit = as.integer(params$maxit %||% 10))
}

.trajectory_slingshot_family <- function(model, accessor, fit_id,
    fit_space_id, family_kind, semantics) {
  values <- as.matrix(accessor(model))
  cells <- rownames(values)
  .trajectory_assert(nrow(values) > 0L && ncol(values) > 0L &&
      !is.null(cells) && !anyDuplicated(cells),
    "Slingshot family has no unique fitted cell axis.")
  keys <- paste0("L", seq_len(ncol(values)))
  retained <- lapply(seq_len(ncol(values)), function(i) {
    index <- which(!is.na(values[, i]))
    list(key = rep(keys[[i]], length(index)), cell_id = cells[index],
      value = as.numeric(values[index, i]))
  })
  key <- unlist(lapply(retained, `[[`, "key"), use.names = FALSE)
  cell_id <- unlist(lapply(retained, `[[`, "cell_id"), use.names = FALSE)
  value <- unlist(lapply(retained, `[[`, "value"), use.names = FALSE)
  .trajectory_assert(length(value) > 0L && all(is.finite(value)),
    "Slingshot family contains no finite values.")
  if (identical(family_kind, "trajectory_membership_weight")) {
    .trajectory_assert(all(value >= -1e-8 & value <= 1 + 1e-8),
      "Slingshot lineage weights must lie in [0, 1].")
    value <- pmax(0, pmin(1, value))
  }
  scworkbench::wb_prim_cell_vector_family(key, cell_id, value,
    value_type = "continuous", meta = list(trajectory_fit_id = fit_id,
      engine = "slingshot", fit_space_id = fit_space_id,
      n_lineages = as.integer(ncol(values)), family_kind = family_kind,
      value_semantics = semantics))
}

.trajectory_slingshot_sparse_layout <- function(writers, output, family) {
  keys <- unique(as.character(family$key))
  cells <- lapply(keys, function(key) as.character(family$cell_id[
    family$key == key]))
  writers$declare_sparse_cell_family_layout(output, keys, cells)
}

.trajectory_slingshot_curves <- function(model, fit_id = NULL,
    fit_space_id = NULL, semantics = NULL, maximum_points = NULL) {
  curves <- slingshot::slingCurves(model)
  .trajectory_assert(length(curves) > 0L,
    "Slingshot model has no fitted lineage curves.")
  rows <- lapply(seq_along(curves), function(i) {
    coords <- as.matrix(curves[[i]]$s)
    .trajectory_assert(nrow(coords) > 0L && ncol(coords) >= 2L &&
        all(is.finite(coords[, 1:2, drop = FALSE])),
      "Slingshot curve is not finite and two-dimensional.")
    if (!is.null(maximum_points) && is.finite(maximum_points) &&
        maximum_points > 0L && nrow(coords) > maximum_points)
      coords <- coords[unique(round(seq.int(1L, nrow(coords),
        length.out = maximum_points))), , drop = FALSE]
    data.frame(key = paste0("L", i), point = seq_len(nrow(coords)),
      x = as.numeric(coords[, 1L]), y = as.numeric(coords[, 2L]))
  })
  rows <- do.call(rbind, rows)
  meta <- if (is.null(fit_id)) list() else list(
    trajectory_fit_id = fit_id, engine = "slingshot",
    fit_space_id = fit_space_id,
    n_lineages = as.integer(length(curves)),
    family_kind = "trajectory_curve_2d", value_semantics = semantics)
  scworkbench::wb_prim_path_2d_family(rows$key, rows$point, rows$x,
    rows$y, meta = meta)
}

.trajectory_slingshot_model_bytes <- function(model) serialize(model, NULL,
  version = 3L)

.trajectory_slingshot_read_model <- function(reader) {
  primary <- reader$primary_path
  .trajectory_assert(identical(primary, "model.rds"),
    "Authenticated Slingshot model primary entry changed.")
  model <- unserialize(reader$read_bytes(primary))
  .trajectory_assert(inherits(model, "SingleCellExperiment") ||
      inherits(model, "SlingshotDataSet"),
    "Authenticated Slingshot model has an unsupported R object type.")
  model
}

.trajectory_slingshot_fit_id <- function(bytes) paste0("trajfit_slingshot_",
  substr(digest::digest(bytes, algo = "sha256", serialize = FALSE), 1L, 20L))

.run_slingshot <- function(readers, writers, params, context) {
  fit <- .trajectory_slingshot_fit_input(readers$space_2d(),
    readers$clusters(), context$effective_cells, params)
  model <- .trajectory_slingshot_fit(fit, params)
  pseudotime <- .trajectory_slingshot_family(model,
    slingshot::slingPseudotime, "", as.character(context$input_ids$space_2d),
    "trajectory_pseudotime", "lineage_specific_pseudotime")
  curves <- .trajectory_slingshot_curves(model)
  .trajectory_slingshot_sparse_layout(writers, "pseudotime", pseudotime)
  writers$write_terminal_primitive("pseudotime", pseudotime)
  writers$write_path_2d_family("curves", curves)
  invisible(list(status = "ok"))
}

.run_slingshot_fit <- function(readers, writers, params, context) {
  space <- readers$fit_space()
  fit <- .trajectory_slingshot_fit_input(space, readers$clusters(),
    context$effective_cells, params)
  writers$declare_derived_fitted_axis(fit$cells)
  model <- .trajectory_slingshot_fit(fit, params)
  bytes <- .trajectory_slingshot_model_bytes(model)
  fit_id <- .trajectory_slingshot_fit_id(bytes)
  fit_space_id <- as.character(context$input_ids$fit_space)[1L]
  pseudotime <- .trajectory_slingshot_family(model,
    slingshot::slingPseudotime, fit_id, fit_space_id,
    "trajectory_pseudotime", "lineage_specific_pseudotime")
  weights <- .trajectory_slingshot_family(model,
    slingshot::slingCurveWeights, fit_id, fit_space_id,
    "trajectory_membership_weight", "lineage_membership_weight")
  writers$write_terminal_primitive("pseudotime", pseudotime)
  writers$write_terminal_primitive("lineage_weights", weights)
  tree <- writers$open_model_tree("model")
  writeBin(bytes, tree$path(tree$primary_path))
  reference <- scworkbench::wb_prim_model_ref(fit_space_id, fit_space_id,
    engine = "slingshot")
  reference$model_type <- "slingshot"
  reference$trajectory_fit_id <- fit_id
  reference$cells <- fit$cells
  reference$n_lineages <- length(unique(pseudotime$key))
  reference$fit_space_id <- fit_space_id
  tree$finish(reference)
  invisible(list(status = "ok"))
}

.run_slingshot_lineage_weights <- function(readers, writers, params,
    context) {
  verified <- readers$fitted_model()
  model <- .trajectory_slingshot_read_model(verified)
  bytes <- verified$read_bytes(verified$primary_path)
  fitted_space <- verified$fitted_source_space
  fit_space_id <- if (is.list(fitted_space))
    as.character(fitted_space$space_id %||% "")[1L] else ""
  weights <- .trajectory_slingshot_family(model,
    slingshot::slingCurveWeights, .trajectory_slingshot_fit_id(bytes),
    fit_space_id, "trajectory_membership_weight",
    "lineage_membership_weight")
  .trajectory_slingshot_sparse_layout(writers, "lineage_weights", weights)
  writers$write_terminal_primitive("lineage_weights", weights)
  invisible(list(status = "ok"))
}

.run_slingshot_project <- function(readers, writers, params, context) {
  verified <- readers$fitted_model()
  model <- .trajectory_slingshot_read_model(verified)
  fitted_cells <- rownames(slingshot::slingPseudotime(model))
  display <- .trajectory_space_matrix(readers$plot_space())
  index <- match(fitted_cells, rownames(display))
  .trajectory_assert(!anyNA(index) && ncol(display) == 2L,
    "Display space is not two-dimensional or lacks fitted cells.")
  display <- display[index, , drop = FALSE]
  projected <- slingshot::embedCurves(model, display,
    stretch = params$stretch %||% NULL,
    shrink = params$shrink %||% NULL,
    approx_points = as.integer(params$approx_points %||% 200L))
  fitted_space <- verified$fitted_source_space
  fit_space_id <- if (is.list(fitted_space))
    as.character(fitted_space$space_id %||% "")[1L] else ""
  curves <- .trajectory_slingshot_curves(projected,
    .trajectory_slingshot_fit_id(verified$read_bytes(verified$primary_path)),
    fit_space_id,
    "projected_lineage_curve",
    as.integer(params$approx_points %||% 200L))
  writers$write_path_2d_family("curves", curves)
  invisible(list(status = "ok"))
}
