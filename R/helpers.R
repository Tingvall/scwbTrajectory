`%||%` <- function(x, y) if (is.null(x)) y else x

.trajectory_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.trajectory_kind <- function(value, kind) {
  .trajectory_assert(inherits(value, "wb_primitive") &&
      identical(as.character(value$kind)[1L], kind),
    paste0("Expected a ", kind, " primitive."))
  invisible(TRUE)
}

.trajectory_space_matrix <- function(space) {
  kind <- as.character(space$kind %||% "")[1L]
  .trajectory_assert(inherits(space, "wb_primitive") &&
      kind %in% c("space", "space_2d"),
    "Expected a space or space_2d primitive.")
  values <- as.matrix(space$values)
  .trajectory_assert(nrow(values) > 0L && ncol(values) >= 2L &&
      !is.null(rownames(values)) && !anyDuplicated(rownames(values)) &&
      all(is.finite(values)),
    "Space coordinates must be finite with a unique non-empty cell axis.")
  values
}

.trajectory_path_matrix <- function(path) {
  .trajectory_kind(path, "space_path_family")
  path_ids <- unique(as.character(path$path_id))
  .trajectory_assert(length(path_ids) == 1L,
    "S1 path operations require exactly one path_id.")
  points <- sort(unique(as.integer(path$point)))
  dimensions <- sort(unique(as.integer(path$dim)))
  .trajectory_assert(length(points) >= 2L && length(dimensions) >= 2L,
    "Analytical paths require at least two points and two dimensions.")
  matrix <- matrix(NA_real_, length(points), length(dimensions),
    dimnames = list(as.character(points), paste0("dim", dimensions)))
  row <- match(as.integer(path$point), points)
  column <- match(as.integer(path$dim), dimensions)
  .trajectory_assert(!anyDuplicated(paste(row, column, sep = "\r")),
    "space_path_family contains duplicate point/dimension rows.")
  matrix[cbind(row, column)] <- as.numeric(path$value)
  .trajectory_assert(all(is.finite(matrix)),
    "space_path_family contains incomplete or non-finite coordinates.")
  list(path_id = path_ids[[1L]], points = matrix,
    dimensions = dimensions)
}

.trajectory_anchor_matrix <- function(path) {
  anchors <- path$meta$anchor_centroids %||% NULL
  if (is.null(anchors)) return(NULL)
  anchors <- as.matrix(anchors)
  storage.mode(anchors) <- "double"
  .trajectory_assert(ncol(anchors) == length(unique(path$dim)) &&
      all(is.finite(anchors)), "Path anchor centroids are malformed.")
  anchors
}

.trajectory_feature_sets <- function(readers, params, matrix_features,
                                      required = FALSE, prefix) {
  explicit <- params$features %||% NULL
  if (!is.null(explicit)) {
    explicit <- unique(as.character(unlist(explicit, use.names = FALSE)))
    explicit <- explicit[nzchar(explicit)]
    .trajectory_assert(length(explicit) > 0L,
      paste0(prefix, ": features is empty after cleaning."))
  }
  set_features <- NULL
  if (!is.null(params$feature_set)) {
    requested <- as.character(params$feature_set)[1L]
    .trajectory_assert(nzchar(requested), paste0(prefix,
      ": feature_set must be a non-empty string when provided."))
    snapshot <- readers$feature_sets()
    sets <- snapshot$sets %||% list()
    match_set <- Filter(function(item) identical(
      as.character(item$feature_set %||% "")[1L], requested), sets)
    .trajectory_assert(length(match_set) == 1L, paste0(prefix,
      ": feature_set '", requested,
      "' does not resolve uniquely in the authenticated registry snapshot."))
    set_features <- vapply(match_set[[1L]]$memberships %||% list(),
      function(item) as.character(item$feature_id %||% "")[1L],
      character(1))
    set_features <- unique(set_features[nzchar(set_features)])
    .trajectory_assert(length(set_features) > 0L, paste0(prefix,
      ": feature_set '", requested, "' has no matrix members."))
  } else {
    # Registry reads are single-use. Consume the authenticated snapshot even
    # when the optional selector is absent so the receipt is complete.
    readers$feature_sets()
  }
  selected <- if (is.null(set_features)) explicit else if (is.null(explicit))
    set_features else intersect(set_features, explicit)
  if (is.null(selected)) selected <- matrix_features
  selected <- selected[selected %in% matrix_features]
  .trajectory_assert(length(selected) > 0L, paste0(prefix,
    if (required) " requires features or feature_set." else
      ": no features remain after selection."))
  selected
}

.trajectory_chunked_matrix <- function(reader, all_features, selected,
                                        read_cells, selected_cells,
                                        chunk_size) {
  chunks <- scworkbench::wb_partition_axis(all_features,
    as.integer(chunk_size))
  retained <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    features <- chunks[[i]]
    block <- reader(features = features, cells = read_cells)$values
    keep <- features %in% selected
    retained[[i]] <- block[keep, , drop = FALSE]
  }
  values <- do.call(rbind, retained)
  values <- values[match(selected, rownames(values)), , drop = FALSE]
  .trajectory_assert(!anyNA(match(selected, rownames(values))),
    "Authenticated matrix chunks did not reconstruct the selected feature axis.")
  index <- match(selected_cells, colnames(values))
  .trajectory_assert(!anyNA(index),
    "Selected cells exceed the authenticated matrix cell axis.")
  values[, index, drop = FALSE]
}

.trajectory_content_identity <- function(bundle, value) {
  utils::getFromNamespace(".wb_prov_payload_identity", "scworkbench")(bundle,
    value)
}
