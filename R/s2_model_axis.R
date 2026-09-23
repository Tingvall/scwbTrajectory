# A self-linked G1 reader returns a verified private copy of the model tree.
# Only its runtime_private_path capability may name the RDS consumed here.
.trajectory_s2_model_axis <- function(reader) {
  .trajectory_assert(is.function(reader),
    "S2 model path requires an authenticated G1 payload reader.")
  payload <- reader()
  fitted <- payload$fitted_source_space %||% NULL
  .trajectory_assert(is.matrix(fitted) && is.numeric(fitted) &&
      nrow(fitted) > 0L && ncol(fitted) > 0L &&
      !is.null(rownames(fitted)) && !anyNA(rownames(fitted)) &&
      all(nzchar(rownames(fitted))) && !anyDuplicated(rownames(fitted)) &&
      !is.null(colnames(fitted)) && !anyNA(colnames(fitted)) &&
      all(nzchar(colnames(fitted))) && !anyDuplicated(colnames(fitted)) &&
      all(is.finite(fitted)),
    "S2 G1 reader returned no valid authenticated fitted source space.")
  primary <- as.character(payload$primary_path %||% "")
  .trajectory_assert(length(primary) == 1L && !is.na(primary) &&
      nzchar(primary) && is.function(payload$runtime_private_path),
    "S2 model payload reader returned no verified primary entry.")
  path <- payload$runtime_private_path(primary)
  .trajectory_assert(is.character(path) && length(path) == 1L &&
      !is.na(path) && nzchar(path),
    "S2 model payload reader returned no private primary path.")
  model <- readRDS(path)
  if (methods::is(model, "SlingshotDataSet")) {
    .trajectory_assert(requireNamespace("slingshot", quietly = TRUE),
      "S2 Slingshot model extraction requires slingshot.")
    pt <- slingshot::slingPseudotime(model)
    .trajectory_assert(is.matrix(pt) && nrow(pt) > 0L && ncol(pt) > 0L &&
        !is.null(rownames(pt)) && !anyNA(rownames(pt)) &&
        all(nzchar(rownames(pt))) && !anyDuplicated(rownames(pt)),
      "S2 G1 Slingshot model has no valid fitted pseudotime axis.")
    .trajectory_assert(identical(rownames(pt), rownames(fitted)),
      "S2 G1 Slingshot pseudotime cells differ from the authenticated fitted axis.")
    axis <- do.call(rbind, lapply(seq_len(ncol(pt)), function(j)
      data.frame(lineage_id = paste0("L", j),
        cell_id = rownames(pt), pseudotime = as.numeric(pt[, j]),
        stringsAsFactors = FALSE)))
    rownames(axis) <- NULL
    return(list(engine = "slingshot", axis = axis))
  }
  if (methods::is(model, "cell_data_set"))
    return(.trajectory_s2_monocle3_axis(model, rownames(fitted)))
  stop("S2 G1 model type needs an authenticated lineage projection.",
    call. = FALSE)
}

.trajectory_s2_model_lineages <- function(source, lineage = NULL,
                                          cells = NULL,
                                          membership = NULL) {
  .trajectory_assert(is.list(source) &&
      as.character(source$engine %||% "")[1L] %in%
        c("slingshot", "monocle3") &&
      is.data.frame(source$axis) &&
      all(c("lineage_id", "cell_id", "pseudotime") %in% names(source$axis)),
    "S2 model projection is missing its fitted lineage axis.")
  axis <- source$axis[, c("lineage_id", "cell_id", "pseudotime"), drop = FALSE]
  axis$lineage_id <- as.character(axis$lineage_id)
  axis$cell_id <- as.character(axis$cell_id)
  axis$pseudotime <- suppressWarnings(as.numeric(axis$pseudotime))
  .trajectory_assert(nrow(axis) > 0L &&
      !anyNA(axis$lineage_id) && all(nzchar(axis$lineage_id)) &&
      !anyNA(axis$cell_id) && all(nzchar(axis$cell_id)) &&
      !anyDuplicated(axis[c("lineage_id", "cell_id")]),
    "S2 model projection has an invalid lineage/cell axis.")
  if (!is.null(membership)) {
    .trajectory_kind(membership, "cell_vector")
    values <- as.character(membership$value)
    selected <- as.character(membership$cell_id)[!is.na(values) &
      nzchar(values)]
    cells <- if (is.null(cells)) selected else
      intersect(as.character(cells), selected)
  }
  if (!is.null(cells))
    axis <- axis[axis$cell_id %in% as.character(cells), , drop = FALSE]
  requested <- NULL
  if (!is.null(lineage)) {
    requested <- unique(as.character(unlist(lineage, use.names = FALSE)))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    .trajectory_assert(length(requested) > 0L &&
        !length(setdiff(requested, axis$lineage_id)),
      "S2 model lineage selector is empty or unknown.")
    axis <- axis[axis$lineage_id %in% requested, , drop = FALSE]
  }
  .trajectory_assert(nrow(axis) > 0L,
    "S2 model projection has no rows after selection.")
  split(axis, axis$lineage_id)
}

.trajectory_s2_model_source <- function(source, lineage, cells, membership) {
  lineages <- .trajectory_s2_model_lineages(source, lineage, cells,
    membership)
  membership_cells <- if (is.null(membership)) NULL else {
    values <- as.character(membership$value)
    unique(as.character(membership$cell_id)[!is.na(values) &
      nzchar(values)])
  }
  key <- as.character(lineage %||% names(lineages))
  list(lineages = lineages, meta = list(
    pseudotime_source_kind = "model_ref", pseudotime_artifact_id = "",
    pseudotime_key = paste(unique(key[!is.na(key) & nzchar(key)]),
      collapse = ","), membership_artifact_id = "",
    n_pseudotime_cells = as.integer(length(unique(unlist(lapply(
      lineages, function(x) as.character(x$cell_id)), use.names = FALSE)))),
    n_membership_cells = if (is.null(membership_cells)) NA_integer_ else
      as.integer(length(membership_cells))))
}
