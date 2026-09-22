# Method-side selection for the two S2 dynamics operations. The caller must
# supply a cell_vector_family obtained through an authenticated typed input;
# this helper has no access to bundle or registry state.
.trajectory_s2_explicit_lineages <- function(pseudotime, pseudotime_key = NULL,
                                             lineage = NULL, cells = NULL,
                                             membership = NULL) {
  .trajectory_kind(pseudotime, "cell_vector_family")
  table <- scworkbench::wb_primitive_to_table(pseudotime,
    include_context = FALSE)
  .trajectory_assert(all(c("key", "cell_id", "value") %in% names(table)),
    "trajectory dynamics: explicit pseudotime table is malformed.")
  table <- data.frame(lineage_id = as.character(table$key),
    cell_id = as.character(table$cell_id),
    pseudotime = suppressWarnings(as.numeric(table$value)),
    stringsAsFactors = FALSE)
  table <- table[!is.na(table$lineage_id) & !is.na(table$cell_id) &
    nzchar(table$lineage_id) & nzchar(table$cell_id) &
    is.finite(table$pseudotime), , drop = FALSE]
  .trajectory_assert(nrow(table) > 0L,
    "trajectory dynamics: explicit pseudotime has no finite rows.")
  .trajectory_assert(!anyDuplicated(table[c("lineage_id", "cell_id")]),
    "trajectory dynamics: explicit pseudotime repeats a key/cell pair.")

  key <- pseudotime_key %||% lineage
  if (is.null(key)) {
    key <- sort(unique(table$lineage_id))
    .trajectory_assert(length(key) == 1L,
      "trajectory dynamics: explicit pseudotime has multiple keys; provide pseudotime_key or lineage.")
  } else {
    key <- unique(as.character(unlist(key, use.names = FALSE)))
    key <- key[nzchar(key)]
    .trajectory_assert(length(key) > 0L,
      "trajectory dynamics: pseudotime_key/lineage is empty.")
    .trajectory_assert(!length(setdiff(key, table$lineage_id)),
      "trajectory dynamics: explicit pseudotime key not found.")
  }
  table <- table[table$lineage_id %in% key, , drop = FALSE]

  if (!is.null(membership)) {
    .trajectory_kind(membership, "cell_vector")
    membership_value <- as.character(membership$value)
    keep <- as.character(membership$cell_id)[!is.na(membership_value) &
      nzchar(membership_value)]
    cells <- if (is.null(cells)) keep else intersect(as.character(cells), keep)
  }
  if (!is.null(cells))
    table <- table[table$cell_id %in% as.character(cells), , drop = FALSE]
  .trajectory_assert(nrow(table) > 0L,
    "trajectory dynamics: no explicit pseudotime rows remain after filters.")
  split(table, table$lineage_id)
}
