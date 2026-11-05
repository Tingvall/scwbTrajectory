# Explicit S2 modes have their own registered entries because a self-linked
# G1 read requires a present model input for each model-mode operation.
.trajectory_s2_call <- function(method, fixed, params,
                                extra_allowed = character()) {
  .trajectory_assert(is.list(params) && !anyDuplicated(names(params)),
    "S2 runner parameters must be a uniquely named list.")
  allowed <- setdiff(names(formals(method)), names(fixed))
  .trajectory_assert(!length(setdiff(names(params),
    c(allowed, extra_allowed))),
    "S2 runner received an undeclared parameter.")
  supplied <- params[intersect(names(params), allowed)]
  do.call(method, c(fixed, supplied))
}

.trajectory_s2_publish <- function(result, writers, terminal_names) {
  .trajectory_assert(identical(names(result$outputs), terminal_names),
    "S2 method returned an incomplete or reordered terminal closure.")
  for (name in terminal_names)
    writers$write_terminal_primitive(name, result$outputs[[name]])
  invisible(list(status = "ok"))
}

.trajectory_s2_source_artifact_id <- function(readers, name) {
  reader <- readers[[name]]
  if (is.null(reader)) return("")
  .trajectory_assert(is.function(reader),
    "S2 artifact ancestry capability is malformed.")
  ancestry <- reader()
  id <- as.character(ancestry$artifact_id %||% "")[1L]
  .trajectory_assert(length(id) == 1L && !is.na(id) && nzchar(id),
    "S2 artifact ancestry has no authenticated artifact ID.")
  id
}

.run_trajectory_geneset_dynamics_explicit <- function(readers, writers,
                                                      params, context) {
  pseudotime <- readers$pseudotime()
  membership <- if (is.function(readers$membership)) readers$membership() else
    NULL
  registry <- readers$feature_sets()
  all_features <- as.character(context$typed_slice_reads$matrix_chunks$features)
  .trajectory_assert(length(all_features) > 0L,
    "trajectory_geneset_dynamics: authenticated matrix feature axis is empty.")
  .trajectory_assert(identical(as.character(unlist(registry$feature_universe,
    use.names = FALSE)), all_features),
    "trajectory_geneset_dynamics: registry and matrix feature axes differ.")
  namespace <- as.character(params$namespace %||% "")[1L]
  .trajectory_assert(nzchar(namespace),
    "trajectory_geneset_dynamics: namespace is required.")
  min_set_size <- as.integer(params$min_set_size %||% 5L)
  max_set_size <- as.integer(params$max_set_size %||% 5000L)
  .trajectory_assert(length(min_set_size) == 1L &&
    is.finite(min_set_size) && min_set_size >= 1L &&
    length(max_set_size) == 1L && is.finite(max_set_size) &&
    max_set_size >= min_set_size,
    "trajectory_geneset_dynamics: invalid feature-set size bounds.")
  selected <- .trajectory_s2_registry_sets(registry, namespace,
    params$feature_sets %||% NULL, all_features, min_set_size,
    max_set_size)
  features <- unique(unlist(selected$members, use.names = FALSE))
  cells <- as.character(context$effective_cells)
  depth_input <- .trajectory_regions_optional_reader(readers, "depth")
  depth <- if (isTRUE(params$adjust_depth %||% TRUE)) {
    .trajectory_assert(!is.null(depth_input), paste0(
      "trajectory_geneset_dynamics: adjust_depth=TRUE requires an ",
      "authenticated depth input."))
    .trajectory_depth_vector(depth_input, cells,
      label = "trajectory_geneset_dynamics")
  } else NULL
  matrix <- as.matrix(.trajectory_chunked_matrix(readers$matrix_chunks,
    all_features, features, cells, cells,
    params$chunk_size %||% 128L))
  result <- .trajectory_s2_call(.trajectory_s2_geneset_dynamics,
    list(matrix = matrix, registry = registry, effective_cells = cells,
      matrix_id = as.character(registry$matrix_id)[1L],
      depth = depth, pseudotime = pseudotime, membership = membership,
      depth_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "depth_ancestry"),
      pseudotime_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "pseudotime_ancestry"),
      membership_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "membership_ancestry")), params,
    extra_allowed = "chunk_size")
  .trajectory_s2_publish(result, writers, c("signature_table",
    "dynamics_table", "top_sets", "driver_table", "driver_curve_table",
    "adjusted_scores", "adjustment_diagnostics", "fit_summary"))
}

.run_trajectory_signal_dynamics_explicit <- function(readers, writers,
                                                     params, context) {
  signal <- readers$signal()
  pseudotime <- readers$pseudotime()
  membership <- if (is.function(readers$membership)) readers$membership() else
    NULL
  result <- .trajectory_s2_call(.trajectory_s2_signal_dynamics,
    list(effective_cells = as.character(context$effective_cells),
      signal = signal, pseudotime = pseudotime, membership = membership,
      pseudotime_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "pseudotime_ancestry"),
      membership_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "membership_ancestry")),
    params)
  .trajectory_s2_publish(result, writers, c("signal_table",
    "dynamics_table", "top_signals", "fit_summary"))
}
