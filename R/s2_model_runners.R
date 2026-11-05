# Grouped S2 model consumers use G1, registry, typed matrix, and terminal
# grants declared in the pack manifest.
.trajectory_s2_model_matrix_inputs <- function(readers, params, context,
                                               required, prefix) {
  snapshot <- readers$feature_sets()
  all_features <- as.character(context$typed_slice_reads$matrix_chunks$features)
  .trajectory_assert(length(all_features) > 0L &&
      !anyDuplicated(all_features),
    paste0(prefix, ": authenticated matrix feature axis is invalid."))
  selected <- .trajectory_s2_model_features(snapshot, all_features,
    feature_set = params$feature_set %||% NULL,
    features = params$features %||% NULL,
    required = required, prefix = prefix)
  cells <- as.character(context$effective_cells)
  matrix <- as.matrix(.trajectory_chunked_matrix(readers$matrix_chunks,
    all_features, selected, cells, cells,
    params$chunk_size %||% 128L))
  list(matrix = matrix, cells = cells,
    matrix_id = as.character(snapshot$matrix_id)[1L])
}

.run_trajectory_de_model <- function(readers, writers, params, context) {
  source <- .trajectory_s2_model_axis(readers$model_payload)
  snapshot <- readers$feature_sets()
  all_features <- as.character(context$typed_slice_reads$matrix_chunks$features)
  selected <- .trajectory_s2_model_features(snapshot, all_features,
    feature_set = params$feature_set %||% NULL,
    features = params$features %||% NULL,
    required = FALSE, prefix = "trajectory_de")
  result <- .trajectory_s2_call(.trajectory_s2_de_stream,
    list(reader = readers$matrix_chunks, all_features = all_features,
      selected_features = selected, source = source,
      effective_cells = as.character(context$effective_cells),
      matrix_id = as.character(snapshot$matrix_id)[1L]), params)
  .trajectory_s2_publish(result, writers,
    c("de_table", "fit_summary"))
}

.run_trajectory_smooth_model <- function(readers, writers, params,
                                         context) {
  source <- .trajectory_s2_model_axis(readers$model_payload)
  input <- .trajectory_s2_model_matrix_inputs(readers, params, context,
    TRUE, "trajectory_smooth")
  result <- .trajectory_s2_call(.trajectory_s2_smooth_model,
    list(matrix = input$matrix, source = source,
      effective_cells = input$cells, matrix_id = input$matrix_id),
    params, extra_allowed = c("chunk_size", "feature_set", "features"))
  .trajectory_s2_publish(result, writers,
    c("trend_table", "fit_summary"))
}

.run_trajectory_signature_smooth_model <- function(readers, writers,
                                                   params, context) {
  source <- .trajectory_s2_model_axis(readers$model_payload)
  input <- .trajectory_s2_model_matrix_inputs(readers, params, context,
    TRUE, "trajectory_signature_smooth")
  result <- .trajectory_s2_call(.trajectory_s2_signature_model,
    list(matrix = input$matrix, source = source,
      effective_cells = input$cells, matrix_id = input$matrix_id),
    params, extra_allowed = c("chunk_size", "features"))
  .trajectory_s2_publish(result, writers,
    c("signature_table", "fit_summary"))
}

.run_trajectory_vector_smooth_model <- function(readers, writers,
                                                params, context) {
  source <- .trajectory_s2_model_axis(readers$model_payload)
  signal <- readers$signal()
  result <- .trajectory_s2_call(.trajectory_s2_vector_model,
    list(signal = signal, source = source,
      effective_cells = as.character(context$effective_cells)), params)
  .trajectory_s2_publish(result, writers,
    c("signal_table", "fit_summary"))
}

.run_trajectory_geneset_dynamics_model <- function(readers, writers,
                                                   params, context) {
  .trajectory_assert(is.null(params$pseudotime_key),
    "Model mode does not accept pseudotime_key.")
  source <- .trajectory_s2_model_axis(readers$model_payload)
  membership <- if (is.function(readers$membership)) readers$membership() else
    NULL
  registry <- readers$feature_sets()
  all_features <- as.character(context$typed_slice_reads$matrix_chunks$features)
  .trajectory_assert(length(all_features) > 0L &&
      identical(as.character(unlist(registry$feature_universe,
        use.names = FALSE)), all_features),
    "trajectory_geneset_dynamics: registry and matrix axes differ.")
  namespace <- as.character(params$namespace %||% "")[1L]
  min_size <- as.integer(params$min_set_size %||% 5L)
  max_size <- as.integer(params$max_set_size %||% 5000L)
  .trajectory_assert(nzchar(namespace) && is.finite(min_size) &&
      min_size >= 1L && is.finite(max_size) && max_size >= min_size,
    "trajectory_geneset_dynamics: invalid namespace or set bounds.")
  selected <- .trajectory_s2_registry_sets(registry, namespace,
    params$feature_sets %||% NULL, all_features, min_size, max_size)
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
      depth = depth, model_source = source, membership = membership,
      depth_artifact_id = .trajectory_s2_source_artifact_id(readers,
        "depth_ancestry")), params,
    extra_allowed = "chunk_size")
  .trajectory_s2_publish(result, writers, c("signature_table",
    "dynamics_table", "top_sets", "driver_table", "driver_curve_table",
    "adjusted_scores", "adjustment_diagnostics", "fit_summary"))
}

.run_trajectory_signal_dynamics_model <- function(readers, writers,
                                                  params, context) {
  .trajectory_assert(is.null(params$pseudotime_key),
    "Model mode does not accept pseudotime_key.")
  source <- .trajectory_s2_model_axis(readers$model_payload)
  signal <- readers$signal()
  membership <- if (is.function(readers$membership)) readers$membership() else
    NULL
  result <- .trajectory_s2_call(.trajectory_s2_signal_dynamics,
    list(effective_cells = as.character(context$effective_cells),
      signal = signal, model_source = source, membership = membership),
    params)
  .trajectory_s2_publish(result, writers,
    c("signal_table", "dynamics_table", "top_signals", "fit_summary"))
}
