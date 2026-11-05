# A small, dependency-free method corpus. The model_ref carries a synthetic
# pseudotime table for the old-core numerical oracle; it is not a G1 reopen
# fixture. A retained S3/S5 model will replace it at integration.
trajectory_s2_fixture <- function() {
  x <- trajectory_fixture()
  x$model <- wb_prim_model_ref(model_space_id = "latent",
    in_space_id = "latent", engine = "slingshot")
  x$model$pseudotime_table <- data.frame(
    lineage_id = rep(c("L1", "L2"), each = length(x$cells)),
    cell_id = rep(x$cells, 2L),
    pseudotime = c(seq(0, 1, length.out = length(x$cells)),
      seq(1, 0, length.out = length(x$cells))),
    stringsAsFactors = FALSE)
  x$explicit <- wb_prim_cell_vector_family(
    rep(c("L1", "L2"), each = length(x$cells)),
    rep(x$cells, 2L), x$model$pseudotime_table$pseudotime,
    value_type = "continuous")
  x$membership <- wb_prim_cell_vector(
    stats::setNames(rep("selected", length(x$cells)), x$cells),
    value_type = "categorical")
  x$signal <- wb_prim_cell_vector_family(
    rep(c("rising", "falling"), each = length(x$cells)),
    rep(x$cells, 2L),
    c(seq(0, 1, length.out = length(x$cells)),
      seq(1, 0, length.out = length(x$cells))),
    value_type = "continuous")
  x
}

trajectory_s2_result_record <- function(bundle, result) {
  stopifnot(is.list(result), is.list(result$outputs),
    is.list(result$params), length(result$outputs) > 0L)
  list(outputs = lapply(result$outputs, function(value) {
    stopifnot(inherits(value, "wb_primitive"))
    list(kind = as.character(value$kind),
      table_type = as.character(value$table_type %||% ""),
      rows = if (is.data.frame(value$table)) nrow(value$table) else NA_integer_,
      identity = trajectory_identity(bundle, value))
  }), params = result$params)
}

trajectory_s2_compare_outputs <- function(bundle, operation, result, reference) {
  expected <- reference$cases[[operation]]$outputs
  if (is.null(expected)) stop("Unknown S2 reference case: ", operation)
  actual <- result$outputs
  allowed_extra <- if (grepl("geneset_dynamics", operation, fixed = TRUE))
    c("adjusted_scores", "adjustment_diagnostics") else character()
  if (!setequal(intersect(names(actual), names(expected)), names(expected)) ||
      length(setdiff(names(actual), names(expected))) > length(allowed_extra) ||
      !setequal(setdiff(names(actual), names(expected)), allowed_extra))
    stop(operation, ": terminal closure differs from the reference")
  for (name in names(expected)) {
    value <- actual[[name]]
    want <- expected[[name]]
    if (!inherits(value, "wb_primitive") ||
        !identical(as.character(value$kind), want$kind) ||
        !identical(as.character(value$table_type %||% ""), want$table_type) ||
        !identical(nrow(value$table), as.integer(want$rows)))
      stop(operation, "/", name, ": terminal shape differs")
    observed <- trajectory_identity(bundle, value)
    for (field in c("structural_fingerprint", "content_sha256",
        "content_identity_scheme", "verification_level",
        "verification_status")) {
      if (!identical(as.character(observed[[field]]), want$identity[[field]]))
        stop(operation, "/", name, ": ", field, " differs")
    }
  }
  invisible(TRUE)
}

trajectory_s2_reference <- function() {
  directory <- if (exists("test_dir", inherits = TRUE))
    get("test_dir", inherits = TRUE) else "tests"
  jsonlite::read_json(file.path(directory, "fixtures",
    "trajectory_s2_core_reference.json"), simplifyVector = FALSE)
}

trajectory_s2_mock_boundary <- function(x, fault_terminal = NULL) {
  declaration <- list(matrix_input = "matrix",
    no_selector_policy = "all_namespaces", empty_match_policy = "allow")
  snapshot <- getFromNamespace(
    ".wb_operation_pack_v2_feature_set_registry_snapshot",
    "scworkbench")(x$bundle, declaration,
      list(matrix = wb_prim_matrix("RNA")), list())
  state <- new.env(parent = emptyenv())
  state$calls <- list()
  state$writes <- list()
  read_once <- function(name, value) {
    force(name); force(value)
    function() {
      state$calls[[name]] <- as.integer(state$calls[[name]] %||% 0L) + 1L
      if (state$calls[[name]] != 1L)
        stop("Typed reader called more than once: ", name)
      value
    }
  }
  read_matrix <- function(features, cells) {
    stopifnot(identical(as.character(cells), x$cells),
      all(features %in% x$features))
    state$calls$matrix_chunks <- c(state$calls$matrix_chunks,
      as.character(features))
    list(values = x$values[features, cells, drop = FALSE])
  }
  readers <- list(pseudotime = read_once("pseudotime", x$explicit),
    membership = read_once("membership", x$membership),
    signal = read_once("signal", x$signal),
    feature_sets = read_once("feature_sets", snapshot),
    matrix_chunks = read_matrix)
  writers <- list(write_terminal_primitive = function(name, value) {
    if (identical(name, fault_terminal)) stop("injected terminal fault")
    state$writes[[name]] <- value
    invisible(TRUE)
  })
  context <- list(effective_cells = x$cells,
    typed_slice_reads = list(matrix_chunks = list(features = x$features)))
  list(readers = readers, writers = writers, context = context,
    state = state, snapshot = snapshot)
}
