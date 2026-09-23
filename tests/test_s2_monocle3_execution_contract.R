source(file.path(test_dir, "test_s2_model_execution_contract.R"),
  local = .GlobalEnv)

.trajectory_s2_branching_cds <- function(cells) {
  assay <- Matrix::Matrix(matrix(1, nrow = 1L, ncol = length(cells),
    dimnames = list("dummy", cells)), sparse = TRUE)
  cds <- monocle3::new_cell_data_set(assay,
    cell_metadata = data.frame(row.names = cells),
    gene_metadata = data.frame(gene_short_name = "dummy",
      row.names = "dummy"))
  graph <- igraph::make_graph(c("A", "B", "B", "C", "B", "D"),
    directed = FALSE)
  closest <- matrix(rep(c("A", "B", "C", "D"),
    times = c(8L, 8L, 4L, 4L)), ncol = 1L,
    dimnames = list(cells, "vertex"))
  methods::slot(cds, "principal_graph") <- S4Vectors::SimpleList(UMAP = graph)
  methods::slot(cds, "principal_graph_aux") <- S4Vectors::SimpleList(list(UMAP = list(
    root_pr_nodes = "A", pr_graph_cell_proj_closest_vertex = closest,
    pseudotime = stats::setNames(seq(0, 1, length.out = length(cells)),
      cells))))
  stopifnot(methods::validObject(cds))
  cds
}

.trajectory_s2_monocle3_core <- function(x, operation, model, params,
                                         cells = x$cells) {
  args <- switch(operation,
    de = list(bundle = x$bundle, matrix_id = "RNA", model = model),
    smooth = list(bundle = x$bundle, matrix_id = "RNA", model = model),
    signature_smooth = list(bundle = x$bundle, matrix_id = "RNA",
      model = model),
    vector_smooth = list(bundle = x$bundle, signal = x$signal,
      model = model),
    geneset_dynamics = list(bundle = x$bundle, matrix_id = "RNA",
      model = model),
    signal_dynamics = list(bundle = x$bundle, signal = x$signal,
      model = model))
  fun <- get(paste0("wb_op_trajectory_", operation),
    envir = asNamespace("scworkbench"))
  params <- params[names(params) %in% names(formals(fun))]
  do.call(fun, c(args, if (!identical(cells, x$cells))
    list(cells = cells) else list(), params))
}

test_trajectory_s2_monocle3_model_oracle_contract <- function() {
  if (!requireNamespace("monocle3", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE))
    return(structure(TRUE, skip = "monocle3 or igraph unavailable"))
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  cds <- .trajectory_s2_branching_cds(x$cells)
  model <- x$model
  model$engine <- "monocle3"
  model$pseudotime_table <- NULL
  model$model_obj <- cds
  reference <- jsonlite::read_json(file.path(test_dir, "fixtures",
    "trajectory_s2_monocle3_reference.json"), simplifyVector = FALSE)
  stopifnot(identical(reference$core_commit,
    "ba24635611756240f21f6817044883b832f12701"),
    identical(reference$comparison_policy, "exact_canonical_value"))
  check_reference <- function(case, expected, actual) {
    trajectory_s2_compare_outputs(x$bundle, case, expected, reference)
    trajectory_s2_compare_outputs(x$bundle, case,
      list(outputs = actual$outputs), reference)
    stopifnot(identical(jsonlite::toJSON(expected$params,
      auto_unbox = TRUE, null = "null", digits = NA),
      jsonlite::toJSON(reference$cases[[case]]$params,
        auto_unbox = TRUE, null = "null", digits = NA)))
  }
  source <- scwbTrajectory:::.trajectory_s2_monocle3_axis(cds, x$cells)
  stopifnot(identical(source$engine, "monocle3"),
    identical(unique(source$axis$lineage_id),
      c("root:A::lineage:C", "root:A::lineage:D")))
  core_axis <- getFromNamespace(".wb_trajectory_de_pt_from_monocle3",
    "scworkbench")(model)
  stopifnot(identical(source$axis, core_axis))
  for (operation in c("de", "smooth", "signature_smooth",
      "geneset_dynamics", "vector_smooth", "signal_dynamics")) {
    params <- .trajectory_s2_model_params(operation)
    params$lineage <- "root:A::lineage:C"
    expected <- .trajectory_s2_monocle3_core(x, operation, model,
      params)
    actual <- .trajectory_s2_run_model_direct(x, operation, params,
      model = cds)
    stopifnot(!inherits(actual$value, "error"),
      identical(names(actual$outputs), names(expected$outputs)))
    check_reference(paste(operation, "C", "full", sep = "__"),
      expected, actual)
    for (name in names(expected$outputs)) {
      stopifnot(identical(actual$outputs[[name]]$kind,
        expected$outputs[[name]]$kind),
        identical(actual$outputs[[name]]$table_type,
          expected$outputs[[name]]$table_type),
        identical(actual$outputs[[name]]$table,
        expected$outputs[[name]]$table),
        identical(actual$outputs[[name]]$meta,
          expected$outputs[[name]]$meta))
      want <- trajectory_identity(x$bundle, expected$outputs[[name]])
      got <- trajectory_identity(x$bundle, actual$outputs[[name]])
      stopifnot(identical(got$content_sha256, want$content_sha256),
        identical(got$structural_fingerprint,
          want$structural_fingerprint))
    }
    alternate <- params
    alternate$lineage <- "root:A::lineage:D"
    expected_alt <- .trajectory_s2_monocle3_core(x, operation, model,
      alternate)
    actual_alt <- .trajectory_s2_run_model_direct(x, operation,
      alternate, model = cds)
    stopifnot(!inherits(actual_alt$value, "error"))
    check_reference(paste(operation, "D", "full", sep = "__"),
      expected_alt, actual_alt)
    for (name in names(expected_alt$outputs)) stopifnot(identical(
      trajectory_identity(x$bundle, actual_alt$outputs[[name]])$
        content_sha256,
      trajectory_identity(x$bundle, expected_alt$outputs[[name]])$
        content_sha256))
    session_cells <- x$cells[-(13:16)]
    expected_subset <- .trajectory_s2_monocle3_core(x, operation,
      model, params, cells = session_cells)
    actual_subset <- .trajectory_s2_run_model_direct(x, operation,
      params, model = cds, effective_cells = session_cells)
    stopifnot(!inherits(actual_subset$value, "error"),
      identical(names(actual_subset$outputs),
        names(expected_subset$outputs)))
    check_reference(paste(operation, "C", "subset", sep = "__"),
      expected_subset, actual_subset)
    for (name in names(expected_subset$outputs)) stopifnot(identical(
      trajectory_identity(x$bundle, actual_subset$outputs[[name]])$
        content_sha256,
      trajectory_identity(x$bundle, expected_subset$outputs[[name]])$
        content_sha256))
  }
  TRUE
}

test_trajectory_s2_monocle3_model_rejection_contract <- function() {
  if (!requireNamespace("monocle3", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE))
    return(structure(TRUE, skip = "monocle3 or igraph unavailable"))
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  cds <- .trajectory_s2_branching_cds(x$cells)
  reordered <- .trajectory_s2_run_model_direct(x, "de", model = cds,
    fitted_cells = rev(x$cells))
  foreign <- .trajectory_s2_run_model_direct(x, "de", model = cds,
    fitted_cells = c(x$cells[-1L], "foreign"))
  broken_graph <- cds
  aux <- methods::slot(broken_graph, "principal_graph")
  aux$UMAP <- igraph::make_ring(4L)
  igraph::V(aux$UMAP)$name <- c("A", "B", "C", "D")
  methods::slot(broken_graph, "principal_graph") <- aux
  invalid_graph <- .trajectory_s2_run_model_direct(x, "de",
    model = broken_graph)
  broken_membership <- cds
  aux <- methods::slot(broken_membership, "principal_graph_aux")
  aux$UMAP$pr_graph_cell_proj_closest_vertex[1L, 1L] <- "foreign"
  methods::slot(broken_membership, "principal_graph_aux") <- aux
  invalid_membership <- .trajectory_s2_run_model_direct(x, "de",
    model = broken_membership)
  for (result in list(reordered, foreign, invalid_graph,
      invalid_membership)) stopifnot(inherits(result$value, "error"),
    length(result$outputs) == 0L)
  TRUE
}
