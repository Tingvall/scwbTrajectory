test_dir <- if (file.exists("helper_s2_fixture.R")) "." else "tests"
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

test_trajectory_s2_monocle3_model_oracle_contract <- function() {
  if (!requireNamespace("monocle3", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE))
    return(structure(TRUE, skip = "monocle3 or igraph unavailable"))
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  cds <- .trajectory_s2_branching_cds(x$cells)
  reference <- jsonlite::read_json(file.path(test_dir, "fixtures",
    "trajectory_s2_monocle3_reference.json"), simplifyVector = FALSE)
  stopifnot(identical(reference$core_commit,
    "ba24635611756240f21f6817044883b832f12701"),
    identical(reference$comparison_policy, "exact_canonical_value"))
  source <- scwbTrajectory:::.trajectory_s2_monocle3_axis(cds, x$cells)
  indices <- list(C = c(seq_len(16L), 17:20),
    D = c(seq_len(16L), 21:24))
  pseudotime <- seq(0, 1, length.out = length(x$cells))
  expected_axis <- do.call(rbind, lapply(names(indices), function(leaf) {
    keep <- indices[[leaf]]
    data.frame(lineage_id = rep(paste0("root:A::lineage:", leaf),
      length(keep)), cell_id = x$cells[keep],
      pseudotime = pseudotime[keep], stringsAsFactors = FALSE)
  }))
  rownames(expected_axis) <- NULL
  stopifnot(identical(source$engine, "monocle3"),
    identical(source$axis, expected_axis))
  check_reference <- function(case, actual) {
    trajectory_s2_compare_outputs(x$bundle, case,
      list(outputs = actual$outputs), reference)
    expected_params <- reference$cases[[case]]$params
    for (terminal in actual$outputs) {
      for (key in intersect(names(terminal$meta), names(expected_params))) {
        want <- expected_params[[key]]
        same <- if (is.null(want)) {
          length(terminal$meta[[key]]) == 1L &&
            is.na(terminal$meta[[key]])
        } else isTRUE(all.equal(terminal$meta[[key]],
          unlist(want, use.names = FALSE), check.attributes = FALSE))
        if (!same) stop(case, "/", terminal$table_type,
          ": effective parameter differs: ", key, call. = FALSE)
      }
    }
  }
  for (operation in c("de", "smooth", "signature_smooth",
      "geneset_dynamics", "vector_smooth", "signal_dynamics")) {
    params <- .trajectory_s2_model_params(operation)
    params$lineage <- "root:A::lineage:C"
    actual <- .trajectory_s2_run_model_direct(x, operation, params,
      model = cds)
    stopifnot(!inherits(actual$value, "error"))
    check_reference(paste(operation, "C", "full", sep = "__"), actual)
    alternate <- params
    alternate$lineage <- "root:A::lineage:D"
    actual_alt <- .trajectory_s2_run_model_direct(x, operation,
      alternate, model = cds)
    stopifnot(!inherits(actual_alt$value, "error"))
    check_reference(paste(operation, "D", "full", sep = "__"),
      actual_alt)
    session_cells <- x$cells[-(13:16)]
    actual_subset <- .trajectory_s2_run_model_direct(x, operation,
      params, model = cds, effective_cells = session_cells)
    stopifnot(!inherits(actual_subset$value, "error"))
    check_reference(paste(operation, "C", "subset", sep = "__"),
      actual_subset)
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
