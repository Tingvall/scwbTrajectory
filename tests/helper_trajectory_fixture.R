trajectory_activate_pack <- function() {
  library(scworkbench)
  if ("scwbtrajectory::space_path_score" %in% wb_list_ops())
    return(invisible(TRUE))
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbTrajectory")
  if (!nzchar(manifest)) stop("scwbTrajectory is not installed")
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
  invisible(TRUE)
}

trajectory_fixture <- function() {
  trajectory_activate_pack()
  root <- tempfile("scwbTrajectory-test-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  cells <- paste0("c", seq_len(24L))
  features <- paste0("g", seq_len(6L))
  wb_define_cells(bundle, cells)
  values <- outer(seq_len(6L), seq_len(24L), function(i, j) i + j %% 5L)
  dimnames(values) <- list(features, cells)
  wb_attach_matrix(bundle, "RNA", Matrix::Matrix(values, sparse = TRUE),
    feature_type = "gene")
  latent <- cbind(d1 = seq(0, 1, length.out = 24L),
    d2 = sin(seq(0, pi, length.out = 24L)),
    d3 = cos(seq(0, pi, length.out = 24L)))
  rownames(latent) <- cells
  plot <- cbind(x = latent[, 1L] + latent[, 2L], y = latent[, 3L])
  rownames(plot) <- cells
  wb_add_space(bundle, "latent", latent)
  wb_add_space(bundle, "plot", plot)
  wb_feature_set_import_list(bundle, "RNA",
    sets = list(selected = features[1:2]), namespace = "trajectory",
    created_by = "test", strict = TRUE, origin = "imported",
    is_readonly = FALSE)
  path <- wb_prim_space_path_family(rep("latent", 12L), rep("p", 12L),
    rep(seq_len(4L), each = 3L), rep(seq_len(3L), 4L),
    as.numeric(t(latent[c(1L, 8L, 16L, 24L), ])),
    meta = list(fit_mode = "linear", checkpoint_n = 2L,
      anchor_centroids = unname(latent[c(1L, 24L), ])))
  pseudotime_family <- wb_prim_cell_vector_family(
    rep(c("A", "B"), each = 24L), rep(cells, 2L),
    c(seq(0, 1, length.out = 24L), seq(1, 0, length.out = 24L)))
  branches <- wb_prim_relation_family(rep(c("A", "B"), each = 12L),
    c(cells[1:12], cells[13:24]), rep(1, 24L))
  pseudotime <- wb_prim_cell_vector(
    setNames(seq(0, 1, length.out = 24L), cells), "continuous")
  route <- wb_prim_relation_family(rep("A", 24L), cells,
    seq(.2, 1, length.out = 24L))
  list(root = root, bundle = bundle, cells = cells, features = features,
    values = values, latent = latent, plot = plot, path = path,
    pseudotime_family = pseudotime_family, branches = branches,
    pseudotime = pseudotime, route = route)
}

trajectory_close <- function(x) {
  try(wb_close_bundle(x$bundle), silent = TRUE)
  unlink(x$root, recursive = TRUE, force = TRUE)
}

trajectory_execute_one <- function(x, operation) {
  run <- function(op, inputs) wb_session_run_op(wb_session_new(), x$bundle,
    op, inputs = inputs)
  switch(operation,
    trajectory_branch_compare = run(
      "scwbtrajectory::trajectory_branch_compare", list(matrix = "RNA",
        pseudotime = x$pseudotime_family,
        branch_membership = x$branches, branch_a = "A", branch_b = "B",
        features = x$features, min_pct = 0, p_adjust_method = "BH",
        min_cells = 1L)),
    trajectory_route_smooth = run(
      "scwbtrajectory::trajectory_route_smooth", list(matrix = "RNA",
        pseudotime = x$pseudotime, branch_membership = x$route,
        route_key = "A", features = x$features[1:2], n_bins = 5L,
        min_cells_per_bin = 1L, weight_mode = "probability", prob_min = 0,
        chunk_size = 128L, min_cells = 1L)),
    space_path_score = run("scwbtrajectory::space_path_score",
      list(space = "latent", path = x$path, score_prefix = "p_",
        min_cells = 1L)),
    space_path_project_2d = run(
      "scwbtrajectory::space_path_project_2d",
      list(path = x$path, source_space = "latent", plot_space = "plot",
        approx_points = 3L, min_cells = 1L)),
    stop("Unknown trajectory operation: ", operation, call. = FALSE))
}

trajectory_run_one <- function(x, operation)
  trajectory_execute_one(x, operation)$run

trajectory_run_all <- function(x) {
  operations <- c("trajectory_branch_compare", "trajectory_route_smooth",
    "space_path_score", "space_path_project_2d")
  stats::setNames(lapply(operations, function(operation)
    trajectory_run_one(x, operation)), operations)
}

trajectory_identity <- function(bundle, value) {
  getFromNamespace(".wb_prov_payload_identity", "scworkbench")(bundle,
    value)
}

trajectory_reference <- function() {
  directory <- if (exists("test_dir", inherits = TRUE))
    get("test_dir", inherits = TRUE) else "tests"
  jsonlite::read_json(file.path(directory, "fixtures",
    "trajectory_old_core_reference.json"), simplifyVector = TRUE)
}
