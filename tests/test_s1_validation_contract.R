test_trajectory_all_operations_adversarial_contract <- function() {
  x <- trajectory_fixture(); on.exit({
    options(scworkbench.operation_pack_v3_output_test_hook = NULL)
    trajectory_close(x)
  }, add = TRUE)
  bad <- list(
    try(wb_session_run_op(wb_session_new(), x$bundle,
      "scwbtrajectory::trajectory_branch_compare", inputs = list(
        matrix = "RNA", pseudotime = x$pseudotime_family,
        branch_membership = x$branches, branch_a = "A", branch_b = "A",
        min_cells = 1L)), silent = TRUE),
    try(wb_session_run_op(wb_session_new(), x$bundle,
      "scwbtrajectory::trajectory_route_smooth", inputs = list(
        matrix = "RNA", pseudotime = x$pseudotime,
        branch_membership = x$route, route_key = "A", features = "g1",
        summary_fun = "median", weight_mode = "probability",
        min_cells_per_bin = 1L, min_cells = 1L)), silent = TRUE),
    try(wb_session_run_op(wb_session_new(), x$bundle,
      "scwbtrajectory::space_path_score", inputs = list(space = "plot",
        path = x$path, min_cells = 1L)), silent = TRUE))
  stopifnot(all(vapply(bad, inherits, logical(1), "try-error")))
  before <- wb_art_list(x$bundle)
  for (operation in c("trajectory_branch_compare", "trajectory_route_smooth",
      "space_path_score", "space_path_project_2d")) {
    options(scworkbench.operation_pack_v3_output_test_hook =
      function(writes, ...) { writes[[1L]] <- NULL; writes })
    failure <- try(trajectory_run_one(x, operation), silent = TRUE)
    options(scworkbench.operation_pack_v3_output_test_hook = NULL)
    stopifnot(inherits(failure, "try-error"),
      identical(before, wb_art_list(x$bundle)))
  }
  TRUE
}

test_trajectory_projection_alignment_and_shape_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  bad_path <- x$path; bad_path$space_id[] <- "plot"
  mismatch <- try(wb_session_run_op(wb_session_new(), x$bundle,
    "scwbtrajectory::space_path_project_2d", inputs = list(path = bad_path,
      source_space = "latent", plot_space = "plot", min_cells = 1L)),
    silent = TRUE)
  bad_plot <- structure(list(kind = "space_2d", space_id = "bad",
    values = x$latent, cell_id = x$cells,
    dimensions = as.integer(dim(x$latent)),
    storage_profile = "space_dense_v1"), class = "wb_primitive")
  dimensions <- try(wb_session_run_op(wb_session_new(), x$bundle,
    "scwbtrajectory::space_path_project_2d", inputs = list(path = x$path,
      source_space = "latent", plot_space = bad_plot, min_cells = 1L)),
    silent = TRUE)
  stopifnot(inherits(mismatch, "try-error"),
    inherits(dimensions, "try-error"))
  TRUE
}

test_trajectory_manifest_semantic_negatives_contract <- function() {
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbTrajectory")
  value <- jsonlite::read_json(manifest, simplifyVector = FALSE)
  check <- function(change) {
    candidate <- change(unserialize(serialize(value, NULL)))
    path <- tempfile(fileext = ".json")
    on.exit(unlink(path, force = TRUE), add = TRUE)
    jsonlite::write_json(candidate, path, auto_unbox = TRUE, null = "null",
      digits = NA)
    try(wb_operation_pack_manifest_validate(path), silent = TRUE)
  }
  bad_prefix <- check(function(x) {
    x$operations[[3L]]$dimension_policies[[1L]]$key_axis$
      prefix_parameter <- "unknown_parameter"; x
  })
  bad_family <- check(function(x) {
    x$operations[[4L]]$outputs[[1L]]$kind <- "obs_table"; x
  })
  stopifnot(inherits(bad_prefix, "try-error"),
    inherits(bad_family, "try-error"))
  TRUE
}
