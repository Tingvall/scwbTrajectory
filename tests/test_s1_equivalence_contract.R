test_trajectory_old_core_exact_equivalence_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  reference <- trajectory_reference()
  stopifnot(identical(reference$core_commit,
    "909b246d81b5c5b6322be99c1e9f2f69a2615f8b"))
  observed <- trajectory_run_all(x)
  for (operation in names(observed)) {
    expected_outputs <- reference$outputs[[operation]]
    stopifnot(identical(sort(names(observed[[operation]]$outputs)),
      sort(names(expected_outputs))))
    for (output in names(expected_outputs)) {
      actual <- trajectory_identity(x$bundle,
        observed[[operation]]$outputs[[output]])
      stopifnot(identical(actual, expected_outputs[[output]]))
    }
  }
  for (operation in c("trajectory_branch_compare",
      "trajectory_route_smooth", "space_path_score")) {
    expected <- reference$params[[operation]]
    actual <- observed[[operation]]$params[names(expected)]
    stopifnot(isTRUE(all.equal(actual, expected, check.attributes = FALSE)))
  }
  project <- observed$space_path_project_2d$params
  stopifnot(identical(project$project_mode, "local_interpolation"),
    identical(project$approx_points, 3L))
  TRUE
}

test_trajectory_score_prefix_and_parameter_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  run <- function(prefix) wb_session_run_op(wb_session_new(), x$bundle,
    "scwbtrajectory::space_path_score", inputs = list(space = "latent",
      path = x$path, score_prefix = prefix, min_cells = 1L))$run
  plain <- run(NULL); prefixed <- run("route_")
  stopifnot(identical(unique(plain$outputs$path_scores$key),
      c("progress", "distance", "nearest_segment", "nearest_anchor")),
    identical(unique(prefixed$outputs$path_scores$key),
      paste0("route_", c("progress", "distance", "nearest_segment",
        "nearest_anchor"))),
    identical(unname(plain$outputs$path_scores$value),
      unname(prefixed$outputs$path_scores$value)),
    identical(prefixed$params$score_mode, "orthogonal_projection"),
    identical(prefixed$params$distance_metric, "euclidean"))
  TRUE
}
