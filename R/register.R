#' Register the provisional scwbTrajectory operation pack
#'
#' Installation is inert. scworkbench invokes this zero-argument entry point
#' only at its explicit operation-pack activation boundary.
#'
#' @return A validated scworkbench operation-pack registration.
#' @export
register_operation_pack <- function() {
  scworkbench::wb_operation_pack_register(
    system.file("scworkbench", "operation-pack.json",
      package = "scwbTrajectory"),
    list(
      `scwbtrajectory::trajectory_branch_compare` =
        .run_trajectory_branch_compare,
      `scwbtrajectory::trajectory_route_smooth` =
        .run_trajectory_route_smooth,
      `scwbtrajectory::space_path_score` = .run_space_path_score,
      `scwbtrajectory::space_path_project_2d` = .run_space_path_project_2d
    )
  )
}
