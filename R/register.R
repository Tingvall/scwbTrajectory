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
      `scwbtrajectory::space_path_project_2d` = .run_space_path_project_2d,
      `scwbtrajectory::space_path_define` = .run_space_path_define,
      `scwbtrajectory::trajectory_segment_from_artifacts` =
        .run_trajectory_segment_from_artifacts,
      `scwbtrajectory::trajectory_space_path_segment` =
        .run_trajectory_space_path_segment,
      `scwbtrajectory::trajectory_segment_join` =
        .run_trajectory_segment_join,
      `scwbtrajectory::trajectory_segment_path_join` =
        .run_trajectory_segment_path_join,
      `scwbtrajectory::trajectory_de` = .run_trajectory_de_model,
      `scwbtrajectory::trajectory_smooth` = .run_trajectory_smooth_model,
      `scwbtrajectory::trajectory_signature_smooth` =
        .run_trajectory_signature_smooth_model,
      `scwbtrajectory::trajectory_geneset_dynamics` =
        .run_trajectory_geneset_dynamics_model,
      `scwbtrajectory::trajectory_vector_smooth` =
        .run_trajectory_vector_smooth_model,
      `scwbtrajectory::trajectory_signal_dynamics` =
        .run_trajectory_signal_dynamics_model,
      `scwbtrajectory::trajectory_geneset_dynamics_explicit` =
        .run_trajectory_geneset_dynamics_explicit,
      `scwbtrajectory::trajectory_signal_dynamics_explicit` =
        .run_trajectory_signal_dynamics_explicit,
      `scwbtrajectory::slingshot` = .run_slingshot,
      `scwbtrajectory::slingshot_fit` = .run_slingshot_fit,
      `scwbtrajectory::slingshot_lineage_weights` =
        .run_slingshot_lineage_weights,
      `scwbtrajectory::slingshot_project` = .run_slingshot_project,
      `scwbtrajectory::monocle3_fit` = .run_monocle3_fit,
      `scwbtrajectory::monocle3_project` = .run_monocle3_project,
      `scwbtrajectory::palantir_fit` = .run_palantir_fit,
      `scwbtrajectory::palantir_project` = .run_palantir_project,
      `scwbtrajectory::trajectory_single_changepoint_regions` =
        .run_trajectory_single_changepoint_regions
    )
  )
}
