occupancy_params <- list(
  pseudotime_key = list(type = "string", minLength = 1L),
  acceptance_output_key = list(type = "string", minLength = 1L),
  n_bins = integer(10L, 2L, 1000L, FALSE),
  analysis_mode = string("descriptive", "descriptive", FALSE))
occupancy_outputs <- c("bin_occupancy", "acceptance_accounting",
  "run_summary")
acceptance_input <- input("mapping_acceptance", "obs_table")
acceptance_input$effective_cell_coverage <- "reference_only"
occupancy_provenance <- lapply(c("pseudotime", "group_id", "sample_id",
    "mapping_acceptance"), function(name) list(
  name = paste0(name, "_ancestry"),
  category = "authenticated_core_provenance_snapshot",
  kind = "artifact_ancestry", reader = "read_artifact_ancestry",
  input = name, activation = "when_resolved_from_artifact",
  lifecycle_policy = "require_current"))

occupancy_operation <- list(
  name = "trajectory_occupancy",
  qualified_name = "scwbtrajectory::trajectory_occupancy",
  title = "Support-aware Trajectory Occupancy",
  description = paste(
    "Describe accepted cell occupancy over deterministic global pseudotime",
    "bins while retaining complete accepted and withheld accounting."),
  runner = ".run_trajectory_occupancy",
  inputs = list(input("pseudotime", "cell_vector_family"),
    input("group_id", "cell_vector"), input("sample_id", "cell_vector", FALSE),
    acceptance_input),
  params_schema = list(type = "object", properties = occupancy_params,
    required = list("pseudotime_key", "acceptance_output_key"),
    additionalProperties = FALSE),
  provenance_reads = occupancy_provenance,
  dimension_policies = lapply(occupancy_outputs, table_policy),
  outputs = list(
    output("bin_occupancy", "obs_table", "trajectory_occupancy",
      "trajectory_occupancy_bins"),
    output("acceptance_accounting", "obs_table", "trajectory_occupancy",
      "trajectory_occupancy_acceptance_accounting"),
    output("run_summary", "obs_table", "trajectory_occupancy",
      "trajectory_occupancy_summary")),
  dependencies = list(dependency("scworkbench", "wb_prim_obs_table",
    ">=0.1.2,<0.2.0")),
  promotion_closure = closure(occupancy_outputs),
  replay_equivalence = equivalence(occupancy_outputs),
  retained_parameters = c(lapply(names(occupancy_params), function(name)
    retained_parameter(name, "effective_parameter_value", required = TRUE)),
    list(retained_row("n_cells_input", "run_summary", "n_cells_input"),
      retained_row("n_cells_accepted", "run_summary", "n_cells_accepted"),
      retained_row("n_cells_withheld", "run_summary", "n_cells_withheld"),
      retained_row("n_groups", "run_summary", "n_groups"),
      retained_row("n_samples", "run_summary", "n_samples"))),
  public_interface = interface(c("trajectory", "occupancy", "descriptive",
      "support_aware"),
    c("pseudotime", "group_id", "sample_id", "mapping_acceptance"), list(
      c("pseudotime_key", "string", "Selected raw pseudotime key."),
      c("acceptance_output_key", "string",
        "Selected mapping-acceptance output key."),
      c("n_bins", "int", "Number of global equal-width pseudotime bins."),
      c("analysis_mode", "choice",
        "Descriptive-only analysis; inference is not admitted.")),
    occupancy_outputs))

occupancy_operations <- list(occupancy_operation)
