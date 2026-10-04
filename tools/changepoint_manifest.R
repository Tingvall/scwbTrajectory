# scwbtrajectory::trajectory_single_changepoint_regions. One interior change
# point per feature along an authenticated pseudotime axis, with the aggregate
# breakpoint density and the caller's directional rank windows. The numeric core
# executes in private Python compute; the ordering, region semantics and every
# terminal write stay in R.
changepoint_region_window <- function() list(
  type = "object",
  properties = list(
    # The parameter-schema profile has no pattern keyword, so the runner enforces
    # the label shape that keeps a region id legible in provenance.
    region_id = list(type = "string", minLength = 1L, maxLength = 64L),
    direction = list(type = "string", enum = list("up", "down")),
    rank_start = list(type = "integer", minimum = 0L),
    rank_end = list(type = "integer", minimum = 0L)),
  required = list("region_id", "direction", "rank_start", "rank_end"),
  additionalProperties = FALSE)

changepoint_params <- list(
  features = strings(),
  feature_set = string(),
  pseudotime_direction = string("ascending", c("ascending", "descending"),
    nullable = FALSE),
  membership_compartment = string(),
  smooth_window = integer(200L, 1L, 1000000L, FALSE),
  window_width = integer(100L, 2L, 1000000L, FALSE),
  window_jump = integer(5L, 1L, 1000000L, FALSE),
  window_min_size = integer(2L, 1L, 1000000L, FALSE),
  # 100 reproduces the notebook's ±100 samples either side of the change point.
  direction_window = integer(100L, 1L, 1000000L, FALSE),
  density_smooth_window = integer(100L, 1L, 1000000L, FALSE),
  density_scale = number(20, 0),
  change_threshold = number(30, 0),
  non_finite_policy = string("error", c("error", "drop"), nullable = FALSE),
  degenerate_bkp_policy = string("error", c("error", "drop"),
    nullable = FALSE),
  chunk_size = integer(256L, 1L, 1000000L, FALSE),
  # Omitted means the whole axis in both directions. A literal empty array is
  # not representable: the core parameter catalog rejects it as an empty list
  # before the runner is entered, so the default must be reachable by omission.
  region_windows = list(type = c("array", "null"),
    items = changepoint_region_window()))

changepoint_output_names <- c("feature_single_changepoints", "trend_regions",
  "trend_density", "run_summary")

# Boundary rank b sits before the cell at 0-based rank b, so the published
# density is shifted by one from the notebook's raw histogram: rank b holds the
# breakpoint count of boundary b, matching the directional delta windows.
changepoint_density_policy <- list(
  output = "trend_density", policy = "authenticated_output_axes_v1",
  shape = "bounded_grouped_cell_key_cartesian",
  cell_axis = list(source = "effective_cells_v1", input = "matrix"),
  key_axis = list(source = "declared_bounded_key_groups_v1",
    required_keys = list("all", "up", "down"), optional_groups = list(),
    selection_owner = "operation_pack_v1", group_order = "manifest_order_v1"),
  row_order = "key_major_then_cell_v1",
  missing_values = "canonical_null_allowed")

changepoint_operation <- list(
  name = "trajectory_single_changepoint_regions",
  qualified_name = "scwbtrajectory::trajectory_single_changepoint_regions",
  title = "Single change-point regions along pseudotime",
  description = paste(
    "Detect one interior L1 change point per matrix feature along an",
    "authenticated pseudotime axis in private Python compute, classify its",
    "direction from the mean signed difference across the boundary, aggregate",
    "breakpoint counts into an aligned density, and publish the features that",
    "exceed the density threshold inside each caller-declared directional rank",
    "window. Rank semantics are half-open and explicit; the notebook raw",
    "density convention is recorded in run_summary."),
  runner = ".run_trajectory_single_changepoint_regions",
  inputs = list(input("matrix", "matrix"), input("pseudotime", "cell_vector"),
    input("membership", "cell_vector", required = FALSE)),
  params_schema = list(type = "object", properties = changepoint_params,
    additionalProperties = FALSE),
  registry_reads = list(feature_registry),
  dimension_policies = list(
    table_policy("feature_single_changepoints"),
    table_policy("trend_regions"), changepoint_density_policy,
    table_policy("run_summary")),
  outputs = list(
    output("feature_single_changepoints", "feature_table",
      "single_changepoint", "trajectory_single_changepoints"),
    output("trend_regions", "obs_table", "single_changepoint",
      "trajectory_single_changepoint_regions"),
    output("trend_density", "cell_vector_family", "single_changepoint"),
    output("run_summary", "obs_table", "single_changepoint",
      "trajectory_single_changepoint_summary")),
  dependencies = list(dependency("scworkbench", "wb_prim_obs_table")),
  promotion_closure = closure(changepoint_output_names),
  replay_equivalence = stats::setNames(
    as.list(rep("exact_canonical_value", length(changepoint_output_names))),
    changepoint_output_names),
  typed_slice_reads = list(matrix_read("chunk_size")),
  retained_parameters = list(
    retained_row("matrix_id", "run_summary", "matrix_id"),
    retained_row("n_cells_analysed", "run_summary", "n_cells_analysed"),
    retained_row("n_features_analysed", "run_summary", "n_features_analysed"),
    retained_row("significant_rank_count", "run_summary",
      "significant_rank_count"),
    retained_parameter("feature_set"),
    retained_parameter("pseudotime_direction", "effective_parameter_value",
      required = TRUE),
    retained_parameter("smooth_window", "effective_parameter_value",
      required = TRUE),
    retained_parameter("window_width", "effective_parameter_value",
      required = TRUE),
    retained_parameter("window_jump", "effective_parameter_value",
      required = TRUE),
    retained_parameter("window_min_size", "effective_parameter_value",
      required = TRUE),
    retained_parameter("direction_window", "effective_parameter_value",
      required = TRUE),
    retained_parameter("density_smooth_window", "effective_parameter_value",
      required = TRUE),
    retained_parameter("density_scale", "effective_parameter_value",
      required = TRUE),
    retained_parameter("change_threshold", "effective_parameter_value",
      required = TRUE),
    retained_parameter("non_finite_policy", "effective_parameter_value",
      required = TRUE),
    retained_parameter("degenerate_bkp_policy", "effective_parameter_value",
      required = TRUE),
    retained_parameter("chunk_size", "effective_parameter_value")),
  public_interface = interface(c("trajectory", "changepoint", "pseudotime"),
    c("matrix", "pseudotime", "membership"), list(
      c("features", "string_vec", "Optional explicit feature IDs."),
      c("feature_set", "string", "Optional matrix feature-set name."),
      c("pseudotime_direction", "choice", "Ascending or descending axis."),
      c("membership_compartment", "string",
        "Optional cell-vector compartment restricting the axis."),
      c("smooth_window", "int", "Centred pre-smoothing window."),
      c("window_width", "int", "Sliding-window cost width."),
      c("window_jump", "int", "Breakpoint quantisation step."),
      c("window_min_size", "int", "Minimum segment length."),
      c("direction_window", "int", "Samples averaged either side."),
      c("density_smooth_window", "int", "Breakpoint histogram smoothing."),
      c("density_scale", "dbl", "Multiplier on the smoothed histogram."),
      c("change_threshold", "dbl",
        "Density value at which a rank is significant."),
      c("non_finite_policy", "choice", "Error or drop on non-finite input."),
      c("degenerate_bkp_policy", "choice",
        "Error or drop when no interior breakpoint exists."),
      c("chunk_size", "int", "Authenticated feature chunk size."),
      c("region_windows", "list",
        "Objects of region_id, direction, rank_start, rank_end; ranks are half-open [rank_start, rank_end). Omit to analyse the whole axis in both directions.")),
    changepoint_output_names))

changepoint_operations <- list(changepoint_operation)
