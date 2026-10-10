# Contract tests for scwbtrajectory::trajectory_single_changepoint_regions.
#
# The numeric core runs in private Python compute, so these tests assert the
# published contract only: rank translation, direction gating (a region is gated by
# the density of its own direction), the aligned density convention, the
# authenticated density-family shape, missing-value handling under both
# non_finite_policy values, the boundary-pseudotime summary statistics, the JSON
# encoding of region member ids, chunk-size invariance and the documented rejection
# paths. Fixture values are planted steps on a pseudotime axis, so every expected
# rank is exact.
#
# Reproduce locally inside the qualified image:
#   TMPDIR=/tmp Rscript tests/run_contract_tests.R
cp_test <- list(
  n = 30L,
  params = list(smooth_window = 3L, window_width = 8L, window_jump = 1L,
    window_min_size = 2L, direction_window = 3L, density_smooth_window = 1L,
    density_scale = 1, change_threshold = 2, chunk_size = 256L,
    degenerate_bkp_policy = "drop", non_finite_policy = "drop"))

cp_step <- function(n, at, low, high) {
  value <- rep(low, n)
  value[seq.int(at, n)] <- high
  value
}

# up1/up2 step up before c16, dn1 steps down before c16, late steps up before
# c26, flat has no breakpoint. density_smooth_window = 1 and density_scale = 1
# make the published density exactly the breakpoint count per cell rank.
cp_fixture <- function(ties = FALSE, n = cp_test$n,
    features = c("up1", "up2", "dn1", "late", "flat")) {
  trajectory_activate_pack()
  root <- tempfile("scwbTrajectory-cp-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  cells <- sprintf("c%02d", seq_len(n))
  time <- (seq_len(n) - 1L) / (n - 1L)
  if (ties) time[c(20L, 21L, 22L)] <- time[20L]
  planted <- list(up1 = cp_step(n, 16L, 0, 9), up2 = cp_step(n, 16L, 1, 7),
    dn1 = cp_step(n, 16L, 9, 0), late = cp_step(n, 26L, 2, 6),
    flat = rep(4, n))
  values <- do.call(rbind, planted[features])
  dimnames(values) <- list(features, cells)
  wb_define_cells(bundle, cells)
  wb_attach_matrix(bundle, "RNA", Matrix::Matrix(values, sparse = TRUE),
    feature_type = "gene")
  list(root = root, bundle = bundle, cells = cells, features = features,
    pseudotime = wb_prim_cell_vector(setNames(time, cells), "continuous"),
    compartment = wb_prim_cell_vector(setNames(ifelse(seq_len(n) <= 25L, "in",
      "out"), cells), "categorical"))
}

cp_run <- function(x, params = list()) {
  inputs <- c(list(matrix = "RNA", pseudotime = x$pseudotime),
    modifyList(cp_test$params, params))
  if (is.character(inputs$membership_compartment))
    inputs$membership <- x$compartment
  outputs <- wb_session_run_op(wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_single_changepoint_regions",
    inputs = inputs)$run$outputs
  stats::setNames(lapply(names(outputs), function(name) {
    terminal <- outputs[[name]]
    if (identical(name, "trend_density"))
      list(key = terminal$key, cell_id = terminal$cell_id,
        value = terminal$value)
    else terminal$table
  }), names(outputs))
}

cp_failed <- function(x, params = list())
  inherits(try(cp_run(x, params), silent = TRUE), "try-error")

cp_window <- function(id, direction, start, end)
  list(region_id = id, direction = direction, rank_start = start,
    rank_end = end)

# Published density values for one family key, in cell-axis order.
cp_key <- function(density, key) density$value[density$key == key]

# Fixture with caller-chosen planted vectors on the same evenly spaced axis, for a
# test that must state its own expected density, ranks and members instead of
# inheriting the shared fixture's feature set.
cp_matrix_fixture <- function(planted, n = cp_test$n) {
  trajectory_activate_pack()
  root <- tempfile("scwbTrajectory-cp-matrix-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  cells <- sprintf("c%02d", seq_len(n))
  time <- (seq_len(n) - 1L) / (n - 1L)
  values <- do.call(rbind, planted)
  dimnames(values) <- list(names(planted), cells)
  wb_define_cells(bundle, cells)
  wb_attach_matrix(bundle, "RNA", Matrix::Matrix(values, sparse = TRUE),
    feature_type = "gene")
  list(root = root, bundle = bundle, cells = cells, features = names(planted),
    pseudotime = wb_prim_cell_vector(setNames(time, cells), "continuous"))
}

# Published density value for one family key at one 0-based cell rank.
cp_at <- function(density, key, rank) cp_key(density, key)[rank + 1L]

# Region member ids, parsed back out of the JSON column. An empty array parses to
# an empty list, so the result is normalized to character(0).
cp_members <- function(regions, direction) {
  parsed <- jsonlite::fromJSON(regions$feature_ids_json[regions$direction ==
    direction])
  if (!length(parsed)) return(character(0))
  as.character(unlist(parsed, use.names = FALSE))
}

# A published double travels through the payload's canonical decimal form, so it
# equals the input value to about 15 significant digits, not bit for bit. Ranks,
# counts and identifiers are asserted with identical() everywhere else.
cp_near <- function(value, target)
  isTRUE(all.equal(as.numeric(value), as.numeric(target), tolerance = 1e-12))

test_trajectory_single_changepoint_default_regions_contract <- function() {
  x <- cp_fixture(); on.exit(trajectory_close(x), add = TRUE)
  out <- cp_run(x)
  features <- out$feature_single_changepoints
  stopifnot(identical(names(out), c("feature_single_changepoints",
      "trend_regions", "trend_density", "run_summary")),
    # Omitted region_windows means the whole axis once per direction.
    identical(nrow(out$trend_regions), 2L),
    identical(out$trend_regions$region_id, c("whole_axis", "whole_axis")),
    identical(out$trend_regions$direction, c("up", "down")),
    identical(out$trend_regions$region_source, rep("default_whole_axis", 2L)),
    identical(out$trend_regions$rank_start, c(0L, 0L)),
    identical(out$trend_regions$rank_end, c(30L, 30L)),
    identical(out$trend_regions$first_included_cell_id, c("c01", "c01")),
    identical(out$trend_regions$last_included_cell_id, c("c30", "c30")),
    # Direction-matched gating: at rank 15 the up density is 2 (up1, up2) and the
    # down density 1 (dn1), so with change_threshold 2 the up region holds both up
    # features and the down region is empty - although the combined density there
    # is 3, which is what used to gate both.
    identical(out$trend_regions$feature_count, c(2L, 0L)),
    identical(out$run_summary$region_window_count, 2L),
    # Planted boundaries: rank 15 for the three steps before c16, rank 25 for
    # the step before c26, none for flat (dropped under the drop policy).
    identical(features$feature_id, x$features),
    identical(features$boundary_rank, c(15L, 15L, 15L, 25L, NA_integer_)),
    identical(features$reason,
      c(NA_character_, NA_character_, NA_character_, NA_character_,
        "degenerate")),
    identical(features$direction, c("up", "up", "down", "up", NA_character_)),
    # Boundary rank b sits before the cell at 0-based rank b.
    identical(features$boundary_before_cell_id[1L], "c15"),
    identical(features$boundary_after_cell_id[1L], "c16"),
    cp_near(features$change_point_pseudotime[1L], 15 / 29),
    # The delta window is [b - w, b) versus [b, b + w) in ranks.
    identical(features$delta_window_first_rank[1L], 12L),
    identical(features$delta_window_last_rank[1L], 17L),
    identical(sign(features$delta), c(1, 1, -1, 1, NA_real_)))
  TRUE
}

test_trajectory_single_changepoint_aligned_density_contract <- function() {
  x <- cp_fixture(); on.exit(trajectory_close(x), add = TRUE)
  out <- cp_run(x)
  density <- cp_key(out$trend_density, "all")
  up <- cp_key(out$trend_density, "up")
  down <- cp_key(out$trend_density, "down")
  # D1: a breakpoint at boundary rank b publishes at cell rank b, not b - 1.
  # Cell c16 carries 0-based rank 15, so it holds the count for boundary 15.
  stopifnot(length(density) == 30L,
    identical(density[c(15L, 16L, 17L)], c(0, 3, 0)),
    identical(which(density > 0) - 1L, c(15L, 25L)),
    identical(up[c(16L, 26L)], c(2, 1)),
    identical(down[c(15L, 16L)], c(0, 1)),
    identical(out$run_summary$aligned_peak_rank_all, 15L),
    # The notebook's raw convention keeps its own peak one slot lower.
    identical(out$run_summary$raw_peak_rank_all, 14L),
    identical(out$run_summary$density_alignment,
      "aligned_boundary_at_cell_rank"),
    identical(out$run_summary$notebook_raw_density_convention,
      "slot_j_carries_boundary_j_plus_1"),
    identical(out$run_summary$density_shift_formula,
      "aligned[j] = raw[j - 1]; aligned[0] = 0"))
  # The threshold applies to the aggregate density, never to one feature.
  stopifnot(identical(which(density >= cp_test$params$change_threshold) - 1L,
      15L),
    identical(out$run_summary$significant_rank_count, 1L),
    identical(out$feature_single_changepoints$density_significant,
      c(TRUE, TRUE, TRUE, FALSE, NA)))
  TRUE
}

test_trajectory_single_changepoint_region_rank_translation_contract <- function() {
  x <- cp_fixture(); on.exit(trajectory_close(x), add = TRUE)
  out <- cp_run(x, list(region_windows = list(
    cp_window("strict", "up", 10L, 20L), cp_window("one", "up", 15L, 16L),
    cp_window("empty", "up", 12L, 12L), cp_window("head", "down", 0L, 1L),
    cp_window("tail", "down", 29L, 30L), cp_window("all", "up", 0L, 30L),
    cp_window("miss", "down", 25L, 26L))))
  regions <- out$trend_regions
  stopifnot(identical(regions$region_id,
      c("strict", "one", "empty", "head", "tail", "all", "miss")),
    # A declared [start, end) covers ranks start..end-1: the cell whose 0-based
    # rank equals start through the cell whose rank equals end - 1.
    identical(regions$rank_count, c(10L, 1L, 0L, 1L, 1L, 30L, 1L)),
    identical(regions$first_included_cell_id,
      c("c11", "c16", NA, "c01", "c30", "c01", "c26")),
    identical(regions$last_included_cell_id,
      c("c20", "c16", NA, "c01", "c30", "c30", "c26")),
    identical(regions$first_included_rank, c(10L, 15L, NA, 0L, 29L, 0L, 25L)),
    identical(regions$last_included_rank, c(19L, 15L, NA, 0L, 29L, 29L, 25L)),
    cp_near(regions$first_included_pseudotime[1L], 10 / 29),
    is.na(regions$first_included_pseudotime[3L]),
    # A zero-length window publishes the empty JSON array, not [""].
    identical(regions$feature_ids_json[3L], "[]"),
    identical(regions$feature_count[c(1L, 2L, 6L)], c(2L, 2L, 2L)),
    identical(regions$feature_count[c(3L, 4L, 5L, 7L)], rep(0L, 4L)),
    all(regions$feature_count == vapply(regions$feature_ids_json,
      function(text) length(jsonlite::fromJSON(text)), integer(1))))
  # A one-cell window over the significant rank keeps both up features; the
  # boundary at rank 25 is present but below the aggregate threshold.
  stopifnot(identical(jsonlite::fromJSON(regions$feature_ids_json[2L]),
      c("up1", "up2")),
    identical(jsonlite::fromJSON(regions$feature_ids_json[6L]),
      c("up1", "up2")))
  # Direction gating: every member carries the window's own direction.
  features <- out$feature_single_changepoints
  direction_of <- stats::setNames(features$direction, features$feature_id)
  for (row in seq_len(nrow(regions))) {
    ids <- jsonlite::fromJSON(regions$feature_ids_json[row])
    if (!length(ids)) next
    stopifnot(identical(unique(unname(direction_of[ids])), regions$direction[row]))
  }
  TRUE
}

test_trajectory_single_changepoint_axis_order_contract <- function() {
  # A tie group makes the ordering test non-trivial. Descending must re-sort on
  # the negated pseudotime with ties broken by ascending cell id; that is not
  # the ascending order read backwards.
  x <- cp_fixture(ties = TRUE); on.exit(trajectory_close(x), add = TRUE)
  window <- list(region_windows = list(cp_window("tie", "up", 8L, 10L)))
  asc <- cp_run(x, window)$trend_regions
  desc <- cp_run(x, modifyList(window,
    list(pseudotime_direction = "descending")))$trend_regions
  stopifnot(
    # Ascending ranks 8 and 9 are c09 and c10.
    identical(asc$first_included_cell_id, "c09"),
    identical(asc$last_included_cell_id, "c10"),
    # Descending: c30..c23 take ranks 0..7, then the tie group takes ranks 8,
    # 9, 10 in ascending cell-id order (c20, c21, c22). The ascending order
    # reversed would publish c22 then c21 here, so this asserts a real re-sort.
    identical(desc$first_included_cell_id, "c20"),
    identical(desc$last_included_cell_id, "c21"))
  # The mirrored boundary translation is the decisive evidence: descending
  # analysis puts the step between c16 and c15, ascending between c15 and c16.
  asc_features <- cp_run(x)$feature_single_changepoints
  desc_features <- cp_run(x,
    list(pseudotime_direction = "descending"))$feature_single_changepoints
  stopifnot(identical(asc_features$boundary_before_cell_id[1L], "c15"),
    identical(asc_features$boundary_after_cell_id[1L], "c16"),
    identical(desc_features$boundary_before_cell_id[1L], "c16"),
    identical(desc_features$boundary_after_cell_id[1L], "c15"),
    identical(desc_features$boundary_rank[1L], 15L))
  TRUE
}

test_trajectory_single_changepoint_membership_contract <- function() {
  x <- cp_fixture(); on.exit(trajectory_close(x), add = TRUE)
  out <- cp_run(x, list(membership_compartment = "in"))
  density <- cp_key(out$trend_density, "all")
  features <- out$feature_single_changepoints
  stopifnot(identical(out$run_summary$n_cells_analysed, 25L),
    identical(out$run_summary$membership_compartment, "in"),
    # The family stays the full key x effective-cell product; cells outside
    # the compartment keep a canonical null in their own slot.
    identical(length(density), 30L),
    identical(which(is.na(density)), 26:30),
    identical(unique(out$trend_density$cell_id), x$cells),
    identical(out$trend_regions$rank_end, c(25L, 25L)),
    # Ranks refer to the restricted axis: the step before c16 stays at rank 15,
    # the step before c26 leaves the 25-cell axis with no interior breakpoint.
    identical(features$boundary_rank[features$feature_id == "up1"], 15L),
    identical(features$reason[features$feature_id == "late"], "degenerate"),
    # A compartment name needs the membership vector and vice versa.
    cp_failed(x, list(membership_compartment = "absent")))
  TRUE
}

test_trajectory_single_changepoint_chunk_invariance_contract <- function() {
  x <- cp_fixture(); on.exit(trajectory_close(x), add = TRUE)
  n <- length(x$cells)
  sizes <- c(1L, 2L, 3L, 256L, n, n + 7L)
  # Five selected features, so the feature-axis chunk plan is known.
  partitions <- c(5L, 3L, 2L, 1L, 1L, 1L)
  reference <- cp_run(x, list(chunk_size = sizes[1L]))
  # run_summary echoes the chunk plan; nothing else may move.
  plan <- c("chunk_size", "partition_count")
  stripped <- function(out) c(
    lapply(names(out), function(name) if (identical(name, "run_summary"))
      out[[name]][, setdiff(names(out[[name]]), plan), drop = FALSE]
    else out[[name]]), list(density = out$trend_density$value))
  for (index in seq_along(sizes)) {
    out <- cp_run(x, list(chunk_size = sizes[index]))
    stopifnot(identical(stripped(out), stripped(reference)),
      identical(out$run_summary$partition_count, partitions[index]))
  }
  # Two runs of the same parameters publish identical payloads.
  stopifnot(identical(stripped(cp_run(x, list(chunk_size = sizes[1L]))),
    stripped(reference)))
  TRUE
}

test_trajectory_single_changepoint_probe_identity_contract <- function() {
  x <- cp_fixture(); on.exit(trajectory_close(x), add = TRUE)
  narrow <- cp_run(x)$run_summary
  wide <- cp_run(x, list(window_width = 60L, window_min_size = 3L,
    window_jump = 7L))$run_summary
  stopifnot(
    # The probe signal must be longer than the analysis window.
    identical(narrow$probe_n, 200L), identical(wide$probe_n, 240L),
    identical(narrow$probe_planted_boundary, 100L),
    identical(wide$probe_planted_boundary, 120L),
    identical(narrow$probe_recovered_boundary, 100L),
    abs(wide$probe_recovered_boundary - 120L) <= 7L,
    identical(narrow$probe_terminal_sample_count, narrow$probe_n),
    identical(wide$probe_terminal_sample_count, wide$probe_n),
    identical(narrow$probe_signal, "deterministic_planted_step"),
    # Recorded runtime identity, not a version gate.
    nzchar(narrow$python_version), nzchar(narrow$changepoint_core_sha256),
    identical(narrow$changepoint_core_sha256, wide$changepoint_core_sha256),
    identical(narrow$ruptures_version, wide$ruptures_version))
  TRUE
}

test_trajectory_single_changepoint_rejection_contract <- function() {
  x <- cp_fixture(); on.exit({
    options(scwbTrajectory.changepoint_python = NULL)
    trajectory_close(x)
  }, add = TRUE)
  bad <- list(
    # The half-open rank domain is [0, n].
    rank_high = list(region_windows = list(cp_window("b", "up", 0L, 31L))),
    reversed = list(region_windows = list(cp_window("b", "up", 10L, 9L))),
    negative = list(region_windows = list(cp_window("b", "up", -1L, 9L))),
    # An empty array is not representable in the core parameter catalog, so the
    # whole-axis default is reached by omission instead.
    empty_array = list(region_windows = list()),
    duplicate = list(region_windows = list(cp_window("d", "up", 0L, 4L),
      cp_window("d", "up", 5L, 9L))),
    illegible_id = list(region_windows = list(cp_window("_b", "up", 0L, 4L))),
    missing_direction = list(region_windows = list(
      list(region_id = "b", rank_start = 0L, rank_end = 4L))),
    unknown_direction = list(region_windows = list(
      list(region_id = "b", direction = "sideways", rank_start = 0L,
        rank_end = 4L))),
    width_too_small = list(window_width = 1L),
    min_size_not_below_width = list(window_width = 8L, window_min_size = 8L),
    degenerate_error = list(degenerate_bkp_policy = "error", features = "flat"),
    no_feature_selected = list(features = "not-a-feature"))
  not_rejected <- names(bad)[!vapply(bad, function(params)
    cp_failed(x, params), logical(1))]
  if (length(not_rejected)) cat("NOT REJECTED:",
    paste(not_rejected, collapse = ", "), "\n")
  stopifnot(!length(not_rejected))
  # The same pair of windows is accepted once each direction is declared.
  stopifnot(!cp_failed(x, list(region_windows = list(cp_window("d", "up", 0L,
      4L), cp_window("d", "down", 0L, 4L)))))
  options(scwbTrajectory.changepoint_python = "/nonexistent/python3")
  stopifnot(cp_failed(x))
  options(scwbTrajectory.changepoint_python = NULL)
  stopifnot(!cp_failed(x))
  TRUE
}

test_trajectory_single_changepoint_non_finite_axis_contract <- function() {
  x <- cp_fixture(); on.exit(trajectory_close(x), add = TRUE)
  time <- (seq_along(x$cells) - 1L) / (length(x$cells) - 1L)
  time[5L] <- NA_real_
  broken <- x
  broken$pseudotime <- wb_prim_cell_vector(setNames(time, x$cells),
    "continuous")
  out <- cp_run(broken)
  density <- cp_key(out$trend_density, "all")
  stopifnot(identical(out$run_summary$n_cells_input, 30L),
    identical(out$run_summary$n_cells_analysed, 29L),
    identical(out$run_summary$n_cells_dropped_non_finite_pseudotime, 1L),
    identical(which(is.na(density)), 5L),
    identical(out$trend_regions$rank_end, c(29L, 29L)),
    # The dropped cell keeps its slot in the authenticated cell axis.
    identical(unique(out$trend_density$cell_id), x$cells),
    identical(out$trend_density$key, rep(c("all", "up", "down"), each = 30L)))
  TRUE
}

test_trajectory_single_changepoint_manifest_contract <- function() {
  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbTrajectory"),
    simplifyVector = FALSE)
  target <- which(vapply(manifest$operations, function(entry) identical(
    entry$qualified_name,
    "scwbtrajectory::trajectory_single_changepoint_regions"), logical(1)))
  stopifnot(identical(length(target), 1L))
  operation <- manifest$operations[[target]]
  policy <- Filter(function(entry) identical(entry$output, "trend_density"),
    operation$dimension_policies)[[1L]]
  stopifnot(
    identical(vapply(operation$outputs, function(entry) entry$name,
      character(1)), c("feature_single_changepoints", "trend_regions",
      "trend_density", "run_summary")),
    identical(vapply(operation$outputs, function(entry) entry$promotion_group,
      character(1)), rep("single_changepoint", 4L)),
    identical(vapply(operation$inputs, function(entry) entry$name,
      character(1)), c("matrix", "pseudotime", "membership")),
    identical(operation$inputs[[3L]]$required, FALSE),
    # The density family is a declared key set times the effective cell axis;
    # this core build has no prefix parameter, so the keys must be bounded.
    identical(policy$shape, "bounded_grouped_cell_key_cartesian"),
    identical(policy$key_axis$source, "declared_bounded_key_groups_v1"),
    identical(unlist(policy$key_axis$required_keys), c("all", "up", "down")),
    identical(policy$row_order, "key_major_then_cell_v1"),
    identical(policy$missing_values, "canonical_null_allowed"),
    # A literal empty region array is unrepresentable, so the parameter has no
    # default value: whole-axis behaviour is reached by omission.
    is.null(operation$params_schema$properties$region_windows$default),
    identical(operation$params_schema$properties$region_windows$type,
      list("array", "null")))
  # Negatives: the density policy shape and null policy are both enforced.
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  tamper <- function(change) {
    broken <- unserialize(serialize(manifest, NULL))
    index <- which(vapply(broken$operations[[target]]$dimension_policies,
      function(entry) identical(entry$output, "trend_density"), logical(1)))
    broken$operations[[target]]$dimension_policies[[index]] <- change(
      broken$operations[[target]]$dimension_policies[[index]])
    jsonlite::write_json(broken, path, auto_unbox = TRUE, null = "null",
      digits = NA)
    try(wb_operation_pack_manifest_validate(path), silent = TRUE)
  }
  stopifnot(inherits(tamper(function(entry) {
    entry$shape <- "cell_key_cartesian"; entry
  }), "try-error"),
    inherits(tamper(function(entry) {
      entry$missing_values <- "silently_zero"; entry
    }), "try-error"),
    inherits(tamper(function(entry) {
      entry$row_order <- "cell_major_then_key"; entry
    }), "try-error"))
  TRUE
}

# The method defaults are the notebook's, published verbatim rather than tuned:
# smoothing 200, detection width 100, jump 5, min_size 2, direction window 100,
# density smoothing 100, scale 20, threshold 30, both policies "error", chunk
# 256. Nothing here is inferred from a run - a default is only a default if the
# manifest says so.
test_trajectory_single_changepoint_published_defaults_contract <- function() {
  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbTrajectory"),
    simplifyVector = FALSE)
  operation <- Filter(function(entry) identical(entry$qualified_name,
    "scwbtrajectory::trajectory_single_changepoint_regions"),
    manifest$operations)[[1L]]
  defaults <- lapply(operation$params_schema$properties,
    function(entry) entry$default)
  notebook <- list(smooth_window = 200, window_width = 100, window_jump = 5,
    window_min_size = 2, direction_window = 100, density_smooth_window = 100,
    density_scale = 20, change_threshold = 30, non_finite_policy = "error",
    degenerate_bkp_policy = "error", chunk_size = 256,
    pseudotime_direction = "ascending")
  # jsonlite hands back whole JSON numbers as integers, so the value is compared
  # numerically and the string defaults exactly.
  for (name in names(notebook)) {
    got <- defaults[[name]]
    want <- notebook[[name]]
    stopifnot(if (is.character(want)) identical(got, want) else
      identical(as.numeric(got), as.numeric(want)))
  }
  # Selection parameters and the region list carry no default: an empty array is
  # unrepresentable in the core parameter catalog, so the whole-axis pair of rows
  # is reached by omission (see the plan's correction C-4).
  without_default <- c("features", "feature_set", "membership_compartment",
    "region_windows")
  stopifnot(identical(setdiff(names(defaults), names(notebook)),
      without_default),
    all(vapply(defaults[without_default], is.null, logical(1))))
  TRUE
}

# With the notebook's own geometry - a 200-wide centred smoother feeding a 100-wide
# L1 window - a noiseless step is smoothed into a 200-sample ramp, and
# Window(model = "l1", width = 100) has nothing left to find. Executed on three
# amplitudes of the same step it returned no interior breakpoint for two of them
# (predict gave only the terminal sample count) and a spurious rank 270 for the
# third. The notebook's silent case is the first: it records [n], the difference
# is NaN and the feature is classified "down". Here the shipped default fails the
# run, and the drop policy publishes the feature as unanalysed, with no boundary,
# no direction and no invented significance.
test_trajectory_single_changepoint_degenerate_policy_contract <- function() {
  trajectory_activate_pack()
  n <- 400L
  root <- tempfile("scwbTrajectory-cp-degenerate-")
  bundle <- wb_init_bundle(root, overwrite = TRUE)
  on.exit(trajectory_close(list(root = root, bundle = bundle)), add = TRUE)
  cells <- sprintf("c%03d", seq_len(n))
  wb_define_cells(bundle, cells)
  values <- rbind(cp_step(n, 201L, 0, 5), cp_step(n, 201L, 5, 0),
    cp_step(n, 201L, 1, 4))
  dimnames(values) <- list(c("f1", "f2", "f3"), cells)
  wb_attach_matrix(bundle, "RNA", Matrix::Matrix(values, sparse = TRUE),
    feature_type = "gene")
  base <- list(matrix = "RNA", chunk_size = 4L,
    pseudotime = wb_prim_cell_vector(setNames((seq_len(n) - 1L) / (n - 1L),
      cells), "continuous"))
  run <- function(extra = list())
    wb_session_run_op(wb_session_new(), bundle,
      "scwbtrajectory::trajectory_single_changepoint_regions",
      inputs = c(base, extra))$run$outputs
  before <- wb_art_list(bundle)
  strict <- try(run(), silent = TRUE)
  stopifnot(inherits(strict, "try-error"), grepl(
    "produced no interior breakpoint",
    conditionMessage(attr(strict, "condition")), fixed = TRUE),
    # The failed promotion wrote nothing to the bundle.
    identical(before, wb_art_list(bundle)))
  features <- run(list(degenerate_bkp_policy = "drop"))
  published <- features$feature_single_changepoints$table
  summary <- features$run_summary$table
  regions <- features$trend_regions$table
  dropped <- !published$analysed
  kept <- published$analysed
  stopifnot(identical(summary$n_cells_analysed, n),
    identical(summary$degenerate_bkp_policy, "drop"),
    identical(published$feature_id, c("f1", "f2", "f3")),
    # Executed fact: under this geometry two of the three amplitudes report no
    # interior breakpoint and the third reports a spurious one (rank 270 for a
    # step at rank 200). Which amplitude lands on which outcome is floating-point
    # luck, so the contract asserted here is what happens to a feature that finds
    # nothing - never a fabricated boundary or direction, which is exactly what
    # the notebook does - and that a surviving boundary is interior.
    any(dropped),
    identical(published$reason[dropped], rep("degenerate", sum(dropped))),
    all(is.na(published$boundary_rank[dropped])),
    all(is.na(published$direction[dropped])),
    all(is.na(published$delta[dropped])),
    all(is.na(published$density_significant[dropped])),
    all(published$boundary_rank[kept] > 0L &
      published$boundary_rank[kept] < n),
    identical(summary$n_features_selected, 3L),
    identical(summary$n_features_analysed, as.integer(sum(kept))),
    identical(summary$n_features_dropped, as.integer(sum(dropped))),
    # Each surviving breakpoint is spread over the 100-wide density smoothing and
    # scaled by 20, so the peak is 20 * count / 100 - far below the threshold of
    # 30. No rank is significant and every region row stays empty.
    identical(summary$significant_rank_count, 0L),
    cp_near(max(features$trend_density$value), 20 * sum(kept) / 100),
    length(features$trend_density$value) == 3L * n,
    identical(regions$feature_ids_json, rep("[]", 2L)),
    identical(regions$feature_count, c(0L, 0L)))
  TRUE
}

# ---------------------------------------------------------------------------
# Correction 1: a region is gated by the density of its own direction.
#
# The fixture plants three UP steps before c16, one DOWN step before c16 and two
# DOWN steps before c26. On the 30-cell axis their boundaries are rank 15 and rank
# 25, and with density_smooth_window 1 and density_scale 1 the published density is
# exactly the breakpoint count of that direction at that rank, so every number
# below is arithmetic on the planted steps rather than a value read from a run:
#
#   rank 15: up 3, down 1, combined 4      rank 25: up 0, down 2, combined 2
#
# Asking each direction about its own density therefore gives, at change_threshold
# 3, the three up features and no down feature; and at change_threshold 2, dn2 and
# dn3 (their own density reaches 2 at rank 25) but still not dn1, whose own density
# is 1 even though rank 15 carries the highest combined density in the run. Both
# regions here are whole-axis, so no rank window can be blamed for an exclusion. A
# rule gated on the combined density admits dn1 at both thresholds.
# ---------------------------------------------------------------------------
test_trajectory_single_changepoint_direction_specific_threshold_contract <- function() {
  n <- cp_test$n
  x <- cp_matrix_fixture(list(up1 = cp_step(n, 16L, 0, 9),
    up2 = cp_step(n, 16L, 1, 7), up3 = cp_step(n, 16L, 2, 8),
    dn1 = cp_step(n, 16L, 9, 0), dn2 = cp_step(n, 26L, 9, 1),
    dn3 = cp_step(n, 26L, 8, 0)))
  on.exit(trajectory_close(x), add = TRUE)
  strict <- cp_run(x, list(change_threshold = 3))
  features <- strict$feature_single_changepoints
  regions <- strict$trend_regions
  stopifnot(
    identical(features$feature_id, x$features),
    identical(features$boundary_rank, c(15L, 15L, 15L, 15L, 25L, 25L)),
    identical(features$direction, c("up", "up", "up", "down", "down", "down")),
    # The density arithmetic itself.
    identical(cp_at(strict$trend_density, "up", 15L), 3),
    identical(cp_at(strict$trend_density, "down", 15L), 1),
    identical(cp_at(strict$trend_density, "all", 15L), 4),
    identical(cp_at(strict$trend_density, "up", 25L), 0),
    identical(cp_at(strict$trend_density, "down", 25L), 2),
    identical(cp_at(strict$trend_density, "all", 25L), 2),
    # up reaches 3 at rank 15, no direction ever reaches 3 at rank 25.
    identical(strict$run_summary$significant_rank_count_up, 1L),
    identical(strict$run_summary$significant_rank_count_down, 0L),
    identical(strict$run_summary$significant_rank_count, 1L),
    identical(cp_members(regions, "up"), c("up1", "up2", "up3")),
    identical(cp_members(regions, "down"), character(0)),
    identical(regions$significance_basis, c("up_density", "down_density")),
    identical(strict$run_summary$region_significance_basis,
      "direction_matched_density_v1"),
    identical(strict$run_summary$combined_density_role, "diagnostic_only"),
    # Per feature, both readings side by side: dn1 is significant on the combined
    # density and not on its own, and that is the difference between sitting in a
    # region and not sitting in one.
    identical(features$density_significant,
      c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)),
    identical(features$direction_density_significant,
      c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)))
  loose <- cp_run(x, list(change_threshold = 2))
  loose_features <- loose$feature_single_changepoints
  stopifnot(
    identical(cp_members(loose$trend_regions, "up"), c("up1", "up2", "up3")),
    identical(cp_members(loose$trend_regions, "down"), c("dn2", "dn3")),
    identical(loose$run_summary$significant_rank_count_up, 1L),
    identical(loose$run_summary$significant_rank_count_down, 1L),
    identical(loose$run_summary$significant_rank_count, 2L),
    identical(loose_features$density_significant, rep(TRUE, 6L)),
    identical(loose_features$direction_density_significant,
      c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)))
  TRUE
}

# ---------------------------------------------------------------------------
# Correction 2: missing matrix values follow non_finite_policy.
#
# patchy carries NA in exactly three cells (c03, c09, c21) and a step before c21;
# clean steps before c16 with no gap; flat is flat. The expectations are the
# planted facts: a missing-value feature is published unanalysed with its own
# count of three and nothing invented for it, the feature beside it is still found
# at exactly the rank its step was planted at, and under "error" the run is refused
# with a message that names the feature and counts its missing values instead of
# quoting a Python conversion failure.
# ---------------------------------------------------------------------------
test_trajectory_single_changepoint_missing_value_policy_contract <- function() {
  n <- cp_test$n
  patchy <- cp_step(n, 21L, 0, 9)
  patchy[c(3L, 9L, 21L)] <- NA_real_
  x <- cp_matrix_fixture(list(clean = cp_step(n, 16L, 0, 9), patchy = patchy,
    flat = rep(4, n)))
  on.exit(trajectory_close(x), add = TRUE)
  before <- wb_art_list(x$bundle)
  refused <- try(cp_run(x, list(non_finite_policy = "error",
    features = c("clean", "patchy"))), silent = TRUE)
  stopifnot(inherits(refused, "try-error"))
  message <- conditionMessage(attr(refused, "condition"))
  stopifnot(grepl("patchy", message, fixed = TRUE),
    grepl("3 missing (null)", message, fixed = TRUE),
    grepl("out of 30 cells", message, fixed = TRUE),
    grepl("non_finite_policy='error'", message, fixed = TRUE),
    !grepl("NoneType", message, fixed = TRUE),
    # A refused promotion wrote nothing into the bundle.
    identical(before, wb_art_list(x$bundle)))
  out <- cp_run(x, list(non_finite_policy = "drop", change_threshold = 1))
  features <- out$feature_single_changepoints
  regions <- out$trend_regions
  row <- function(id) features[features$feature_id == id, , drop = FALSE][1L, ]
  stopifnot(
    identical(features$feature_id, c("clean", "patchy", "flat")),
    identical(features$missing_value_count, c(0L, 3L, 0L)),
    identical(row("patchy")$analysed, FALSE),
    identical(row("patchy")$reason, "non_finite"),
    is.na(row("patchy")$boundary_rank),
    is.na(row("patchy")$delta),
    is.na(row("patchy")$direction),
    is.na(row("patchy")$density_significant),
    is.na(row("patchy")$direction_density_significant),
    # The missing values cost only their own feature.
    identical(row("clean")$boundary_rank, 15L),
    identical(row("clean")$missing_value_count, 0L),
    identical(row("clean")$direction, "up"),
    identical(row("flat")$reason, "degenerate"),
    identical(row("flat")$missing_value_count, 0L),
    identical(out$run_summary$n_features_selected, 3L),
    identical(out$run_summary$n_features_analysed, 1L),
    identical(out$run_summary$n_features_dropped, 2L),
    identical(out$run_summary$n_features_with_missing_values, 1L),
    # Only clean contributes a breakpoint, so the up density is 1 at rank 15 and
    # nowhere else; that one rank reaches the threshold of 1 and its sole member is
    # the analysed feature.
    identical(which(cp_key(out$trend_density, "up") >= 1), 16L),
    identical(out$run_summary$significant_rank_count_up, 1L),
    identical(cp_members(regions, "up"), "clean"),
    identical(cp_members(regions, "down"), character(0)),
    identical(regions$feature_count, c(1L, 0L)))
  TRUE
}

# ---------------------------------------------------------------------------
# Correction 3: the boundary-pseudotime summary is about the boundaries.
#
# Steps planted before c16, c21 and c26 on an axis whose cell i has pseudotime
# (i - 1) / 29 put their boundaries at 15/29, 20/29 and 25/29: min 15/29, median
# 20/29 (the middle of three distinct values), max 25/29. flat finds nothing and
# contributes nothing. The cell axis runs the whole way from 0 to 1, so statistics
# taken from the axis cannot land on those three values at once - which is what
# makes this assertion a check rather than a restatement.
# ---------------------------------------------------------------------------
test_trajectory_single_changepoint_boundary_pseudotime_summary_contract <- function() {
  n <- cp_test$n
  x <- cp_matrix_fixture(list(e15 = cp_step(n, 16L, 0, 9),
    e20 = cp_step(n, 21L, 0, 9), e25 = cp_step(n, 26L, 0, 9),
    flat = rep(4, n)))
  on.exit(trajectory_close(x), add = TRUE)
  out <- cp_run(x)
  features <- out$feature_single_changepoints
  summary <- out$run_summary
  published <- stats::na.exclude(features$change_point_pseudotime)
  stopifnot(identical(features$boundary_rank, c(15L, 20L, 25L, NA_integer_)),
    identical(summary$n_features_analysed, 3L),
    cp_near(summary$change_point_pseudotime_min, 15 / 29),
    cp_near(summary$change_point_pseudotime_median, 20 / 29),
    cp_near(summary$change_point_pseudotime_max, 25 / 29),
    # Same quantity as the per-feature column a caller can read; the tolerance is
    # the published double's decimal round trip, not slack in the statistic.
    cp_near(summary$change_point_pseudotime_min, min(published)),
    cp_near(summary$change_point_pseudotime_median, stats::median(published)),
    cp_near(summary$change_point_pseudotime_max, max(published)),
    # Not the cell axis, which is where these columns used to come from.
    !cp_near(summary$change_point_pseudotime_min, 0),
    !cp_near(summary$change_point_pseudotime_max, 1),
    !cp_near(summary$change_point_pseudotime_median,
      stats::median((seq_len(n) - 1L) / (n - 1L))))
  # With no detected boundary at all the three statistics are NA: an unanalysed
  # axis must not borrow the cell axis's own spread.
  none <- cp_run(x, list(features = "flat"))$run_summary
  stopifnot(identical(none$n_features_analysed, 0L),
    is.na(none$change_point_pseudotime_min),
    is.na(none$change_point_pseudotime_median),
    is.na(none$change_point_pseudotime_max))
  TRUE
}

# ---------------------------------------------------------------------------
# Correction 4: region member ids are encoded by a JSON writer.
#
# feature_ids_json is one JSON array inside a character column, so an id holding a
# double quote, a backslash, a newline or a tab has to be escaped or the array
# stops being parseable. Both the encoder and the published column are checked: the
# encoder on hand-built ids, then the real operation on a matrix whose features
# carry those ids, because the ids a caller can plant are the ones that matter.
# ---------------------------------------------------------------------------
test_trajectory_single_changepoint_region_ids_json_contract <- function() {
  nasty <- c('quote"id', "back\\slash", "new\nline", "tab\there")
  encode <- scwbTrajectory:::.trajectory_cp_ids_json
  text <- encode(nasty)
  parsed <- jsonlite::fromJSON(text)
  stopifnot(identical(length(parsed), 4L), all(parsed == nasty),
    # One line: a raw newline inside a JSON string is not valid JSON.
    identical(nchar(text), nchar(gsub("\n", "", text, fixed = TRUE))),
    # Still an array for one member, and [] for none - auto_unbox off.
    identical(encode("solo"), "[\"solo\"]"),
    identical(encode(character(0)), "[]"))
  x <- cp_matrix_fixture(setNames(lapply(seq_along(nasty), function(i)
    cp_step(cp_test$n, 16L, i, i + 9)), nasty))
  on.exit(trajectory_close(x), add = TRUE)
  out <- cp_run(x, list(change_threshold = 1))
  features <- out$feature_single_changepoints
  regions <- out$trend_regions
  published <- jsonlite::fromJSON(regions$feature_ids_json[1L])
  stopifnot(identical(features$feature_id, nasty),
    identical(regions$direction, c("up", "down")),
    identical(length(published), 4L), all(published == nasty),
    identical(regions$feature_count[1L], length(published)),
    identical(regions$feature_count, c(4L, 0L)),
    identical(length(cp_members(regions, "up")), 4L),
    all(cp_members(regions, "up") == nasty),
    # The published column itself is a single line of JSON for four awkward ids.
    identical(nchar(regions$feature_ids_json[1L]),
      nchar(gsub("\n", "", regions$feature_ids_json[1L], fixed = TRUE))))
  TRUE
}
