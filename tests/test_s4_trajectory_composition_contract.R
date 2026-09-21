.trajectory_s4_method <- function(name)
  getFromNamespace(name, "scwbTrajectory")

.trajectory_s4_pt <- function(cells, values, fit_id, key = "L1")
  wb_prim_cell_vector_family(rep(key, length(cells)), cells, values,
    value_type = "continuous", meta = list(trajectory_fit_id = fit_id,
      engine = "slingshot", family_kind = "trajectory_pseudotime"))

.trajectory_s4_path <- function(fit_id, key = "L1")
  wb_prim_path_2d_family(rep(key, 3L), 1:3,
    seq(0, 1, length.out = 3L), seq(0, 1, length.out = 3L),
    meta = list(trajectory_fit_id = fit_id, engine = "slingshot"))

test_trajectory_s4_old_core_payload_equivalence_contract <- function() {
  trajectory_activate_pack()
  root <- tempfile("scwb-trajectory-s4-equivalence-")
  b <- wb_init_bundle(root, overwrite = TRUE)
  on.exit({try(wb_close_bundle(b), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
  cells <- paste0("c", 1:9); wb_define_cells(b, cells)
  fitted_method <- .trajectory_s4_method(
    ".trajectory_segment_from_artifacts")
  bridge_method <- .trajectory_s4_method(".trajectory_space_path_segment")
  join_method <- .trajectory_s4_method(".trajectory_segment_join")
  path_join_method <- .trajectory_s4_method(".trajectory_segment_path_join")
  fitted_1 <- fitted_method(pseudotime = .trajectory_s4_pt(
      c("c1", "c2", "c3", "c4"), 0:3, "fit_s1"), lineage = "L1",
    segment_id = "fitted_1", segment_order = 1L,
    path = .trajectory_s4_path("fit_s1"), source_engine = "slingshot",
    start_label = "STHSC/root", end_label = "DN3")
  stopifnot(identical(trajectory_identity(b,
      fitted_1$outputs$pseudotime)$content_sha256,
      "1b8dbc4d798b76d57142a59374005c7dd3a21acfa0daf4729da3871f7d116cdb"),
    identical(trajectory_identity(b,
      fitted_1$outputs$confidence)$content_sha256,
      "5d122dc574449369e4f8df52826132ef79aea563f6bdae850fea3c2dbff61e5f"))
  scores <- wb_prim_cell_vector_family(rep(c("progress", "distance",
      "nearest_segment", "nearest_anchor"), each = 3L),
    rep(c("c4", "c5", "c6"), 4L),
    c(0, .5, 1, .05, .10, .05, 1, 1, 1, 1, 1, 2), "continuous")
  bridge <- bridge_method(path_scores = scores, segment_id = "bridge_dn3",
    segment_order = 2L, distance_threshold = .11, start_label = "DN3",
    end_label = "post_bridge")
  stopifnot(identical(trajectory_identity(b,
      bridge$outputs$pseudotime)$content_sha256,
      "94faf7c8323ba4f6635dbd330b4dec3cf73a6502ae77b50659deefd4060b7918"),
    identical(trajectory_identity(b,
      bridge$outputs$confidence)$content_sha256,
      "5e1aff4fcaa5aaf228f0bf2a8c5ad233b20e9503f25d7be68aab7df22181caa6"))
  fitted_3 <- fitted_method(pseudotime = .trajectory_s4_pt(
      c("c6", "c7", "c8", "c9"), 0:3, "fit_s3"), lineage = "L1",
    segment_id = "fitted_3", segment_order = 3L,
    path = .trajectory_s4_path("fit_s3"), source_engine = "slingshot")
  joined <- join_method(segments = list(
    list(segment_id = "fitted_1", order = 1L, source_kind = "fitted_model",
      source_engine = "slingshot", pseudotime = fitted_1$outputs$pseudotime,
      lineage = "L1"),
    list(segment_id = "bridge_dn3", order = 2L,
      source_kind = "space_path_bridge", source_engine = "space_path",
      pseudotime = bridge$outputs$pseudotime, lineage = "bridge_dn3",
      membership = bridge$outputs$membership,
      confidence = bridge$outputs$confidence,
      confidence_key = "bridge_dn3"),
    list(segment_id = "fitted_3", order = 3L, source_kind = "fitted_model",
      source_engine = "slingshot", pseudotime = fitted_3$outputs$pseudotime,
      lineage = "L1")), composite_id = "mixed_contract",
    composite_key = "mixed_contract")
  stopifnot(identical(trajectory_identity(b,
      joined$outputs$composite_pseudotime_unit_segment)$content_sha256,
      "4d66bf32b9d75d9d3184212651941ccd85a99a280c19fb76ab42865c54b582e0"),
    identical(trajectory_identity(b,
      joined$outputs$composite_pseudotime_cell_weighted_0_1)$content_sha256,
      "a776c36f80e825d4621de2803ec3d37b9e603fc6c3a33400dc91a22fe326dce0"),
    identical(joined$outputs$composite_membership$cell_id,
      sort(unique(cells), method = "radix")))
  paths <- path_join_method(paths = list(
    list(segment_id = "fitted_1", order = 1L, source_kind = "fitted_model",
      source_engine = "slingshot", path = .trajectory_s4_path("fit_s1"),
      lineage = "L1"),
    list(segment_id = "fitted_3", order = 2L, source_kind = "fitted_model",
      source_engine = "slingshot", path = .trajectory_s4_path("fit_s3"),
      lineage = "L1")), composite_id = "display_only",
    composite_key = "display_only")
  stopifnot(identical(paths$outputs$curves$kind, "path_2d_family"),
    is.null(paths$outputs$curves$meta$trajectory_fit_id),
    grepl("not a merged model", paths$outputs$curves$meta$caveats,
      fixed = TRUE))
  TRUE
}

test_trajectory_s4_g7_lifecycle_and_rollback_contract <- function() {
  trajectory_activate_pack()
  root <- tempfile("scwb-trajectory-s4-lifecycle-")
  b <- wb_init_bundle(root, overwrite = TRUE)
  on.exit({options(scworkbench.operation_pack_v2_test_hook = NULL)
    try(wb_close_bundle(b), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
  cells <- c("c2", "c10", "c1", "c3", "c4")
  wb_define_cells(b, cells)
  space <- cbind(d1 = seq_along(cells), d2 = seq_along(cells)^2)
  rownames(space) <- cells; wb_add_space(b, "representation", space)
  wb_cell_set_replace(b, "set:start", cells[1:2], source = "s4-test")
  wb_cell_set_replace(b, "set:end", cells[4:5], source = "s4-test")
  path <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::space_path_define", inputs = list(
      space = "representation", checkpoint_sets = c("set:start", "set:end"),
      n_path_points = 5L, path_id_out = "representation_path",
      min_cells = 1L))
  score <- wb_session_run_op(path$session, b,
    "scwbtrajectory::space_path_score", inputs = list(
      space = "representation", path = wb_ref_run(path$run$run_id, "path"),
      min_cells = 1L))
  score_promotion <- wb_session_promote_current(b, score$session,
    outputs = c("path_scores", "fit_summary"), created_by = "s4-g7")
  score_id <- score_promotion$promoted$artifact_id[
    score_promotion$promoted$output == "path_scores"][[1L]]
  bridge <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::trajectory_space_path_segment", inputs = list(
      path_scores = wb_ref_artifact(score_id, "path_scores"),
      segment_id = "bridge", segment_order = 2L,
      distance_quantile = 1, pseudotime_key = "custom_pt",
      confidence_key = "custom_conf", min_cells = 1L))
  receipt_path <- file.path(getFromNamespace(".wb_operation_pack_v2_store",
      "scworkbench")(bridge$session$staging_bundle),
    bridge$run$operation_pack_transition$transition_id, "receipt.json")
  receipt <- jsonlite::read_json(receipt_path, simplifyVector = FALSE)
  stopifnot(length(receipt$derived_cell_output_axis_observations) == 1L,
    identical(bridge$run$outputs$pseudotime$key,
      rep("custom_pt", length(cells))),
    identical(bridge$run$outputs$confidence$key,
      rep("custom_conf", length(cells))))
  promoted <- wb_session_promote_current(b, bridge$session,
    outputs = names(bridge$run$outputs), created_by = "s4-g7")
  stopifnot(all(promoted$promoted$status == "ok"))
  artifact_ids <- promoted$promoted$artifact_id
  before <- wb_art_list(b)
  before_state <- getFromNamespace(".wb_session_source_identity",
    "scworkbench")(b)$scientific_state_fingerprint
  options(scworkbench.operation_pack_v2_test_hook = function(event, ...)
    if (identical(event, "before_publication")) stop("s4 fault"))
  failure <- try(wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::trajectory_segment_from_artifacts", inputs = list(
      pseudotime = .trajectory_s4_pt(cells, seq(0, 1,
        length.out = length(cells)), "fit_fault"), lineage = "L1",
      segment_id = "fault", min_cells = 1L)), silent = TRUE)
  options(scworkbench.operation_pack_v2_test_hook = NULL)
  stopifnot(inherits(failure, "try-error"), identical(before, wb_art_list(b)),
    identical(before_state, getFromNamespace(".wb_session_source_identity",
      "scworkbench")(b)$scientific_state_fingerprint))
  stopifnot(isTRUE(wb_finalize_recipient_bundle(b)$ready))
  wb_close_bundle(b)
  b <- wb_open_bundle(root, read_only = TRUE)
  stopifnot(all(vapply(artifact_ids, function(id)
    inherits(wb_art_load(b, id), "wb_primitive"), logical(1))))
  wb_close_bundle(b); b <- NULL
  TRUE
}
