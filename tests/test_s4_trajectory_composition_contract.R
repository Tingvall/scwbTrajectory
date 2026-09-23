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

.trajectory_s4_sha <- function(bundle, value)
  trajectory_identity(bundle, value)$content_sha256

.trajectory_s4_assert_family <- function(value, cells, key, expected,
    value_semantics, tolerance = 0) {
  stopifnot(identical(value$kind, "cell_vector_family"),
    identical(value$cell_id, cells),
    identical(value$key, rep(key, length(cells))),
    identical(value$value_type, "continuous"),
    isTRUE(all.equal(as.numeric(value$value), as.numeric(expected),
      tolerance = tolerance)),
    identical(sum(names(value$meta) == "value_semantics"), 1L),
    identical(value$meta$value_semantics, value_semantics))
  invisible(TRUE)
}

.trajectory_s4_identity_ledger <- function() {
  path <- system.file("scworkbench",
    "s4-metadata-identity-change-ledger.csv", package = "scwbTrajectory")
  stopifnot(nzchar(path), file.exists(path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

test_trajectory_s4_scientific_equivalence_with_metadata_identity_exception_contract <- function() {
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
  stopifnot(identical(names(fitted_1$outputs), c("pseudotime",
      "membership", "confidence", "segment_manifest")),
    identical(fitted_1$outputs$membership$cell_id,
      c("c1", "c2", "c3", "c4")),
    identical(fitted_1$outputs$membership$value,
      rep("fitted_1", 4L)),
    identical(fitted_1$outputs$membership$value_type, "categorical"),
    identical(.trajectory_s4_sha(b, fitted_1$outputs$pseudotime),
      "1b8dbc4d798b76d57142a59374005c7dd3a21acfa0daf4729da3871f7d116cdb"),
    identical(.trajectory_s4_sha(b, fitted_1$outputs$membership),
      "8bd936db1e47a3b3b594c4c30a5f93d5b300b0f8b2fa1bf17d960eb270a568c5"),
    identical(.trajectory_s4_sha(b, fitted_1$outputs$segment_manifest),
      "8eded9390d2ebf08cdae7ec72a6a78a07f0e48e491fb293df4b79ba2dc899326"),
    identical(.trajectory_s4_sha(b, fitted_1$outputs$confidence),
      "131b91351453b4a8ca7d76733815906391025be0c78c3ca23f1a72111e95acd1"),
    identical(fitted_1$params$segment_id, "fitted_1"),
    identical(fitted_1$params$segment_order, 1),
    identical(fitted_1$params$pseudotime_key, "L1"),
    identical(fitted_1$params$confidence_key, "L1"),
    identical(fitted_1$params$trajectory_fit_id, "fit_s1"))
  .trajectory_s4_assert_family(fitted_1$outputs$pseudotime,
    c("c1", "c2", "c3", "c4"), "L1", 0:3,
    "trajectory_segment_fitted_pseudotime")
  .trajectory_s4_assert_family(fitted_1$outputs$confidence,
    c("c1", "c2", "c3", "c4"), "L1", rep(1, 4L),
    "segment_confidence_or_membership_weight")
  scores <- wb_prim_cell_vector_family(rep(c("progress", "distance",
      "nearest_segment", "nearest_anchor"), each = 3L),
    rep(c("c4", "c5", "c6"), 4L),
    c(0, .5, 1, .05, .10, .05, 1, 1, 1, 1, 1, 2), "continuous")
  bridge <- bridge_method(path_scores = scores, segment_id = "bridge_dn3",
    segment_order = 2L, distance_threshold = .11, start_label = "DN3",
    end_label = "post_bridge")
  stopifnot(identical(names(bridge$outputs), c("pseudotime",
      "membership", "confidence", "segment_manifest",
      "display_endpoint_qc")),
    identical(.trajectory_s4_sha(b, bridge$outputs$pseudotime),
      "94faf7c8323ba4f6635dbd330b4dec3cf73a6502ae77b50659deefd4060b7918"),
    identical(.trajectory_s4_sha(b, bridge$outputs$membership),
      "71d4730c82138ad8f248a855c081e76e0948b8ee1c78ddf49c261e56fd956697"),
    identical(.trajectory_s4_sha(b, bridge$outputs$segment_manifest),
      "f0fa4078e1ebf381372906315053d49151b858d626b1a34c084856a6b92801ae"),
    identical(.trajectory_s4_sha(b, bridge$outputs$confidence),
      "ad339b4be4b727c393512cd7e61dda75a88a1595cffe624febfa3bcfcda92bcc"),
    identical(bridge$params$segment_id, "bridge_dn3"),
    identical(bridge$params$segment_order, 2),
    identical(bridge$params$pseudotime_key, "bridge_dn3"),
    identical(bridge$params$confidence_key, "bridge_dn3"))
  .trajectory_s4_assert_family(bridge$outputs$pseudotime,
    c("c4", "c5", "c6"), "bridge_dn3", c(0, .5, 1),
    "trajectory_segment_path_progress")
  .trajectory_s4_assert_family(bridge$outputs$confidence,
    c("c4", "c5", "c6"), "bridge_dn3",
    1 - c(.05, .10, .05) / .11,
    "distance_derived_segment_confidence", tolerance = 1e-14)
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
  joined_cells <- sort(unique(cells), method = "radix")
  stopifnot(identical(names(joined$outputs), c(
      "composite_pseudotime_unit_segment",
      "composite_pseudotime_cell_weighted_0_1", "segment_local_table",
      "composite_manifest", "composite_membership")),
    identical(joined$outputs$composite_membership$cell_id, joined_cells),
    identical(joined$outputs$composite_membership$value,
      c("fitted_1", "fitted_1", "fitted_1", "waypoint",
        "bridge_dn3", "waypoint", "fitted_3", "fitted_3",
        "fitted_3")),
    identical(.trajectory_s4_sha(b, joined$outputs$segment_local_table),
      "b099313ee7eb882bc21af5a1b87a2ca3476df70903fcd4a81330f316f812137d"),
    identical(.trajectory_s4_sha(b, joined$outputs$composite_manifest),
      "f2f1148686bc3b6c18f55e3c08e90955bbde828e834b93c0affd12660f10785e"),
    identical(.trajectory_s4_sha(b, joined$outputs$composite_membership),
      "63174309eb1226c60755c0705916a46720f9cbb0c6195fcc1fa25ca2bc0caf27"),
    identical(.trajectory_s4_sha(b,
      joined$outputs$composite_pseudotime_unit_segment),
      "e56ed73cca558739fe825a92a9512dec64e35c2097919d718ed93886c54dc3d4"),
    identical(.trajectory_s4_sha(b,
      joined$outputs$composite_pseudotime_cell_weighted_0_1),
      "6c4e433bab08ce0494884a45737b5c423075bb76c3c1a23f215f828d59e37a7b"),
    identical(joined$params$trajectory_composite_id, "mixed_contract"),
    identical(joined$params$composite_key, "mixed_contract"),
    identical(joined$params$join_policy, "exact_boundary_collapse"),
    identical(joined$params$parent_segment_ids,
      c("fitted_1", "bridge_dn3", "fitted_3")),
    identical(joined$params$parent_trajectory_fit_ids,
      c("fit_s1", "fit_s3")))
  .trajectory_s4_assert_family(
    joined$outputs$composite_pseudotime_unit_segment, joined_cells,
    "mixed_contract", c(0, 1 / 3, 2 / 3, 1, 1.5, 2, 7 / 3,
      8 / 3, 3), "composite_pseudotime_unit_segment", tolerance = 1e-14)
  .trajectory_s4_assert_family(
    joined$outputs$composite_pseudotime_cell_weighted_0_1, joined_cells,
    "mixed_contract", c(0, 4 / 33, 8 / 33, 12 / 33, .5,
      21 / 33, 25 / 33, 29 / 33, 1),
    "composite_pseudotime_cell_weighted_0_1", tolerance = 1e-14)
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

test_trajectory_s4_metadata_identity_change_ledger_contract <- function() {
  trajectory_activate_pack()
  ledger <- .trajectory_s4_identity_ledger()
  stopifnot(nrow(ledger) == 4L,
    identical(ledger$change_class, rep("metadata_only", 4L)),
    all(ledger$old_content_sha256 != ledger$new_content_sha256),
    identical(unique(ledger$old_pack_version), "0.1.0"),
    identical(unique(ledger$new_pack_version), "0.1.1"),
    all(grepl("Duplicate-named value_semantics metadata was normalized",
      ledger$reason, fixed = TRUE)))
  stopifnot(identical(ledger$old_content_sha256, c(
      "5d122dc574449369e4f8df52826132ef79aea563f6bdae850fea3c2dbff61e5f",
      "5e1aff4fcaa5aaf228f0bf2a8c5ad233b20e9503f25d7be68aab7df22181caa6",
      "4d66bf32b9d75d9d3184212651941ccd85a99a280c19fb76ab42865c54b582e0",
      "a776c36f80e825d4621de2803ec3d37b9e603fc6c3a33400dc91a22fe326dce0")),
    identical(ledger$new_content_sha256, c(
      "131b91351453b4a8ca7d76733815906391025be0c78c3ca23f1a72111e95acd1",
      "ad339b4be4b727c393512cd7e61dda75a88a1595cffe624febfa3bcfcda92bcc",
      "e56ed73cca558739fe825a92a9512dec64e35c2097919d718ed93886c54dc3d4",
      "6c4e433bab08ce0494884a45737b5c423075bb76c3c1a23f215f828d59e37a7b")))
  manifest_path <- system.file("scworkbench", "operation-pack.json",
    package = "scwbTrajectory")
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  identity <- getFromNamespace(".wb_operation_pack_v3_id",
    "scworkbench")
  manifest_identity <- identity(
    "scworkbench.operation_pack.manifest.v3", manifest)$sha256
  stopifnot(identical(manifest$package$version, "0.1.1"),
    identical(manifest_identity,
      "3d72ae2afdc762f759ee1b00dfa5a0206aed7805314685398ec174d4d55da64e"),
    identical(as.character(utils::packageVersion("scwbTrajectory")),
      "0.1.1"))
  expected <- c(
    trajectory_segment_from_artifacts =
      "5b255d889af5efabe7178faefbb22e854b91f54f1b0146e8d91f226471bc5052",
    trajectory_space_path_segment =
      "c70f22e162680e32b0be147e1325b6964580ce3760b4df7056299ae3e37ffb06",
    trajectory_segment_join =
      "331a6b99f09ce0a1560b060a9bd548aaa0ed04cbbb129102b958299ef16a6464",
    trajectory_segment_path_join =
      "17a0c2038db02110bfe7f621660c8ec3a261e30cff35e3345ba453a2b514457f")
  for (name in names(expected)) {
    operation <- manifest$operations[[match(name, vapply(
      manifest$operations, `[[`, character(1), "name"))]]
    stopifnot(identical(identity(
      "scworkbench.operation_pack.definition.v3", operation)$sha256,
      unname(expected[[name]])))
  }
  affected <- manifest$operations[vapply(manifest$operations,
    function(operation) operation$name %in% names(expected)[1:3],
    logical(1))]
  stopifnot(all(vapply(affected, function(operation)
    grepl("canonical_value_semantics_v2", operation$description,
      fixed = TRUE), logical(1))))
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
  receipt_outputs <- stats::setNames(receipt$outputs,
    vapply(receipt$outputs, `[[`, character(1), "name"))
  stopifnot(length(receipt$derived_cell_output_axis_observations) == 1L,
    identical(bridge$run$outputs$pseudotime$key,
      rep("custom_pt", length(cells))),
    identical(bridge$run$outputs$confidence$key,
      rep("custom_conf", length(cells))),
    identical(sum(names(bridge$run$outputs$confidence$meta) ==
      "value_semantics"), 1L),
    identical(bridge$run$outputs$confidence$meta$value_semantics,
      "distance_derived_segment_confidence"),
    identical(receipt_outputs$confidence$descriptor$content_identity$scheme,
      "scworkbench.cell_vector_family.content.v1"))
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

test_trajectory_s4_assembled_composite_contract <- function() {
  x <- trajectory_fixture(); b <- x$bundle
  on.exit({
    options(scworkbench.promotion_workspace_publish_hook = NULL)
    trajectory_close(x)
  }, add = TRUE)
  wb_cell_set_replace(b, "set:start", x$cells[1:4], source = "s4-composite")
  wb_cell_set_replace(b, "set:end", x$cells[21:24], source = "s4-composite")

  promote_all <- function(result) wb_session_promote_current(b,
    result$session, outputs = names(result$run$outputs),
    created_by = "s4-composite")
  artifact_id <- function(promotion, operation, output) {
    rows <- promotion$promoted$op_name == operation &
      promotion$promoted$output == output &
      promotion$promoted$status == "ok"
    stopifnot(sum(rows) == 1L)
    promotion$promoted$artifact_id[rows][[1L]]
  }

  path <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::space_path_define", inputs = list(space = "latent",
      checkpoint_sets = c("set:start", "set:end"), n_path_points = 8L,
      path_id_out = "s4_path", min_cells = 1L))
  path_promotion <- promote_all(path)
  path_id <- artifact_id(path_promotion,
    "scwbtrajectory::space_path_define", "path")
  projected_1 <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::space_path_project_2d", inputs = list(
      path = wb_ref_artifact(path_id, "path"), source_space = "latent",
      plot_space = "plot", approx_points = 4L, min_cells = 1L))
  projection_1_promotion <- promote_all(projected_1)
  display_1_id <- artifact_id(projection_1_promotion,
    "scwbtrajectory::space_path_project_2d", "curves")

  projected_2 <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::space_path_project_2d", inputs = list(
      path = wb_ref_artifact(path_id, "path"), source_space = "latent",
      plot_space = "plot", approx_points = 5L, min_cells = 1L))
  display_2_promotion <- promote_all(projected_2)
  display_2_id <- artifact_id(display_2_promotion,
    "scwbtrajectory::space_path_project_2d", "curves")
  scored <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::space_path_score", inputs = list(space = "latent",
      path = wb_ref_artifact(path_id, "path"), min_cells = 1L))
  score_promotion <- promote_all(scored)
  score_id <- artifact_id(score_promotion,
    "scwbtrajectory::space_path_score", "path_scores")

  fitted_1 <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::trajectory_segment_from_artifacts", inputs = list(
      pseudotime = .trajectory_s4_pt(x$cells[1:8],
        seq(0, 1, length.out = 8L), "fit_s1"), lineage = "L1",
      segment_id = "fitted_1", segment_order = 1L,
      source_engine = "slingshot",
      start_label = "root", end_label = "bridge", min_cells = 1L))
  fitted_1_promotion <- promote_all(fitted_1)
  bridge <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::trajectory_space_path_segment", inputs = list(
      path_scores = wb_ref_artifact(score_id, "path_scores"),
      path = wb_ref_artifact(display_1_id, "curves"),
      display_path = wb_ref_artifact(display_1_id, "curves"),
      segment_id = "bridge", segment_order = 2L,
      candidate_cells = x$cells[8:17], distance_quantile = 1,
      display_space_id = "plot", start_endpoint = "set:start",
      end_endpoint = "set:end", endpoint_qc_action = "none",
      start_label = "bridge", end_label = "terminal",
      min_cells = 1L))
  bridge_promotion <- promote_all(bridge)
  fitted_3 <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::trajectory_segment_from_artifacts", inputs = list(
      pseudotime = .trajectory_s4_pt(x$cells[17:24],
        seq(0, 1, length.out = 8L), "fit_s3"), lineage = "L1",
      segment_id = "fitted_3", segment_order = 3L,
      source_engine = "slingshot",
      start_label = "bridge", end_label = "terminal", min_cells = 1L))
  fitted_3_promotion <- promote_all(fitted_3)

  segment_spec <- function(promotion, id, order, lineage, source_kind,
      source_engine) list(segment_id = id, order = order,
    source_kind = source_kind, source_engine = source_engine,
    lineage = lineage,
    pseudotime = wb_ref_artifact(artifact_id(promotion,
      if (identical(source_kind, "space_path_bridge"))
        "scwbtrajectory::trajectory_space_path_segment" else
        "scwbtrajectory::trajectory_segment_from_artifacts", "pseudotime"),
      "pseudotime"),
    membership = wb_ref_artifact(artifact_id(promotion,
      if (identical(source_kind, "space_path_bridge"))
        "scwbtrajectory::trajectory_space_path_segment" else
        "scwbtrajectory::trajectory_segment_from_artifacts", "membership"),
      "membership"),
    confidence = wb_ref_artifact(artifact_id(promotion,
      if (identical(source_kind, "space_path_bridge"))
        "scwbtrajectory::trajectory_space_path_segment" else
        "scwbtrajectory::trajectory_segment_from_artifacts", "confidence"),
      "confidence"),
    segment_manifest = wb_ref_artifact(artifact_id(promotion,
      if (identical(source_kind, "space_path_bridge"))
        "scwbtrajectory::trajectory_space_path_segment" else
        "scwbtrajectory::trajectory_segment_from_artifacts",
      "segment_manifest"), "segment_manifest"))
  joined <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::trajectory_segment_join", inputs = list(segments = list(
      segment_spec(fitted_1_promotion, "fitted_1", 1L, "L1", "fitted_model",
        "slingshot"),
      segment_spec(bridge_promotion, "bridge", 2L, "bridge", "space_path_bridge",
        "space_path"),
      segment_spec(fitted_3_promotion, "fitted_3", 3L, "L1", "fitted_model",
        "slingshot")), composite_id = "s4_composite",
      composite_key = "s4_composite", min_cells = 1L))
  path_spec <- function(artifact, id, order, points) list(segment_id = id,
    order = order, source_kind = "space_path_bridge",
    source_engine = "space_path", lineage = "s4_path",
    path = wb_ref_artifact(artifact, "curves"),
    display_space_id = "plot", start_endpoint = "set:start",
    end_endpoint = "set:end", endpoint_qc_action = "none",
    start_label = paste0("path_", points, "_start"),
    end_label = paste0("path_", points, "_end"))
  path_joined <- wb_session_run_op(wb_session_new(), b,
    "scwbtrajectory::trajectory_segment_path_join", inputs = list(
      paths = list(path_spec(display_1_id, "display_1", 1L, 4L),
        path_spec(display_2_id, "display_2", 2L, 5L)),
      composite_id = "s4_display", composite_key = "s4_display",
      display_space_id = "plot", endpoint_qc_action = "none",
      min_cells = 1L))

  stopifnot(identical(names(joined$run$outputs), c(
      "composite_pseudotime_unit_segment",
      "composite_pseudotime_cell_weighted_0_1", "segment_local_table",
      "composite_manifest", "composite_membership")),
    identical(joined$run$outputs$composite_membership$cell_id,
      sort(x$cells, method = "radix")),
    identical(as.character(
      joined$run$outputs$composite_manifest$table$segment_id),
      c("fitted_1", "bridge", "fitted_3")),
    identical(unique(as.character(
      joined$run$outputs$composite_manifest$table$trajectory_fit_id)[
        nzchar(as.character(joined$run$outputs$composite_manifest$table$
          trajectory_fit_id))]),
      c("fit_s1", "fit_s3")),
    nrow(bridge$run$outputs$display_endpoint_qc$table) == 2L,
    nrow(path_joined$run$outputs$endpoint_qc$table) == 4L,
    is.null(path_joined$run$outputs$curves$meta$trajectory_fit_id),
    grepl("not a merged fitted model",
      path_joined$run$outputs$curves$meta$caveat, fixed = TRUE))

  transition_root <- function(result) file.path(
    getFromNamespace(".wb_operation_pack_v2_store", "scworkbench")(
      result$session$staging_bundle),
    result$run$operation_pack_transition$transition_id, "receipt.json")
  join_receipt <- jsonlite::read_json(transition_root(joined),
    simplifyVector = FALSE)
  path_receipt <- jsonlite::read_json(transition_root(path_joined),
    simplifyVector = FALSE)
  stopifnot(length(join_receipt$derived_cell_output_axis_observations) == 1L,
    identical(join_receipt$derived_cell_output_axis_observations[[1L]]$
      derivation_mode, "declared_ordered_unique_union_v1"),
    length(path_receipt$derived_cell_output_axis_observations) == 0L)

  before <- wb_art_list(b)
  before_state <- getFromNamespace(".wb_session_source_identity",
    "scworkbench")(b)$scientific_state_fingerprint
  options(scworkbench.promotion_workspace_publish_hook = function(...)
    stop("injected S4 composite publication failure"))
  failed <- try(wb_session_promote_current(b, path_joined$session,
    outputs = names(path_joined$run$outputs), created_by = "s4-composite"),
    silent = TRUE)
  options(scworkbench.promotion_workspace_publish_hook = NULL)
  stopifnot(inherits(failed, "try-error"), identical(before, wb_art_list(b)),
    identical(before_state, getFromNamespace(".wb_session_source_identity",
      "scworkbench")(b)$scientific_state_fingerprint))
  path_join_promotion <- wb_session_promote_current(b, path_joined$session,
    outputs = names(path_joined$run$outputs), created_by = "s4-composite")
  join_promotion <- promote_all(joined)
  stopifnot(all(path_join_promotion$promoted$status == "ok"),
    all(join_promotion$promoted$status == "ok"))
  artifact_ids <- wb_art_list(b)$artifact_id
  stopifnot(isTRUE(wb_finalize_recipient_bundle(b)$ready))

  registry <- getFromNamespace(".wb_operation_pack_registry", "scworkbench")
  saved_pack <- registry$packs$scwbtrajectory
  saved_operations <- registry$operations[saved_pack$operation_names]
  on.exit({registry$packs$scwbtrajectory <- saved_pack
    registry$operations[names(saved_operations)] <- saved_operations},
    add = TRUE)
  wb_close_bundle(b); x$bundle <- NULL
  for (name in saved_pack$operation_names) registry$operations[[name]] <- NULL
  registry$packs$scwbtrajectory <- NULL
  b <- wb_open_bundle(x$root, read_only = TRUE); x$bundle <- b
  recipient_report <- wb_bundle_validation_report(x$root, "recipient")
  stopifnot(all(vapply(artifact_ids, function(id)
    inherits(wb_art_load(b, id), "wb_primitive"), logical(1))),
    isTRUE(recipient_report$ok))
  TRUE
}
