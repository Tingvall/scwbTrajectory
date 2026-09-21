test_trajectory_promotion_reopen_and_terminal_policy_contract <- function() {
  x <- trajectory_fixture(); b <- x$bundle
  on.exit(trajectory_close(x), add = TRUE)
  session <- wb_session_new()
  branch <- wb_session_run_op(session, b,
    "scwbtrajectory::trajectory_branch_compare", inputs = list(
      matrix = "RNA", pseudotime = x$pseudotime_family,
      branch_membership = x$branches, branch_a = "A", branch_b = "B",
      features = x$features, min_cells = 1L))
  promoted <- wb_session_promote_current(b, branch$session,
    outputs = names(branch$run$outputs),
    created_by = "trajectory-contract")
  stopifnot(setequal(promoted$promoted$output,
      c("branch_de_table", "branch_summary")),
    all(promoted$promoted$status == "ok"))
  details <- wb_prov_run_details(b, branch$run$run_id,
    promoted$provenance_snapshot_id)
  stopifnot(all(details$outputs$replay_equivalence_policy ==
      "exact_canonical_value"),
    all(details$outputs$verification_status == "content_verified"))
  artifacts <- promoted$promoted$artifact_id
  stopifnot(all(vapply(artifacts, function(id)
    inherits(wb_art_load(b, id), "wb_primitive"), logical(1))))
  stopifnot(isTRUE(wb_finalize_recipient_bundle(b)$ready))
  wb_close_bundle(b)
  registry <- getFromNamespace(".wb_operation_pack_registry", "scworkbench")
  saved_pack <- registry$packs$scwbtrajectory
  saved_operations <- registry$operations[saved_pack$operation_names]
  on.exit({registry$packs$scwbtrajectory <- saved_pack
    registry$operations[names(saved_operations)] <- saved_operations},
    add = TRUE)
  for (name in saved_pack$operation_names) registry$operations[[name]] <- NULL
  registry$packs$scwbtrajectory <- NULL
  stopifnot(!any(startsWith(wb_list_ops(), "scwbtrajectory::")))
  b <- wb_open_bundle(x$root, read_only = TRUE)
  x$bundle <- b
  stopifnot(all(vapply(artifacts, function(id)
    inherits(wb_art_load(b, id), "wb_primitive"), logical(1))))
  wb_close_bundle(b)
  x$bundle <- NULL
  stopifnot(isTRUE(wb_bundle_validation_report(x$root, "recipient")$ok))
  TRUE
}

test_trajectory_all_operations_promotion_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  operations <- c("trajectory_branch_compare", "trajectory_route_smooth",
    "space_path_score", "space_path_project_2d")
  for (operation in operations) {
    result <- trajectory_execute_one(x, operation)
    promoted <- wb_session_promote_current(x$bundle, result$session,
      outputs = names(result$run$outputs),
      created_by = "trajectory-all-operation-contract")
    details <- wb_prov_run_details(x$bundle, result$run$run_id,
      promoted$provenance_snapshot_id)
    stopifnot(all(promoted$promoted$status == "ok"),
      setequal(promoted$promoted$output, names(result$run$outputs)),
      all(details$outputs$verification_status == "content_verified"),
      all(details$outputs$replay_equivalence_policy ==
        "exact_canonical_value"))
  }
  TRUE
}

test_trajectory_publication_fault_rollback_and_retry_contract <- function() {
  x <- trajectory_fixture(); b <- x$bundle
  on.exit({options(scworkbench.promotion_workspace_publish_hook = NULL)
    trajectory_close(x)}, add = TRUE)
  operations <- c("trajectory_branch_compare", "trajectory_route_smooth",
    "space_path_score", "space_path_project_2d")
  for (operation in operations) {
    result <- trajectory_execute_one(x, operation)
    outputs <- names(result$run$outputs)
    before <- wb_art_list(b)
    before_state <- utils::getFromNamespace(".wb_session_source_identity",
      "scworkbench")(b)$scientific_state_fingerprint
    options(scworkbench.promotion_workspace_publish_hook = function(...)
      stop("injected trajectory publication failure"))
    failure <- try(wb_session_promote_current(b, result$session,
      outputs = outputs), silent = TRUE)
    options(scworkbench.promotion_workspace_publish_hook = NULL)
    stopifnot(inherits(failure, "try-error"),
      identical(before, wb_art_list(b)),
      identical(before_state, utils::getFromNamespace(
        ".wb_session_source_identity",
        "scworkbench")(b)$scientific_state_fingerprint))
    retry <- wb_session_promote_current(b, result$session,
      outputs = outputs)
    stopifnot(setequal(retry$promoted$output, outputs),
      all(retry$promoted$status == "ok"))
  }
  TRUE
}

test_trajectory_pack_ownership_contract <- function() {
  trajectory_activate_pack()
  qualified <- paste0("scwbtrajectory::", c("trajectory_branch_compare",
    "trajectory_route_smooth", "space_path_score",
    "space_path_project_2d", "space_path_define",
    "trajectory_segment_from_artifacts", "trajectory_space_path_segment",
    "trajectory_segment_join", "trajectory_segment_path_join"))
  stopifnot(all(qualified %in% wb_list_ops()),
    !any(sub("^scwbtrajectory::", "", qualified) %in% wb_list_ops()),
    all(!vapply(c("wb_op_trajectory_branch_compare",
      "wb_op_trajectory_route_smooth", "wb_op_space_path_score",
      "wb_op_space_path_project_2d", "wb_op_space_path_define",
      "wb_op_trajectory_segment_from_artifacts",
      "wb_op_trajectory_space_path_segment", "wb_op_trajectory_segment_join",
      "wb_op_trajectory_segment_path_join"), exists, logical(1),
      envir = asNamespace("scworkbench"), inherits = FALSE)))
  TRUE
}
