trajectory_s5_oracle <- function() {
  filename <- "trajectory_s5_monocle3_old_core_oracle.json"
  directory <- if (file.exists(file.path("fixtures", filename)))
    "." else "tests"
  jsonlite::read_json(file.path(directory, "fixtures", filename),
    simplifyVector = FALSE)
}

trajectory_s5_fixture <- function() {
  trajectory_activate_pack()
  oracle <- trajectory_s5_oracle()
  cells <- as.character(unlist(oracle$fit$source_cells))
  t <- seq(0.15, 1.2, length.out = 8L)
  coordinates <- rbind(cbind(t, rep(0, 8L)),
    cbind(-t, 0.25 * t), cbind(rep(0, 8L), -t))
  rownames(coordinates) <- cells
  colnames(coordinates) <- c("d1", "d2")
  display <- cbind(x = 1.5 * coordinates[, 1L] +
      0.2 * coordinates[, 2L],
    y = 0.8 * coordinates[, 2L] - 0.1 * coordinates[, 1L])
  rownames(display) <- cells
  root <- tempfile("trajectory-s5-m3-")
  bundle <- scworkbench::wb_init_bundle(root)
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_add_space(bundle, "latent", coordinates)
  scworkbench::wb_add_space(bundle, "plot", display)
  list(root = root, bundle = bundle, cells = cells,
    coordinates = coordinates, display = display, oracle = oracle)
}

trajectory_s5_close <- function(x) {
  try(scworkbench::wb_close_bundle(x$bundle), silent = TRUE)
  unlink(x$root, recursive = TRUE, force = TRUE)
}

trajectory_s5_fit <- function(x, extra = list()) {
  inputs <- utils::modifyList(list(cell_set = "all", cell_ids = x$cells,
    fit_space = "latent", dims = as.integer(1:2),
    root_cells = x$cells[[1L]], k = 5L, random_seed = 13L,
    use_partition = FALSE, close_loop = FALSE), extra,
    keep.null = TRUE)
  scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbtrajectory::monocle3_fit", inputs = inputs)
}

trajectory_s5_state <- function(bundle) list(
  artifacts = scworkbench::wb_art_list(bundle),
  fingerprint = getFromNamespace(".wb_session_source_identity",
    "scworkbench")(bundle)$scientific_state_fingerprint)

trajectory_s5_identity <- function(value) {
  ident <- getFromNamespace(".wb_operation_pack_v2_identity_plain",
    "scworkbench")
  ident(paste0("scworkbench.", value$kind, ".content.v1"),
    unclass(value))$sha256
}

trajectory_s5_assert_recipient <- function(path) {
  report <- scworkbench::wb_bundle_validation_report(path,
    level = "recipient")
  if (!isTRUE(report$ok)) stop(paste(c("Recipient validation failed:",
    capture.output(print(report$diagnostics))), collapse = "\n"),
    call. = FALSE)
  TRUE
}

trajectory_s5_promote_fit <- function(x, fit) {
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("pseudotime", "branch_membership", "model"))
  stopifnot(nrow(promoted$promoted) == 3L,
    all(promoted$promoted$status == "ok"))
  as.character(promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]])
}

trajectory_s5_project <- function(bundle, model_id, extra = list()) {
  inputs <- utils::modifyList(list(
    model = scworkbench::wb_ref_artifact(model_id, "model"),
    plot_space = "plot", approx_points = 200L), extra,
    keep.null = TRUE)
  scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    bundle, "scwbtrajectory::monocle3_project", inputs = inputs)
}

trajectory_s5_assert_fit_oracle <- function(run, oracle) {
  outputs <- run$run$outputs
  stopifnot(identical(names(outputs),
      c("pseudotime", "branch_membership", "model")),
    identical(outputs$pseudotime$kind, "cell_vector"),
    identical(outputs$branch_membership$kind, "relation_family"),
    identical(outputs$model$kind, "model_ref"),
    identical(outputs$model$engine, "monocle3"),
    identical(outputs$model$model_type, "monocle3_cds"),
    identical(trajectory_s5_identity(outputs$pseudotime),
      oracle$fit$terminal_content_identities$pseudotime$direct$sha256),
    identical(trajectory_s5_identity(outputs$branch_membership),
      oracle$fit$content_identity$sha256))
  relation <- outputs$branch_membership
  rows <- oracle$fit$ordered_rows
  stopifnot(identical(as.character(relation$key),
      vapply(rows, `[[`, character(1), "key")),
    identical(as.character(relation$cell_id),
      vapply(rows, `[[`, character(1), "cell_id")),
    identical(sprintf("%a", as.numeric(relation$weight)),
      vapply(rows, `[[`, character(1), "weight_hex")),
    is.null(relation$value), is.null(relation$distance),
    identical(unique(as.character(relation$key)),
      as.character(unlist(oracle$fit$ordered_keys))))
  params <- run$run$params
  expected <- oracle$project$fit_params
  stopifnot(identical(as.integer(params$k), 5L),
    identical(as.integer(params$random_seed), 13L),
    identical(as.integer(params$dims), 1:2),
    identical(as.character(params$root_cells),
      as.character(oracle$fit$source_cells[[1L]])),
    identical(as.character(params$cluster_method),
      expected$cluster_method),
    identical(as.integer(params$num_iter),
      as.integer(expected$num_iter)),
    identical(as.numeric(params$partition_qval),
      as.numeric(expected$partition_qval)),
    identical(isTRUE(params$weight), expected$weight),
    identical(isTRUE(params$use_partition), expected$use_partition),
    identical(isTRUE(params$close_loop), expected$close_loop),
    identical(as.character(outputs$model$root_pr_nodes),
      as.character(unlist(expected$root_pr_nodes))),
    identical(as.character(outputs$model$lineage_keys),
      as.character(unlist(expected$lineage_keys))),
    identical(lapply(outputs$model$lineage_paths, as.character),
      lapply(expected$lineage_paths, function(path)
        as.character(unlist(path)))))
  TRUE
}

trajectory_s5_assert_curves <- function(curves, oracle) {
  expected <- oracle$project$curves
  stopifnot(identical(curves$kind, "path_2d_family"),
    identical(as.character(curves$key),
      vapply(expected, `[[`, character(1), "key")),
    identical(as.integer(curves$point),
      vapply(expected, `[[`, integer(1), "point")))
  absolute_tolerance <- as.numeric(oracle$project$comparison$
    numeric_absolute_tolerance)
  relative_tolerance <- as.numeric(oracle$project$comparison$
    numeric_relative_tolerance)
  for (axis in c("x", "y")) {
    observed <- as.numeric(curves[[axis]])
    reference <- vapply(expected, `[[`, numeric(1), axis)
    stopifnot(length(observed) == length(reference),
      all(is.finite(observed)), all(is.finite(reference)),
      all(abs(observed - reference) <= absolute_tolerance +
        relative_tolerance * abs(reference)))
  }
  TRUE
}

test_trajectory_s5_monocle3_curve_comparator_policy_contract <- function() {
  oracle <- trajectory_s5_oracle()
  rows <- oracle$project$curves
  curves <- scworkbench::wb_prim_path_2d_family(
    key = vapply(rows, `[[`, character(1), "key"),
    point = vapply(rows, `[[`, integer(1), "point"),
    x = vapply(rows, `[[`, numeric(1), "x"),
    y = vapply(rows, `[[`, numeric(1), "y"))
  stopifnot(trajectory_s5_assert_curves(curves, oracle))
  altered <- curves
  altered$x[[1L]] <- altered$x[[1L]] + 1e-4
  stopifnot(inherits(try(trajectory_s5_assert_curves(altered,
    oracle), silent = TRUE), "try-error"))
  TRUE
}

test_trajectory_s5_monocle3_registered_oracle_reopen_contract <- function() {
  x <- trajectory_s5_fixture()
  on.exit(trajectory_s5_close(x), add = TRUE)
  fit <- trajectory_s5_fit(x)
  trajectory_s5_assert_fit_oracle(fit, x$oracle)
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("pseudotime", "branch_membership", "model"))
  stopifnot(nrow(promoted$promoted) == 3L,
    all(promoted$promoted$status == "ok"))
  model_id <- as.character(promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]])
  relation_id <- as.character(promoted$promoted$artifact_id[
    promoted$promoted$output == "branch_membership"][[1L]])
  scworkbench::wb_session_discard(fit$session)
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    x$bundle)$ready))
  scworkbench::wb_close_bundle(x$bundle)
  x$bundle <- NULL
  relocated <- tempfile("trajectory-s5-m3-relocated-")
  fs::dir_copy(x$root, relocated)
  on.exit(unlink(relocated, recursive = TRUE, force = TRUE), add = TRUE)
  trajectory_s5_assert_recipient(relocated)
  reopened <- scworkbench::wb_open_bundle(relocated)
  on.exit(if (!is.null(reopened))
    try(scworkbench::wb_close_bundle(reopened), silent = TRUE),
    add = TRUE)
  model <- scworkbench::wb_art_load(reopened, model_id,
    expected_kind = "model_ref")
  stopifnot(identical(model$engine, "monocle3"),
    identical(model$model_type, "monocle3_cds"),
    identical(model$model_source_linkage$source_space_id, "latent"),
    file.exists(file.path(model$model_tree_root,
      model$model_primary_path)))
  observed_tree <- getFromNamespace(
    ".wb_operation_pack_v3_model_tree_manifest", "scworkbench")(
      model$model_tree_root, model$engine, model$model_type,
      model$model_primary_path)
  stopifnot(identical(observed_tree$manifest_identity,
    model$model_tree_manifest$manifest_identity))
  relation <- scworkbench::wb_art_load(reopened, relation_id,
    expected_kind = "relation_family")
  stopifnot(identical(trajectory_s5_identity(relation),
    x$oracle$fit$reopened_content_identity$sha256))
  project <- scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), reopened,
    "scwbtrajectory::monocle3_project", inputs = list(
      model = scworkbench::wb_ref_artifact(model_id, "model"),
      plot_space = "plot", approx_points = 200L, min_cells = 1L))
  stopifnot(identical(names(project$run$outputs), "curves"),
    identical(as.integer(project$run$params$approx_points), 200L),
    identical(scworkbench::wb_path_2d_family_space_id(reopened,
      scworkbench::wb_ref_run(project$run$run_id, "curves"),
      runs = project$session$runs), "plot"))
  trajectory_s5_assert_curves(project$run$outputs$curves, x$oracle)
  curves_promotion <- scworkbench::wb_session_promote_current(reopened,
    project$session, outputs = "curves")
  stopifnot(nrow(curves_promotion$promoted) == 1L,
    identical(curves_promotion$promoted$status[[1L]], "ok"))
  curves_id <- curves_promotion$promoted$artifact_id[[1L]]
  scworkbench::wb_session_discard(project$session)
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    reopened)$ready))
  scworkbench::wb_close_bundle(reopened)
  reopened <- NULL
  trajectory_s5_assert_recipient(relocated)
  recipient <- scworkbench::wb_open_bundle(relocated, read_only = TRUE)
  on.exit(scworkbench::wb_close_bundle(recipient), add = TRUE)
  trajectory_s5_assert_curves(scworkbench::wb_art_load(recipient,
    curves_id, expected_kind = "path_2d_family"), x$oracle)
  TRUE
}

test_trajectory_s5_monocle3_exact_name_ownership_contract <- function() {
  script <- tempfile("trajectory-s5-ownership-", fileext = ".R")
  log <- tempfile("trajectory-s5-ownership-", fileext = ".log")
  on.exit(unlink(c(script, log)), add = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(capture.output(dput(.libPaths())),
      collapse = ""), ")"),
    "library(scworkbench)",
    "names <- c('monocle3_fit', 'monocle3_project')",
    "before <- wb_list_ops()",
    "stopifnot(all(!names %in% before),",
    "  all(!paste0('scwbtrajectory::', names) %in% before))",
    "manifest <- system.file('scworkbench', 'operation-pack.json',",
    "  package = 'scwbTrajectory')",
    "wb_operation_pack_activate(manifest, lib.loc = .libPaths())",
    "after <- wb_list_ops()",
    "stopifnot(all(!names %in% after),",
    "  all(vapply(paste0('scwbtrajectory::', names), function(x)",
    "    sum(after == x) == 1L, logical(1))))",
    "cat('RESULT S5 two qualified names once; no core owner\\n')"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    args = shQuote(script), stdout = log, stderr = log)
  if (!identical(as.integer(status), 0L))
    stop(paste(readLines(log, warn = FALSE), collapse = "\n"),
      call. = FALSE)
  stopifnot(any(readLines(log, warn = FALSE) ==
    "RESULT S5 two qualified names once; no core owner"))
  TRUE
}

test_trajectory_s5_monocle3_invalid_input_zero_residue_contract <- function() {
  x <- trajectory_s5_fixture()
  on.exit(trajectory_s5_close(x), add = TRUE)
  before <- trajectory_s5_state(x$bundle)
  invalid <- list(
    list(root_cells = "foreign_cell"),
    list(dims = c(1L, 3L)),
    list(cell_ids = x$cells[1:9]),
    list(root_pr_nodes = "foreign_vertex"))
  for (extra in invalid) {
    rejected <- try(trajectory_s5_fit(x, extra), silent = TRUE)
    stopifnot(inherits(rejected, "try-error"),
      identical(before, trajectory_s5_state(x$bundle)))
  }
  stopifnot(nrow(scworkbench::wb_art_list(x$bundle)) == 0L)
  TRUE
}

test_trajectory_s5_monocle3_rollback_retry_contract <- function() {
  x <- trajectory_s5_fixture()
  on.exit({
    options(scworkbench.promotion_workspace_publish_hook = NULL)
    trajectory_s5_close(x)
  }, add = TRUE)
  fit <- trajectory_s5_fit(x)
  before <- trajectory_s5_state(x$bundle)
  options(scworkbench.promotion_workspace_publish_hook = function(...)
    stop("injected S5 publication failure"))
  rejected <- try(trajectory_s5_promote_fit(x, fit), silent = TRUE)
  options(scworkbench.promotion_workspace_publish_hook = NULL)
  stopifnot(inherits(rejected, "try-error"),
    grepl("injected S5 publication failure", as.character(rejected),
      fixed = TRUE),
    identical(before, trajectory_s5_state(x$bundle)))
  model_id <- trajectory_s5_promote_fit(x, fit)
  stopifnot(nzchar(model_id),
    nrow(scworkbench::wb_art_list(x$bundle)) == 3L)
  scworkbench::wb_session_discard(fit$session)
  project <- trajectory_s5_project(x$bundle, model_id)
  before_project <- trajectory_s5_state(x$bundle)
  options(scworkbench.promotion_workspace_publish_hook = function(...)
    stop("injected S5 curve publication failure"))
  rejected_project <- try(scworkbench::wb_session_promote_current(
    x$bundle, project$session, outputs = "curves"), silent = TRUE)
  options(scworkbench.promotion_workspace_publish_hook = NULL)
  stopifnot(inherits(rejected_project, "try-error"),
    identical(before_project, trajectory_s5_state(x$bundle)))
  retry <- scworkbench::wb_session_promote_current(x$bundle,
    project$session, outputs = "curves")
  stopifnot(nrow(retry$promoted) == 1L,
    identical(retry$promoted$status[[1L]], "ok"),
    nrow(scworkbench::wb_art_list(x$bundle)) == 4L)
  scworkbench::wb_session_discard(project$session)
  TRUE
}

test_trajectory_s5_monocle3_model_tamper_and_missing_pack_contract <- function() {
  x <- trajectory_s5_fixture()
  on.exit(trajectory_s5_close(x), add = TRUE)
  fit <- trajectory_s5_fit(x)
  model_id <- trajectory_s5_promote_fit(x, fit)
  scworkbench::wb_session_discard(fit$session)
  model <- scworkbench::wb_art_load(x$bundle, model_id,
    expected_kind = "model_ref")
  cds <- readRDS(file.path(model$model_tree_root,
    model$model_primary_path))
  lineage <- getFromNamespace(".trajectory_monocle3_lineage_paths",
    "scwbTrajectory")(cds)
  lineage$lineage_paths <- list()
  before_zero <- trajectory_s5_state(x$bundle)
  empty_relation <- try(getFromNamespace(
    ".trajectory_monocle3_branch_membership_prim",
    "scwbTrajectory")(cds, lineage), silent = TRUE)
  empty_curves <- try(getFromNamespace(
    ".trajectory_monocle3_lineage_curves_2d_prim",
    "scwbTrajectory")(cds, lineage, x$display), silent = TRUE)
  stopifnot(inherits(empty_relation, "try-error"),
    inherits(empty_curves, "try-error"),
    identical(before_zero, trajectory_s5_state(x$bundle)))
  source_path <- file.path(x$root,
    attr(scworkbench::wb_space(x$bundle, "latent"),
      "scworkbench_payload_read")$path)
  source_values <- rhdf5::h5read(source_path, "coordinates")
  altered_values <- source_values
  altered_values[1L, 1L] <- altered_values[1L, 1L] + 0.25
  rhdf5::h5write(altered_values, source_path, "coordinates")
  before_source <- trajectory_s5_state(x$bundle)
  source_rejected <- try(trajectory_s5_project(x$bundle, model_id),
    silent = TRUE)
  stopifnot(inherits(source_rejected, "try-error"),
    identical(before_source, trajectory_s5_state(x$bundle)))
  rhdf5::h5write(source_values, source_path, "coordinates")
  paths <- c(file.path(model$model_tree_root, model$model_primary_path),
    file.path(model$model_tree_root,
      model$model_source_linkage$selected_axis_evidence$relative_path))
  for (path in paths) {
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    mode <- file.info(path)$mode
    Sys.chmod(path, mode = "0644")
    writeBin(c(bytes, as.raw(0x7f)), path)
    before <- trajectory_s5_state(x$bundle)
    rejected <- try(trajectory_s5_project(x$bundle, model_id),
      silent = TRUE)
    stopifnot(inherits(rejected, "try-error"),
      identical(before, trajectory_s5_state(x$bundle)))
    writeBin(bytes, path)
    Sys.chmod(path, mode = mode)
  }
  registry <- getFromNamespace(".wb_operation_pack_registry",
    "scworkbench")
  saved_pack <- registry$packs$scwbtrajectory
  saved_operations <- registry$operations[saved_pack$operation_names]
  on.exit({
    registry$packs$scwbtrajectory <- saved_pack
    registry$operations[names(saved_operations)] <- saved_operations
  }, add = TRUE)
  for (name in saved_pack$operation_names)
    registry$operations[[name]] <- NULL
  registry$packs$scwbtrajectory <- NULL
  before <- trajectory_s5_state(x$bundle)
  missing <- try(trajectory_s5_project(x$bundle, model_id),
    silent = TRUE)
  stopifnot(inherits(missing, "try-error"),
    identical(before, trajectory_s5_state(x$bundle)),
    identical(scworkbench::wb_art_load(x$bundle, model_id,
      expected_kind = "model_ref")$engine, "monocle3"),
    is.list(scworkbench::wb_prov_bundle_dag(x$bundle)))
  registry$packs$scwbtrajectory <- saved_pack
  registry$operations[names(saved_operations)] <- saved_operations
  retry <- trajectory_s5_project(x$bundle, model_id)
  trajectory_s5_assert_curves(retry$run$outputs$curves, x$oracle)
  scworkbench::wb_session_discard(retry$session)
  TRUE
}
