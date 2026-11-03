trajectory_s3_fixture <- function() {
  trajectory_activate_pack()
  root <- tempfile("trajectory-s3-")
  bundle <- scworkbench::wb_init_bundle(root)
  cells <- paste0("cell_", seq_len(120L))
  groups <- rep(c("root", "left", "right"), each = 40L)
  coords <- cbind(
    PC1 = rep(c(0, -2, 2), each = 40L) +
      rep(seq(-.2, .2, length.out = 40L), 3L),
    PC2 = rep(c(0, 2, 2), each = 40L) +
      sin(seq_len(120L) * .37) * .12,
    PC3 = rep(c(0, -1, 1), each = 40L) +
      cos(seq_len(120L) * .23) * .11)
  rownames(coords) <- cells
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_add_space(bundle, "fit", coords)
  plot <- coords[, 1:2, drop = FALSE]
  colnames(plot) <- c("x", "y")
  scworkbench::wb_add_space(bundle, "plot", plot)
  cluster <- scworkbench::wb_prim_cell_vector(stats::setNames(groups,
    cells), value_type = "categorical")
  list(root = root, bundle = bundle, cells = cells, groups = groups,
    coords = coords, plot = plot, cluster = cluster)
}

trajectory_s3_close <- function(x) {
  try(scworkbench::wb_close_bundle(x$bundle), silent = TRUE)
  unlink(x$root, recursive = TRUE, force = TRUE)
}

trajectory_s3_fit <- function(x, extra = list()) {
  scworkbench::wb_session_run_op(scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::slingshot_fit", inputs = c(list(
      fit_space = "fit", clusters = x$cluster,
      start.clus = "root", end.clus = c("left", "right"),
      min_cells = 1L), extra))
}

trajectory_s3_promote_fit <- function(x, fit) {
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("model", "pseudotime", "lineage_weights"))
  stopifnot(nrow(promoted$promoted) == 3L,
    all(promoted$promoted$status == "ok"))
  promoted$promoted$artifact_id[promoted$promoted$output == "model"][[1L]]
}

trajectory_s3_state <- function(bundle) {
  list(artifacts = scworkbench::wb_art_list(bundle),
    fingerprint = getFromNamespace(".wb_session_source_identity",
      "scworkbench")(bundle)$scientific_state_fingerprint)
}

trajectory_s3_contained_model <- function(result) {
  reference <- result$run$outputs$model
  path <- file.path(result$session$staging_bundle$dir,
    reference$model_tree_root, reference$model_primary_path)
  stopifnot(file.exists(path))
  readRDS(path)
}

test_trajectory_s3_registered_fit_contract <- function() {
  trajectory_activate_pack()
  manifest <- scworkbench::wb_operation_pack_inspect(
    "scwbtrajectory")$manifest
  operations <- stats::setNames(manifest$operations,
    vapply(manifest$operations, `[[`, "", "name"))
  for (name in c("slingshot", "slingshot_fit",
    "slingshot_lineage_weights", "slingshot_project")) {
    entry <- operations[[name]]
    outputs <- vapply(entry$outputs, `[[`, "", "name")
    stopifnot(all(vapply(entry$promotion_closure, function(closure)
      setequal(unlist(closure), outputs), logical(1))),
      identical(names(entry$replay_equivalence), outputs))
  }
  root <- tempfile("trajectory-s3-fit-")
  bundle <- scworkbench::wb_init_bundle(root)
  on.exit({
    try(scworkbench::wb_close_bundle(bundle), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  cells <- paste0("cell_", seq_len(120L))
  groups <- rep(c("root", "left", "right"), each = 40L)
  coords <- cbind(
    PC1 = rep(c(0, -2, 2), each = 40L) +
      rep(seq(-.2, .2, length.out = 40L), 3L),
    PC2 = rep(c(0, 2, 2), each = 40L) +
      sin(seq_len(120L) * .37) * .12,
    PC3 = rep(c(0, -1, 1), each = 40L) +
      cos(seq_len(120L) * .23) * .11)
  rownames(coords) <- cells
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_add_space(bundle, "fit", coords)
  plot <- coords[, 1:2, drop = FALSE]
  colnames(plot) <- c("x", "y")
  scworkbench::wb_add_space(bundle, "plot", plot)
  cluster <- scworkbench::wb_prim_cell_vector(stats::setNames(groups,
    cells), value_type = "categorical")
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    bundle, "scwbtrajectory::slingshot_fit", inputs = list(
      fit_space = "fit", clusters = cluster,
      start.clus = "root", end.clus = c("left", "right"),
      min_cells = 1L))
  stopifnot(identical(sort(names(fit$run$outputs)),
    c("lineage_weights", "model", "pseudotime")),
    length(unique(fit$run$outputs$pseudotime$cell_id)) <= 120L)
  promoted <- scworkbench::wb_session_promote_current(bundle, fit$session,
    outputs = c("model", "pseudotime", "lineage_weights"))
  stopifnot(nrow(promoted$promoted) == 3L,
    all(promoted$promoted$status == "ok"))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  scworkbench::wb_session_discard(fit$session)
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(bundle)$ready))
  scworkbench::wb_close_bundle(bundle)
  relocated <- tempfile("trajectory-s3-relocated-")
  fs::dir_copy(root, relocated)
  on.exit(unlink(relocated, recursive = TRUE, force = TRUE), add = TRUE)
  reopened <- scworkbench::wb_open_bundle(relocated)
  on.exit(try(scworkbench::wb_close_bundle(reopened), silent = TRUE),
    add = TRUE)
  stopifnot(isTRUE(scworkbench::wb_validate_bundle(relocated)))
  model_ref <- scworkbench::wb_ref_artifact(model_id, "model")
  weights <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    reopened, "scwbtrajectory::slingshot_lineage_weights",
    inputs = list(model = model_ref, min_cells = 1L))
  stopifnot(length(weights$run$outputs$lineage_weights$key) > 0L,
    identical(as.character(weights$run$outputs$lineage_weights$meta$fit_space_id),
      "fit"))
  project <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    reopened, "scwbtrajectory::slingshot_project", inputs = list(
      model = model_ref, plot_space = "plot", approx_points = 20L,
      min_cells = 1L))
  stopifnot(length(project$run$outputs$curves$key) > 0L,
    identical(as.character(project$run$outputs$curves$meta$fit_space_id),
      "fit"))
  scworkbench::wb_session_discard(weights$session)
  scworkbench::wb_session_discard(project$session)
  TRUE
}

test_trajectory_s3_exact_name_ownership_contract <- function() {
  script <- tempfile("trajectory-s3-ownership-", fileext = ".R")
  log <- tempfile("trajectory-s3-ownership-", fileext = ".log")
  on.exit(unlink(c(script, log)), add = TRUE)
  writeLines(c(
    paste0(".libPaths(", paste(capture.output(dput(.libPaths())),
      collapse = ""), ")"),
    "library(scworkbench)",
    "names <- c('slingshot', 'slingshot_fit',",
    "  'slingshot_lineage_weights', 'slingshot_project')",
    "before <- wb_list_ops()",
    "stopifnot(all(!names %in% before),",
    "  all(!paste0('scwbtrajectory::', names) %in% before),",
    "  sum(before == 'slingshot_adopt') == 1L)",
    "manifest <- system.file('scworkbench', 'operation-pack.json',",
    "  package = 'scwbTrajectory')",
    "wb_operation_pack_activate(manifest, lib.loc = .libPaths())",
    "after <- wb_list_ops()",
    "fit_interface <- wb_op_interface('scwbtrajectory::slingshot_fit')",
    "stopifnot(all(c('start.clus', 'end.clus', 'exclude_cells',",
    "  'include_clusters', 'exclude_clusters', 'stretch', 'reweight',",
    "  'shrink', 'maxit', 'omega', 'omega_scale') %in% vapply(fit_interface$params,",
    "  `[[`, character(1), 'name')))",
    "stopifnot(all(!names %in% after),",
    "  all(vapply(paste0('scwbtrajectory::', names), function(x)",
    "    sum(after == x) == 1L, logical(1))),",
    "  sum(after == 'slingshot_adopt') == 1L)",
    "cat('RESULT S3 four pack names once; core adoption once\\n')"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    args = shQuote(script), stdout = log, stderr = log)
  if (!identical(as.integer(status), 0L))
    stop(paste(readLines(log, warn = FALSE), collapse = "\n"),
      call. = FALSE)
  stopifnot(any(readLines(log, warn = FALSE) ==
    "RESULT S3 four pack names once; core adoption once"))
  TRUE
}

test_trajectory_s3_omega_interface_and_default_contract <- function() {
  x <- trajectory_s3_fixture(); on.exit(trajectory_s3_close(x), add = TRUE)
  interface <- scworkbench::wb_op_interface(
    "scwbtrajectory::slingshot_fit")
  parameters <- stats::setNames(interface$params,
    vapply(interface$params, `[[`, character(1), "name"))
  stopifnot(identical(parameters$omega$type, "bool_or_dbl"),
    identical(parameters$omega$min, 0),
    identical(parameters$omega_scale$type, "dbl"),
    identical(parameters$omega_scale$min, 0))

  omitted <- trajectory_s3_fit(x)
  explicit <- trajectory_s3_fit(x, list(omega = FALSE,
    omega_scale = 1.5))
  omitted_model <- trajectory_s3_contained_model(omitted)
  explicit_model <- trajectory_s3_contained_model(explicit)
  same_family <- function(left, right) identical(left$key, right$key) &&
    identical(left$cell_id, right$cell_id) &&
    identical(left$value, right$value) &&
    identical(left$value_type, right$value_type)
  stopifnot(
    same_family(omitted$run$outputs$pseudotime,
      explicit$run$outputs$pseudotime),
    same_family(omitted$run$outputs$lineage_weights,
      explicit$run$outputs$lineage_weights),
    identical(slingshot::slingLineages(omitted_model),
      slingshot::slingLineages(explicit_model)),
    identical(slingshot::slingParams(omitted_model)$omega, FALSE),
    identical(slingshot::slingParams(omitted_model)$omega_scale, 1.5),
    identical(slingshot::slingParams(explicit_model)$omega, FALSE),
    identical(slingshot::slingParams(explicit_model)$omega_scale, 1.5))

  requested <- trajectory_s3_fit(x, list(omega = 0.8,
    omega_scale = 3.2))
  requested_params <- slingshot::slingParams(
    trajectory_s3_contained_model(requested))
  stopifnot(identical(requested_params$omega, 0.8),
    identical(requested_params$omega_scale, 3.2))

  before <- trajectory_s3_state(x$bundle)
  for (extra in list(list(omega = -0.1), list(omega_scale = -1),
      list(omega = "TRUE"), list(omega = NA))) {
    rejected <- try(trajectory_s3_fit(x, extra), silent = TRUE)
    stopifnot(inherits(rejected, "try-error"),
      identical(before, trajectory_s3_state(x$bundle)))
  }
  invisible(scworkbench::wb_session_discard(omitted$session))
  invisible(scworkbench::wb_session_discard(explicit$session))
  invisible(scworkbench::wb_session_discard(requested$session))
  TRUE
}

test_trajectory_s3_omega_outgroup_fixture_contract <- function() {
  trajectory_activate_pack()
  root <- tempfile("trajectory-s3-outgroup-")
  bundle <- scworkbench::wb_init_bundle(root)
  on.exit({
    try(scworkbench::wb_close_bundle(bundle), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  centers <- rbind(root = c(0, 0), mid = c(1, 0), end = c(2, 0),
    outgroup = c(30, 30))
  angles <- seq(0, 2 * pi, length.out = 13L)[-13L]
  coordinates <- do.call(rbind, lapply(seq_len(nrow(centers)), function(i)
    cbind(D1 = centers[i, 1L] + 0.03 * cos(angles),
      D2 = centers[i, 2L] + 0.02 * sin(angles))))
  cells <- paste0("cell_", seq_len(nrow(coordinates)))
  groups <- rep(rownames(centers), each = length(angles))
  rownames(coordinates) <- cells
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_add_space(bundle, "fit", coordinates)
  clusters <- scworkbench::wb_prim_cell_vector(
    stats::setNames(groups, cells), value_type = "categorical")
  run <- function(omega, omega_scale = NULL) {
    inputs <- list(fit_space = "fit", clusters = clusters,
      start.clus = "root", omega = omega, min_cells = 1L)
    if (!is.null(omega_scale)) inputs$omega_scale <- omega_scale
    scworkbench::wb_session_run_op(scworkbench::wb_session_new(), bundle,
      "scwbtrajectory::slingshot_fit", inputs = inputs)
  }
  connected <- run(FALSE)
  separated <- run(TRUE, 1.5)
  connected_model <- trajectory_s3_contained_model(connected)
  separated_model <- trajectory_s3_contained_model(separated)
  stopifnot(
    identical(unname(slingshot::slingLineages(connected_model)[[1L]]),
      c("root", "mid", "end", "outgroup")),
    identical(unname(slingshot::slingLineages(separated_model)[[1L]]),
      c("root", "mid", "end")),
    identical(unname(slingshot::slingLineages(separated_model)[[2L]]),
      "outgroup"),
    identical(slingshot::slingParams(separated_model)$omega, TRUE),
    identical(slingshot::slingParams(separated_model)$omega_scale, 1.5))
  invisible(scworkbench::wb_session_discard(connected$session))
  invisible(scworkbench::wb_session_discard(separated$session))
  TRUE
}

test_trajectory_s3_scientific_axes_and_values_contract <- function() {
  x <- trajectory_s3_fixture(); on.exit(trajectory_s3_close(x), add = TRUE)
  legacy <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbtrajectory::slingshot", inputs = list(
      space_2d = "plot", clusters = x$cluster,
      start.clus = "root", end.clus = c("left", "right"),
      min_cells = 1L))
  fit <- trajectory_s3_fit(x, list(exclude_cells = x$cells[1:5]))
  stopifnot(identical(sort(names(legacy$run$outputs)),
    c("curves", "pseudotime")),
    identical(sort(names(fit$run$outputs)),
      c("lineage_weights", "model", "pseudotime")))
  fit_cells <- x$cells[-(1:5)]
  model <- fit$run$outputs$model
  stopifnot(identical(as.character(model$cells), fit_cells))
  for (name in c("pseudotime", "lineage_weights")) {
    family <- fit$run$outputs[[name]]
    stopifnot(identical(family$kind, "cell_vector_family"),
      identical(family$value_type, "continuous"),
      identical(unique(as.character(family$key)), c("L1", "L2")),
      all(family$cell_id %in% fit_cells),
      identical(as.character(family$cell_id[family$key == "L1"]),
        fit_cells[fit_cells %in% family$cell_id[family$key == "L1"]]))
  }
  pt <- fit$run$outputs$pseudotime
  stopifnot(any(vapply(unique(pt$key), function(key)
      sum(pt$key == key) < length(fit_cells), logical(1))),
    all(is.finite(pt$value)))
  weight <- fit$run$outputs$lineage_weights
  stopifnot(all(is.finite(weight$value)),
    all(weight$value >= 0 & weight$value <= 1))
  curves <- legacy$run$outputs$curves
  stopifnot(identical(curves$kind, "path_2d_family"),
    identical(unique(curves$key), c("L1", "L2")),
    all(is.finite(curves$x)), all(is.finite(curves$y)))
  scworkbench::wb_session_discard(legacy$session)
  scworkbench::wb_session_discard(fit$session)
  TRUE
}

test_trajectory_s3_publication_rollback_retry_contract <- function() {
  x <- trajectory_s3_fixture(); on.exit({
    options(scworkbench.promotion_workspace_publish_hook = NULL)
    trajectory_s3_close(x)
  }, add = TRUE)
  fit <- trajectory_s3_fit(x)
  before <- trajectory_s3_state(x$bundle)
  options(scworkbench.promotion_workspace_publish_hook = function(...)
    stop("injected S3 publication failure"))
  rejected <- try(trajectory_s3_promote_fit(x, fit), silent = TRUE)
  options(scworkbench.promotion_workspace_publish_hook = NULL)
  stopifnot(inherits(rejected, "try-error"),
    grepl("injected S3 publication failure", as.character(rejected),
      fixed = TRUE),
    identical(before, trajectory_s3_state(x$bundle)))
  model_id <- trajectory_s3_promote_fit(x, fit)
  stopifnot(nzchar(model_id), nrow(scworkbench::wb_art_list(x$bundle)) == 3L)
  scworkbench::wb_session_discard(fit$session)
  TRUE
}

test_trajectory_s3_model_and_axis_tampering_contract <- function() {
  x <- trajectory_s3_fixture(); on.exit(trajectory_s3_close(x), add = TRUE)
  fit <- trajectory_s3_fit(x)
  model_id <- trajectory_s3_promote_fit(x, fit)
  scworkbench::wb_session_discard(fit$session)
  model <- scworkbench::wb_art_load(x$bundle, model_id,
    expected_kind = "model_ref")
  reference <- scworkbench::wb_ref_artifact(model_id, "model")
  inputs <- list(model = reference, min_cells = 1L)
  paths <- c(file.path(model$model_tree_root, model$model_primary_path),
    file.path(model$model_tree_root,
      model$model_source_linkage$selected_axis_evidence$relative_path))
  for (path in paths) {
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    mode <- file.info(path)$mode
    Sys.chmod(path, mode = "0644")
    writeBin(c(bytes, as.raw(0x7f)), path)
    mutated <- trajectory_s3_state(x$bundle)
    for (name in c("slingshot_lineage_weights", "slingshot_project")) {
      request <- inputs
      if (identical(name, "slingshot_project"))
        request <- c(request, list(plot_space = "plot"))
      rejected <- try(scworkbench::wb_session_run_op(
        scworkbench::wb_session_new(), x$bundle,
        paste0("scwbtrajectory::", name), inputs = request), silent = TRUE)
      stopifnot(inherits(rejected, "try-error"),
        identical(mutated, trajectory_s3_state(x$bundle)))
    }
    writeBin(bytes, path)
    Sys.chmod(path, mode = mode)
  }
  retry <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbtrajectory::slingshot_lineage_weights",
    inputs = inputs)
  stopifnot(nrow(scworkbench::wb_art_list(x$bundle)) == 3L,
    length(retry$run$outputs$lineage_weights$key) > 0L)
  scworkbench::wb_session_discard(retry$session)
  TRUE
}

test_trajectory_s3_pack_free_retained_inspection_contract <- function() {
  x <- trajectory_s3_fixture(); on.exit(trajectory_s3_close(x), add = TRUE)
  fit <- trajectory_s3_fit(x)
  model_id <- trajectory_s3_promote_fit(x, fit)
  scworkbench::wb_session_discard(fit$session)
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
  before <- trajectory_s3_state(x$bundle)
  missing <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::slingshot_lineage_weights",
    inputs = list(model = scworkbench::wb_ref_artifact(model_id, "model"))),
    silent = TRUE)
  stopifnot(inherits(missing, "try-error"),
    identical(before, trajectory_s3_state(x$bundle)),
    identical(scworkbench::wb_art_load(x$bundle, model_id,
      expected_kind = "model_ref")$engine, "slingshot"),
    is.list(scworkbench::wb_prov_bundle_dag(x$bundle)),
    "slingshot_adopt" %in% scworkbench::wb_list_ops())
  TRUE
}
