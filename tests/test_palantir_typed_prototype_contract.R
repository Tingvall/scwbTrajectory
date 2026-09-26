palantir_proto_fixture <- function() {
  root <- tempfile("palantir-typed-prototype-")
  bundle <- scworkbench::wb_init_bundle(root)
  cells <- c(sprintf("c%02d", 0:9), sprintf("a%02d", 0:9),
    sprintf("b%02d", 0:9))
  coordinates <- t(vapply(seq_along(cells), function(j) {
    i <- (j - 1L) %% 10L
    if (j <= 10L) return(c(i / 9, 0.02 * sin(j), 0.02 * cos(j)))
    time <- (i + 1) / 10
    sign <- if (j <= 20L) 1 else -1
    c(1 + time, sign * time + 0.02 * sin(j), 0.02 * cos(j))
  }, numeric(3)))
  dimnames(coordinates) <- list(cells, c("d1", "d2", "d3"))
  plot <- coordinates[, 1:2, drop = FALSE]
  colnames(plot) <- c("x", "y")
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_add_space(bundle, "latent", coordinates)
  scworkbench::wb_add_space(bundle, "plot", plot)
  list(root = root, bundle = bundle, cells = cells,
    coordinates = coordinates, plot = plot)
}

palantir_proto_close <- function(x) {
  try(scworkbench::wb_close_bundle(x$bundle), silent = TRUE)
  unlink(x$root, recursive = TRUE, force = TRUE)
}

palantir_proto_fit_inputs <- function() list(fit_space = "latent",
  root_cell = "c00", terminal_cells = c("a09", "b09"),
  dims = as.integer(1:3), knn = 5L, n_diffusion_components = 3L,
  n_waypoints = 10L, use_early_cell_as_start = TRUE, seed = 7L)

palantir_proto_project_inputs <- function(model) list(model = model,
  plot_space = "plot", prob_min = 0.05, smoothness = 1,
  min_cells_per_branch = 3L, approx_points = 20L,
  pickle_execution_trust = "operator_verified_source_v1")

palantir_proto_state <- function(x) {
  con <- x$bundle$con
  tables <- sort(DBI::dbListTables(con))
  counts <- stats::setNames(vapply(tables, function(table) as.integer(
    DBI::dbGetQuery(con, paste0("SELECT COUNT(*) AS n FROM ",
      as.character(DBI::dbQuoteIdentifier(con, table))))$n[[1L]]),
    integer(1)), tables)
  files <- list.files(x$root, recursive = TRUE, all.files = TRUE,
    include.dirs = FALSE)
  db_file <- basename(DBI::dbGetInfo(con)$dbname)
  files <- sort(files[!basename(files) %in% c("manifest.json", db_file,
    paste0(db_file, c("-shm", "-wal", "-journal")))])
  payload <- stats::setNames(vapply(files, function(relative)
    digest::digest(file = file.path(x$root, relative), algo = "sha256"),
    character(1)), files)
  fingerprint <- tryCatch(getFromNamespace(
    ".wb_session_source_identity", "scworkbench")(
      x$bundle)$scientific_state_fingerprint,
    error = function(e) paste0("rejected:", conditionMessage(e)))
  list(artifacts = scworkbench::wb_art_list(x$bundle),
    fingerprint = fingerprint,
    counts = counts, payload = payload)
}

palantir_proto_tamper <- function(path, code) {
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  mode <- file.info(path)$mode
  on.exit({
    writeBin(bytes, path)
    Sys.chmod(path, mode = mode)
  }, add = TRUE)
  Sys.chmod(path, mode = "0644")
  writeBin(c(bytes, as.raw(0x7f)), path)
  force(code)
}

palantir_proto_oracle <- function() {
  path <- file.path(if (file.exists("tests/fixtures/phase5_s5_old_core_relation_oracle.json"))
    "tests" else ".", "fixtures", "phase5_s5_old_core_relation_oracle.json")
  stopifnot(file.exists(path))
  jsonlite::read_json(path, simplifyVector = FALSE)$oracles$palantir
}

palantir_proto_numeric_oracle <- function() {
  path <- file.path(if (file.exists("tests/fixtures/phase5_palantir_old_core_numeric_30.json"))
    "tests" else ".", "fixtures", "phase5_palantir_old_core_numeric_30.json")
  stopifnot(file.exists(path))
  jsonlite::read_json(path, simplifyVector = TRUE)
}

palantir_proto_numeric_equal <- function(actual, expected,
    tolerance = 1e-8) isTRUE(all.equal(as.numeric(actual),
      as.numeric(expected), tolerance = tolerance, check.attributes = FALSE))

test_palantir_typed_prototype_transient_fixture_contract <- function() {
  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbTrajectory"),
    simplifyVector = FALSE)
  names <- vapply(manifest$operations, function(operation)
    as.character(operation$name), character(1))
  stopifnot(identical(sum(names == "palantir_fit"), 1L),
    identical(sum(names == "palantir_project"), 1L),
    !any(c("fit_typed", "project_typed") %in% names),
    !any(c("palantir_fit", "palantir_project") %in%
      scworkbench::wb_list_ops()))
  TRUE
}

test_palantir_typed_prototype_fit_project_oracle_contract <-
    function(proto = NULL) {
  if (is.null(proto)) {
    proto <- palantir_proto_activate()
    on.exit(proto$cleanup(), add = TRUE)
  }
  current <- palantir_proto_fixture()
  on.exit(palantir_proto_close(current), add = TRUE)
  oracle <- palantir_proto_oracle()
  numeric_oracle <- palantir_proto_numeric_oracle()
  stopifnot(identical(as.character(unlist(oracle$source_cells,
    use.names = FALSE)), current$cells),
    identical(as.integer(oracle$row_count), 60L),
    !identical(oracle$terminal_content_identities$branch_summary$direct,
      oracle$terminal_content_identities$branch_summary$reopened))
  typed <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    current$bundle, proto$fit, inputs = palantir_proto_fit_inputs())
  closure <- sort(c("pseudotime", "branch_membership",
    "branch_summary", "entropy", "model"))
  stopifnot(identical(sort(names(typed$run$outputs)), closure))
  outputs <- typed$run$outputs
  for (name in c("pseudotime", "entropy")) {
    value <- outputs[[name]]
    stopifnot(identical(as.character(value$cell_id), current$cells),
      identical(as.character(value$value_type), "continuous"),
      palantir_proto_numeric_equal(value$value, numeric_oracle[[name]], 1e-8))
  }
  relation <- outputs$branch_membership
  rows <- oracle$ordered_rows
  stopifnot(length(relation$key) == 60L,
    identical(as.character(relation$key),
      vapply(rows, function(row) as.character(row$key), character(1))),
    identical(as.character(relation$cell_id),
      vapply(rows, function(row) as.character(row$cell_id), character(1))),
    palantir_proto_numeric_equal(relation$weight,
      as.numeric(vapply(rows, function(row)
        as.character(row$weight_hex), character(1))), 1e-10),
    all(is.finite(relation$weight)),
    all(relation$weight >= 0 & relation$weight <= 1))
  summary <- outputs$branch_summary$table
  expected_summary <- numeric_oracle$branch_summary
  stopifnot(identical(as.character(summary$relation_key),
      as.character(expected_summary$relation_key)),
    identical(as.character(summary$branch),
      as.character(expected_summary$branch)),
    identical(as.integer(summary$n_cells_prob_gt_0),
      as.integer(expected_summary$n_cells_prob_gt_0)),
    palantir_proto_numeric_equal(summary$mean_prob,
      expected_summary$mean_prob),
    palantir_proto_numeric_equal(summary$max_prob,
      expected_summary$max_prob),
    identical(as.character(summary$top_cell_id),
      as.character(expected_summary$top_cell_id)))
  expected <- palantir_proto_fit_inputs()
  stopifnot(identical(as.character(typed$run$params$root_cell),
      expected$root_cell),
    identical(as.character(typed$run$params$terminal_cells),
      expected$terminal_cells),
    identical(as.integer(typed$run$params$dims), expected$dims),
    identical(as.integer(typed$run$params$knn), expected$knn),
    identical(as.integer(typed$run$params$n_diffusion_components),
      expected$n_diffusion_components),
    identical(as.integer(typed$run$params$n_waypoints),
      expected$n_waypoints),
    identical(as.integer(typed$run$params$seed), expected$seed))
  stopifnot(identical(oracle$terminal_content_identities$branch_summary$direct$sha256,
      "2732769edcd0876158aa7fa3cda10522fa8681d8434cafa278ebfab4571626e6"),
    identical(oracle$terminal_content_identities$branch_summary$reopened$sha256,
      "5c2d4477cc88907c98093649497649ee322e91f59e621f7d712f7dbbf846cf05"))
  before <- palantir_proto_state(current)
  prior_hook <- getOption("scworkbench.promotion_workspace_publish_hook")
  on.exit(options(scworkbench.promotion_workspace_publish_hook =
    prior_hook), add = TRUE)
  options(scworkbench.promotion_workspace_publish_hook = function(...)
    stop("injected prototype publication failure"))
  rejected <- try(scworkbench::wb_session_promote_current(
    current$bundle, typed$session, outputs = closure), silent = TRUE)
  options(scworkbench.promotion_workspace_publish_hook = NULL)
  stopifnot(inherits(rejected, "try-error"),
    identical(before, palantir_proto_state(current)))
  promoted <- scworkbench::wb_session_promote_current(current$bundle,
    typed$session, outputs = closure)
  fit_details <- scworkbench::wb_prov_run_details(current$bundle,
    typed$run$run_id, promoted$provenance_snapshot_id)
  stopifnot(nrow(promoted$promoted) == 5L,
    all(promoted$promoted$status == "ok"),
    nrow(fit_details$outputs) == 5L,
    all(fit_details$outputs$verification_status == "content_verified"),
    !identical(before, palantir_proto_state(current)))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  model <- scworkbench::wb_art_load(current$bundle, model_id,
    expected_kind = "model_ref")
  stopifnot(identical(model$engine, "palantir"),
    identical(model$model_type, "palantir_fit"),
    identical(as.character(model$model_primary_path), "model.pkl"),
    file.exists(file.path(model$model_tree_root,
      model$model_primary_path)),
    identical(as.character(model$model_source_linkage$source_space_id),
      "latent"))
  stopifnot(startsWith(normalizePath(model$model_tree_root),
    paste0(normalizePath(current$root), .Platform$file.sep)))
  scworkbench::wb_session_discard(typed$session)
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    current$bundle)$ready))
  scworkbench::wb_close_bundle(current$bundle)
  current$bundle <- NULL
  relocated <- tempfile("palantir-typed-relocated-")
  fs::dir_copy(current$root, relocated)
  on.exit(unlink(relocated, recursive = TRUE, force = TRUE), add = TRUE)
  reopened <- scworkbench::wb_open_bundle(relocated)
  on.exit(try(scworkbench::wb_close_bundle(reopened), silent = TRUE),
    add = TRUE)
  stopifnot(isTRUE(scworkbench::wb_validate_bundle(relocated)))
  loaded <- scworkbench::wb_art_load(reopened, model_id,
    expected_kind = "model_ref")
  stopifnot(identical(loaded$engine, "palantir"),
    file.exists(file.path(loaded$model_tree_root,
      loaded$model_primary_path)))
  model_ref <- scworkbench::wb_ref_artifact(model_id, "model")
  projection <- scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), reopened, proto$project,
    inputs = palantir_proto_project_inputs(model_ref))
  curves <- projection$run$outputs$curves
  old_curves <- numeric_oracle$curves
  stopifnot(identical(names(projection$run$outputs), "curves"),
    identical(as.character(curves$key), as.character(old_curves$key)),
    identical(as.integer(curves$point), as.integer(old_curves$point)),
    palantir_proto_numeric_equal(curves$x, old_curves$x, 1e-10),
    palantir_proto_numeric_equal(curves$y, old_curves$y, 1e-10),
    identical(as.integer(projection$run$params$approx_points), 20L),
    identical(as.integer(projection$run$params$min_cells_per_branch), 3L),
    identical(as.character(projection$run$params$pickle_execution_trust),
      "operator_verified_source_v1"),
    length(unique(projection$run$outputs$curves$key)) == 2L,
    all(is.finite(projection$run$outputs$curves$x)),
    all(is.finite(projection$run$outputs$curves$y)))
  projection_promoted <- scworkbench::wb_session_promote_current(reopened,
    projection$session, outputs = "curves")
  project_details <- scworkbench::wb_prov_run_details(reopened,
    projection$run$run_id, projection_promoted$provenance_snapshot_id)
  stopifnot(nrow(projection_promoted$promoted) == 1L,
    identical(projection_promoted$promoted$status, "ok"),
    nrow(project_details$outputs) == 1L,
    identical(project_details$outputs$verification_status,
      "content_verified"))
  curve_id <- projection_promoted$promoted$artifact_id[[1L]]
  scworkbench::wb_session_discard(projection$session)
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(reopened)$ready))
  scworkbench::wb_close_bundle(reopened)
  reopened <- scworkbench::wb_open_bundle(relocated, read_only = TRUE)
  stopifnot(isTRUE(scworkbench::wb_bundle_validation_report(relocated,
    "recipient")$ok),
    inherits(scworkbench::wb_art_load(reopened, curve_id,
      expected_kind = "path_2d_family"), "wb_primitive"))
  TRUE
}

test_palantir_typed_prototype_private_failure_zero_residue_contract <-
    function(proto = NULL) {
  if (is.null(proto)) {
    proto <- palantir_proto_activate()
    on.exit(proto$cleanup(), add = TRUE)
  }
  x <- palantir_proto_fixture()
  on.exit(palantir_proto_close(x), add = TRUE)
  before <- palantir_proto_state(x)
  prior <- getOption("scwbTrajectory.palantir_typed_python")
  on.exit(options(scwbTrajectory.palantir_typed_python = prior), add = TRUE)
  options(scwbTrajectory.palantir_typed_python =
    "/nonexistent/palantir-python")
  rejected <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, proto$fit,
    inputs = palantir_proto_fit_inputs()), silent = TRUE)
  stopifnot(inherits(rejected, "try-error"),
    identical(before, palantir_proto_state(x)),
    nrow(scworkbench::wb_art_list(x$bundle)) == 0L)
  TRUE
}

test_palantir_typed_prototype_source_model_tamper_contract <-
    function(proto = NULL) {
  if (is.null(proto)) {
    proto <- palantir_proto_activate()
    on.exit(proto$cleanup(), add = TRUE)
  }
  x <- palantir_proto_fixture()
  on.exit(palantir_proto_close(x), add = TRUE)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, proto$fit, inputs = palantir_proto_fit_inputs())
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("pseudotime", "branch_membership",
      "branch_summary", "entropy", "model"))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  scworkbench::wb_session_discard(fit$session)
  model <- scworkbench::wb_art_load(x$bundle, model_id,
    expected_kind = "model_ref")
  reference <- scworkbench::wb_ref_artifact(model_id, "model")
  request <- palantir_proto_project_inputs(reference)
  model_path <- file.path(model$model_tree_root,
    model$model_primary_path)
  palantir_proto_tamper(model_path, {
    before <- palantir_proto_state(x)
    rejected <- try(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), x$bundle, proto$project,
      inputs = request), silent = TRUE)
    stopifnot(inherits(rejected, "try-error"),
      !grepl("Private Palantir computation failed", as.character(rejected),
        fixed = TRUE), identical(before, palantir_proto_state(x)))
  })
  source_row <- DBI::dbGetQuery(x$bundle$con,
    "SELECT path FROM space_payloads WHERE space_id = ?;",
    params = list("latent"))
  stopifnot(nrow(source_row) == 1L)
  source_path <- file.path(x$root, source_row$path[[1L]])
  palantir_proto_tamper(source_path, {
    before <- palantir_proto_state(x)
    rejected <- try(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), x$bundle, proto$project,
      inputs = request), silent = TRUE)
    stopifnot(inherits(rejected, "try-error"),
      !grepl("Private Palantir computation failed", as.character(rejected),
        fixed = TRUE), identical(before, palantir_proto_state(x)))
  })
  prior <- getOption("scwbTrajectory.palantir_typed_python")
  on.exit(options(scwbTrajectory.palantir_typed_python = prior), add = TRUE)
  options(scwbTrajectory.palantir_typed_python =
    "/nonexistent/palantir-python")
  before <- palantir_proto_state(x)
  rejected <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, proto$project,
    inputs = request), silent = TRUE)
  stopifnot(inherits(rejected, "try-error"),
    identical(before, palantir_proto_state(x)))
  options(scwbTrajectory.palantir_typed_python = prior)
  retry <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, proto$project, inputs = request)
  stopifnot(length(retry$run$outputs$curves$key) == 40L,
    nrow(scworkbench::wb_art_list(x$bundle)) == 5L)
  scworkbench::wb_session_discard(retry$session)
  TRUE
}
