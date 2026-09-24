test_dir <- if (file.exists("helper_s2_fixture.R")) "." else "tests"
source(file.path(test_dir, "test_s2_model_execution_contract.R"),
  local = .GlobalEnv)

# This fixture uses the accepted core G1 test pack to seal and reopen an
# actual Slingshot RDS. It does not register any S2 consumer.
.trajectory_s2_g1_fixture_source <- function(directory) {
  fixture_source <- file.path(directory, "fixtures",
    "s2-g1-derived-axis")
  fixture_files <- c("DESCRIPTION", "NAMESPACE", "R/register.R",
    "inst/scworkbench/operation-pack.json")
  fixture_sha256 <- c(
    "0c4b651150ceb54ebe978bf6c27f4a13294386829cccf7939e9bd3b020e9366f",
    "b1f21ca084303169f23123f59d2ed35079ff820fcd99ab5bbfb906233410e72e",
    "340de2b30eea159118f8ab7b858d8b21d3fb4ccbafb660cf03f4232b9b5bebaa",
    "af1c672ef2f4d5e75b77ec4811ac4970b3d5b802690995d09694770ca12c8fa8")
  paths <- file.path(fixture_source, fixture_files)
  if (!all(file.exists(paths)))
    stop("Real G1 proof fixture is missing: ", fixture_source,
      call. = FALSE)
  observed_sha256 <- vapply(paths, digest::digest, character(1),
    algo = "sha256", file = TRUE)
  if (!identical(unname(observed_sha256), fixture_sha256))
    stop("Real G1 proof fixture identity differs from the accepted G1 pin.",
      call. = FALSE)
  fixture_source
}

test_trajectory_s2_real_g1_fixture_pin_contract <- function() {
  absent <- tempfile("forge2-g1-absent-"); dir.create(absent)
  on.exit(unlink(absent, recursive = TRUE, force = TRUE), add = TRUE)
  missing <- tryCatch(.trajectory_s2_g1_fixture_source(absent),
    error = identity)
  stopifnot(inherits(missing, "error"),
    grepl("fixture is missing", conditionMessage(missing), fixed = TRUE))
  fixture <- file.path(absent, "fixtures", "s2-g1-derived-axis")
  dir.create(dirname(fixture), recursive = TRUE)
  fs::dir_copy(file.path(test_dir, "fixtures", "s2-g1-derived-axis"),
    fixture)
  write("tamper", file = file.path(fixture, "R", "register.R"),
    append = TRUE)
  changed <- tryCatch(.trajectory_s2_g1_fixture_source(absent),
    error = identity)
  stopifnot(inherits(changed, "error"),
    grepl("fixture identity differs", conditionMessage(changed),
      fixed = TRUE))
  TRUE
}

test_trajectory_s2_real_g1_reopen_contract <- function() {
  core <- asNamespace("scworkbench")
  required <- c(".wb_operation_pack_v3_prepare_self_linked_model_reads",
    ".wb_operation_pack_v3_model_adapters")
  if (!all(vapply(required, exists, logical(1), envir = core,
      inherits = FALSE)))
    stop("Real G1 proof requires the accepted derived fitted-axis core.",
      call. = FALSE)
  fixture_source <- .trajectory_s2_g1_fixture_source(test_dir)
  if (!requireNamespace("slingshot", quietly = TRUE) ||
      !requireNamespace("fs", quietly = TRUE))
    stop("Real G1 proof requires slingshot and fs.", call. = FALSE)
  library <- tempfile("forge2-g1-lib-"); dir.create(library)
  install_log <- tempfile("forge2-g1-install-", fileext = ".log")
  status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL",
    paste0("--library=", shQuote(library)), shQuote(fixture_source)),
    stdout = install_log, stderr = install_log)
  if (!identical(status, 0L)) stop(paste(readLines(install_log,
    warn = FALSE), collapse = "\n"))
  old_paths <- .libPaths(); .libPaths(c(library, old_paths))
  root <- tempfile("forge2-g1-bundle-")
  relocated <- tempfile("forge2-g1-relocated-")
  stage <- tempfile("forge2-g1-private-")
  bundle <- NULL; reopened <- NULL; x <- NULL
  on.exit({
    try(scworkbench::wb_operation_pack_unload("scwbv3phase5"), silent = TRUE)
    try(scworkbench::wb_close_bundle(reopened), silent = TRUE)
    try(scworkbench::wb_close_bundle(bundle), silent = TRUE)
    if (!is.null(x)) try(trajectory_close(x), silent = TRUE)
    .libPaths(old_paths)
    unlink(c(library, root, relocated, stage, install_log),
      recursive = TRUE, force = TRUE)
  }, add = TRUE)
  loadNamespace("scwbV3Phase5", lib.loc = library)
  fit_runner <- function(readers, writers, params, context) {
    source <- readers$source()
    fitted <- as.character(params$fitted_cells)
    writers$declare_derived_fitted_axis(fitted)
    writers$write_terminal_primitive("embedding",
      scworkbench::wb_prim_cell_vector(stats::setNames(seq_along(fitted),
        fitted), value_type = "continuous"))
    writers$write_terminal_primitive("weights",
      scworkbench::wb_prim_cell_vector_family(
        rep("lineage_a", length(fitted)), fitted,
        rep(1, length(fitted)), value_type = "continuous"))
    tree <- writers$open_model_tree("model")
    saveRDS(.trajectory_s2_synthetic_slingshot(fitted),
      tree$path(tree$primary_path))
    model <- scworkbench::wb_prim_model_ref("fixture.derived",
      as.character(context$input_ids$source)[1L], engine = "fixture")
    model$model_type <- "derived_axis_model"
    tree$finish(model)
    invisible(list(status = "ok"))
  }
  assignInNamespace(".fit_derived_axis_model", fit_runner,
    ns = "scwbV3Phase5")
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbV3Phase5", lib.loc = library)
  scworkbench::wb_operation_pack_activate(manifest, lib.loc = library)
  x <- trajectory_s2_fixture()
  fitted <- x$cells[1:12]
  bundle <- scworkbench::wb_init_bundle(root)
  scworkbench::wb_define_cells(bundle, x$cells)
  scworkbench::wb_add_space(bundle, "pca", matrix(
    seq_len(length(x$cells) * 3L), nrow = length(x$cells),
    dimnames = list(x$cells, c("PC1", "PC2", "PC3"))))
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    bundle, "scwbv3phase5::fit_derived_axis_model", inputs = list(
      source = "pca", dims = c(1L, 3L), fitted_cells = fitted))
  promoted <- scworkbench::wb_session_promote_current(bundle, fit$session,
    outputs = c("embedding", "model", "weights"),
    created_by = "forge2-g1-reopen")
  stopifnot(nrow(promoted$promoted) == 3L,
    isTRUE(scworkbench::wb_finalize_recipient_bundle(bundle)$ready))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  scworkbench::wb_close_bundle(bundle); bundle <- NULL
  fs::dir_copy(root, relocated)
  reopened <- scworkbench::wb_open_bundle(relocated, read_only = TRUE)
  durable <- scworkbench::wb_art_load(reopened, model_id,
    expected_kind = "model_ref")
  registry <- get(".wb_operation_pack_registry", envir = core)
  consume <- registry$operations[[
    "scwbv3phase5::consume_self_linked_model"]]$definition$v3_definition
  descriptor <- getFromNamespace(
    ".wb_operation_pack_v2_input_descriptor", "scworkbench")(
      "model", durable, consume$inputs[[1L]])$descriptor
  dir.create(stage)
  prepare <- getFromNamespace(
    ".wb_operation_pack_v3_prepare_self_linked_model_reads",
    "scworkbench")
  prepared <- prepare(reopened, consume, list(model = durable),
    list(model = descriptor), stage)
  stopifnot(identical(rownames(prepared$values[[1L]]$fitted_space),
    fitted))
  request <- list(outputs = list(),
    self_linked_model_payload_grants = unname(prepared$grants))
  attr(request, "operation_pack_v3_self_linked_model_read_values") <-
    prepared$values
  state <- new.env(parent = emptyenv()); state$used <- FALSE
  adapter <- getFromNamespace(".wb_operation_pack_v3_model_adapters",
    "scworkbench")(list(state = state, readers = list(), writers = list()),
      request, consume, stage)
  observed <- .trajectory_s2_run_model_direct(x, "smooth",
    params = list(lineage = "L1", features = x$features[1:3],
      n_bins = 4L, min_cells_per_bin = 1L),
    effective_cells = fitted,
    model_payload_reader = adapter$readers$fitted_model)
  stopifnot(!inherits(observed$value, "error"),
    identical(names(observed$outputs), c("trend_table", "fit_summary")),
    identical(state$model_payload_calls$fitted_model, 1L))
  reference <- jsonlite::read_json(file.path(test_dir, "fixtures",
    "trajectory_s2_parameter_reference.json"), simplifyVector = FALSE)
  trajectory_s2_compare_outputs(x$bundle, "smooth_reopened",
    list(outputs = observed$outputs), reference)
  stopifnot(length(durable$model_tree_root) == 1L,
    length(prepared$values[[1L]]$manifest$primary_path) == 1L,
    length(durable$model_source_linkage$selected_axis_evidence$relative_path) == 1L)
  retained_state <- function() {
    tree <- durable$model_tree_root
    entries <- sort(list.files(tree, recursive = TRUE,
      all.files = TRUE, no.. = TRUE))
    paths <- file.path(tree, entries)
    files <- paths[!file.info(paths)$isdir]
    list(artifacts = scworkbench::wb_art_list(reopened),
      payload = stats::setNames(vapply(files, digest::digest,
        character(1), algo = "sha256", file = TRUE),
        substring(files, nchar(tree) + 2L)))
  }
  # Both attacks target the relocated retained model, before a private read.
  test_attack <- function(path) {
    before <- retained_state()
    before_calls <- state$model_payload_calls
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    mode <- file.info(path)$mode
    Sys.chmod(path, mode = "0644")
    restored <- FALSE
    restore <- function() {
      connection <- file(path, "wb")
      writeBin(bytes, connection); close(connection)
      Sys.chmod(path, mode = mode)
      restored <<- TRUE
    }
    on.exit(if (!restored) restore(), add = TRUE)
    connection <- file(path, "ab")
    writeBin(charToRaw("tamper"), connection); close(connection)
    tampered <- retained_state()
    stopifnot(identical(before$artifacts, tampered$artifacts),
      !identical(before$payload, tampered$payload))
    reject_stage <- tempfile("forge2-g1-reject-")
    on.exit(unlink(reject_stage, recursive = TRUE, force = TRUE), add = TRUE)
    rejected <- tryCatch(prepare(reopened, consume,
      list(model = durable), list(model = descriptor),
      reject_stage), error = identity)
    stopifnot(inherits(rejected, "error"),
      identical(tampered, retained_state()),
      identical(before_calls, state$model_payload_calls))
    restore()
    stopifnot(identical(before, retained_state()),
      identical(file.info(path)$mode, mode))
  }
  test_attack(file.path(durable$model_tree_root,
    prepared$values[[1L]]$manifest$primary_path))
  test_attack(file.path(durable$model_tree_root,
    durable$model_source_linkage$selected_axis_evidence$relative_path))
  TRUE
}
