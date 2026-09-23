source(file.path(test_dir, "test_s2_model_execution_contract.R"),
  local = .GlobalEnv)

# This fixture uses the accepted core G1 test pack to seal and reopen an
# actual Slingshot RDS. It does not register any S2 consumer.
test_trajectory_s2_real_g1_reopen_contract <- function() {
  core <- asNamespace("scworkbench")
  required <- c(".wb_operation_pack_v3_prepare_self_linked_model_reads",
    ".wb_operation_pack_v3_model_adapters")
  if (!all(vapply(required, exists, logical(1), envir = core,
      inherits = FALSE)))
    return(structure(TRUE, skip = "accepted G1 core unavailable"))
  fixture_source <- file.path(
    "/home/dirk/Projects/.scworkbench-worktrees/phase5-derived-fitted-axis-forge-20260922",
    "tests/fixtures/operation-pack-v3-phase5")
  stopifnot(dir.exists(fixture_source),
    requireNamespace("slingshot", quietly = TRUE),
    requireNamespace("fs", quietly = TRUE))
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
  old_model <- x$model
  old_model$pseudotime_table <- NULL
  old_model$model_obj <- .trajectory_s2_synthetic_slingshot(fitted)
  expected <- scworkbench::wb_op_trajectory_smooth(x$bundle, "RNA",
    old_model, cells = fitted, lineage = "L1",
    features = x$features[1:3], n_bins = 4L,
    min_cells_per_bin = 1L)$outputs
  for (name in names(expected)) stopifnot(identical(
    trajectory_identity(x$bundle, observed$outputs[[name]])$content_sha256,
    trajectory_identity(x$bundle, expected[[name]])$content_sha256))
  stopifnot(length(durable$model_tree_root) == 1L,
    length(prepared$values[[1L]]$manifest$primary_path) == 1L,
    length(durable$model_source_linkage$selected_axis_evidence$relative_path) == 1L)
  # Both attacks target the relocated retained model, before a private read.
  test_attack <- function(path) {
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    mode <- file.info(path)$mode
    Sys.chmod(path, mode = "0644")
    on.exit({
      connection <- file(path, "wb")
      writeBin(bytes, connection); close(connection)
      Sys.chmod(path, mode = mode)
    }, add = TRUE)
    connection <- file(path, "ab")
    writeBin(charToRaw("tamper"), connection); close(connection)
    before <- scworkbench::wb_art_list(x$bundle)
    rejected <- tryCatch(prepare(reopened, consume,
      list(model = durable), list(model = descriptor),
      tempfile("forge2-g1-reject-")), error = identity)
    stopifnot(inherits(rejected, "error"),
      identical(before, scworkbench::wb_art_list(x$bundle)))
  }
  test_attack(file.path(durable$model_tree_root,
    prepared$values[[1L]]$manifest$primary_path))
  test_attack(file.path(durable$model_tree_root,
    durable$model_source_linkage$selected_axis_evidence$relative_path))
  TRUE
}
