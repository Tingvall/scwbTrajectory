# Build a temporary test producer from the pack-owned, pinned generic G1
# fixture. Only the test package's model engine/type policy is specialized;
# S2 registration and the accepted core remain untouched.
.trajectory_s2_registered_producer_state <- new.env(parent = emptyenv())
.trajectory_s2_registered_producer_state$next_id <- 0L

trajectory_s2_g1_test_producer <- function(engine, model_type, model_factory) {
  stopifnot(engine %in% c("slingshot", "monocle3"),
    is.function(model_factory))
  root <- tempfile("s2-g1-producer-"); dir.create(root)
  source <- file.path(root, "source")
  library <- file.path(root, "library"); dir.create(library)
  fs::dir_copy(file.path(test_dir, "fixtures", "s2-g1-derived-axis"),
    source)
  .trajectory_s2_registered_producer_state$next_id <-
    .trajectory_s2_registered_producer_state$next_id + 1L
  suffix <- paste0(if (identical(engine, "slingshot"))
    "Slingshot" else "Monocle",
    .trajectory_s2_registered_producer_state$next_id)
  package_name <- paste0("scwbV3Phase5", suffix)
  namespace <- paste0("scwbv3phase5", tolower(suffix))
  for (relative in c("DESCRIPTION", "R/register.R",
      "inst/scworkbench/operation-pack.json")) {
    path <- file.path(source, relative)
    lines <- readLines(path, warn = FALSE)
    lines <- gsub("scwbV3Phase5", package_name, lines, fixed = TRUE)
    lines <- gsub("scwbv3phase5", namespace, lines, fixed = TRUE)
    writeLines(lines, path)
  }
  manifest <- file.path(source, "inst", "scworkbench",
    "operation-pack.json")
  value <- paste(readLines(manifest, warn = FALSE), collapse = "\n")
  old_policy <- paste0('"policy": ',
    '"staged_contained_model_sibling_closure_v1", "engine": "fixture"')
  new_policy <- sub('"engine": "fixture"',
    paste0('"engine": "', engine, '"'), old_policy, fixed = TRUE)
  stopifnot(grepl(old_policy, value, fixed = TRUE))
  value <- sub(old_policy, new_policy, value, fixed = TRUE)
  old_type <- '"model_type": "derived_axis_model", "primary_path": "model.bin"'
  new_type <- sub("derived_axis_model", model_type, old_type,
    fixed = TRUE)
  stopifnot(grepl(old_type, value, fixed = TRUE))
  value <- sub(old_type, new_type, value, fixed = TRUE)
  writeLines(value, manifest)
  log <- file.path(root, "install.log")
  status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL",
    paste0("--library=", shQuote(library)), shQuote(source)),
    stdout = log, stderr = log)
  if (!identical(status, 0L)) stop(paste(readLines(log,
    warn = FALSE), collapse = "\n"), call. = FALSE)
  previous <- .libPaths(); .libPaths(c(library, previous))
  loadNamespace(package_name, lib.loc = library)
  fit_runner <- function(readers, writers, params, context) {
    source <- readers$source()
    fitted <- as.character(params$fitted_cells)
    stopifnot(all(fitted %in% as.character(source$cell_id)))
    writers$declare_derived_fitted_axis(fitted)
    writers$write_terminal_primitive("embedding",
      scworkbench::wb_prim_cell_vector(stats::setNames(seq_along(fitted),
        fitted), value_type = "continuous"))
    writers$write_terminal_primitive("weights",
      scworkbench::wb_prim_cell_vector_family(
        rep(c("rising", "falling"), each = length(fitted)),
        rep(fitted, 2L),
        c(seq(0, 1, length.out = length(fitted)),
          seq(1, 0, length.out = length(fitted))),
        value_type = "continuous"))
    tree <- writers$open_model_tree("model")
    saveRDS(model_factory(fitted), tree$path(tree$primary_path))
    model_space_id <- if (model_type %in% c("slingshot", "monocle3_cds"))
      "fixture.derived" else paste0("fixture.derived.", model_type)
    model <- scworkbench::wb_prim_model_ref(model_space_id,
      as.character(context$input_ids$source)[1L], engine = engine)
    model$model_type <- model_type
    tree$finish(model)
    invisible(list(status = "ok"))
  }
  assignInNamespace(".fit_derived_axis_model", fit_runner,
    ns = package_name)
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = package_name, lib.loc = library)
  scworkbench::wb_operation_pack_activate(manifest, lib.loc = library)
  list(root = root, library = library, previous = previous,
    fit_op = paste0(namespace, "::fit_derived_axis_model"),
    cleanup = function() {
      try(scworkbench::wb_operation_pack_unload(namespace),
        silent = TRUE)
      if (package_name %in% loadedNamespaces())
        try(unloadNamespace(package_name), silent = TRUE)
      .libPaths(previous)
      unlink(root, recursive = TRUE, force = TRUE)
    })
}

trajectory_s2_promote_test_model <- function(x, cells = x$cells) {
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, producer$fit_op, inputs = list(
      source = "latent", dims = c(1L, 3L), fitted_cells = cells))
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("embedding", "model", "weights"),
    created_by = "s2-registration-contract")
  stopifnot(nrow(promoted$promoted) == 3L,
    all(promoted$promoted$status == "ok"))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  list(fit = fit, promoted = promoted, model_id = model_id)
}
