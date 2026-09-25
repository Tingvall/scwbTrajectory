# Execute in separate R installations: old core (`emit old PATH`), current
# core plus activated scwbTrajectory (`emit pack PATH`), then `compare A B`.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(3L, 4L), args[[1L]] %in% c("emit", "compare"))

family_rows <- function(x) data.frame(key = as.character(x$key),
  cell_id = as.character(x$cell_id), value = as.numeric(x$value))
path_rows <- function(x) data.frame(key = as.character(x$key),
  point = as.integer(x$point), x = as.numeric(x$x), y = as.numeric(x$y))

if (identical(args[[1L]], "emit")) {
  mode <- args[[2L]]
  stopifnot(mode %in% c("old", "pack"))
  library(scworkbench)
  if (identical(mode, "pack")) {
    manifest <- system.file("scworkbench", "operation-pack.json",
      package = "scwbTrajectory")
    stopifnot(nzchar(manifest))
    wb_operation_pack_activate(manifest, lib.loc = .libPaths())
  }
  prefix <- if (identical(mode, "pack")) "scwbtrajectory::" else ""
  root <- tempfile("s3-legacy-comparator-")
  b <- wb_init_bundle(root)
  on.exit({try(wb_close_bundle(b), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
  cells <- paste0("cell_", seq_len(120L))
  groups <- rep(c("root", "left", "right"), each = 40L)
  coordinates <- cbind(
    PC1 = rep(c(0, -2, 2), each = 40L) +
      rep(seq(-.2, .2, length.out = 40L), 3L),
    PC2 = rep(c(0, 2, 2), each = 40L) +
      sin(seq_len(120L) * .37) * .12,
    PC3 = rep(c(0, -1, 1), each = 40L) +
      cos(seq_len(120L) * .23) * .11)
  rownames(coordinates) <- cells
  wb_define_cells(b, cells)
  wb_add_space(b, "fit", coordinates)
  plot <- coordinates[, 1:2, drop = FALSE]
  colnames(plot) <- c("x", "y")
  wb_add_space(b, "plot", plot)
  clusters <- wb_prim_cell_vector(stats::setNames(groups, cells),
    value_type = "categorical")
  shared <- list(clusters = clusters, start.clus = "root",
    end.clus = c("left", "right"), min_cells = 1L)
  legacy <- wb_session_run_op(wb_session_new(), b,
    paste0(prefix, "slingshot"),
    inputs = c(list(space_2d = "plot"), shared))
  fit <- wb_session_run_op(wb_session_new(), b,
    paste0(prefix, "slingshot_fit"),
    inputs = c(list(fit_space = "fit"), shared))
  model <- wb_ref_run(fit$run$run_id, "model")
  weights <- wb_session_run_op(fit$session, b,
    paste0(prefix, "slingshot_lineage_weights"),
    inputs = list(model = model, min_cells = 1L))
  project <- wb_session_run_op(fit$session, b,
    paste0(prefix, "slingshot_project"),
    inputs = list(model = model, plot_space = "plot",
      approx_points = 20L, min_cells = 1L))
  result <- list(
    terminal_names = list(legacy = sort(names(legacy$run$outputs)),
      fit = sort(names(fit$run$outputs)),
      weights = sort(names(weights$run$outputs)),
      project = sort(names(project$run$outputs))),
    legacy_pseudotime = family_rows(legacy$run$outputs$pseudotime),
    legacy_curves = path_rows(legacy$run$outputs$curves),
    fit_pseudotime = family_rows(fit$run$outputs$pseudotime),
    fit_weights = family_rows(fit$run$outputs$lineage_weights),
    extracted_weights = family_rows(weights$run$outputs$lineage_weights),
    projected_curves = path_rows(project$run$outputs$curves),
    fit_cells = as.character(fit$run$outputs$model$cells),
    fit_engine = as.character(fit$run$outputs$model$engine),
    fit_lineages = as.integer(fit$run$outputs$model$n_lineages),
    effective_parameters = list(
      start.clus = fit$run$params[["start.clus"]],
      end.clus = fit$run$params[["end.clus"]],
      stretch = fit$run$params$stretch,
      reweight = fit$run$params$reweight,
      shrink = fit$run$params$shrink,
      maxit = fit$run$params$maxit,
      approx_points = project$run$params$approx_points))
  stopifnot(identical(as.character(result$effective_parameters$start.clus),
      "root"),
    identical(as.character(result$effective_parameters$end.clus),
      c("left", "right")),
    isTRUE(all.equal(as.numeric(result$effective_parameters$stretch), 2)),
    identical(as.logical(result$effective_parameters$reweight), TRUE),
    isTRUE(all.equal(as.numeric(result$effective_parameters$shrink), 1)),
    identical(as.integer(result$effective_parameters$maxit), 10L),
    identical(as.integer(result$effective_parameters$approx_points), 20L))
  saveRDS(result, args[[3L]], version = 3L)
  cat("EMIT", mode, "terminals=4 records=",
    sum(vapply(result[c("legacy_pseudotime", "legacy_curves",
      "fit_pseudotime", "fit_weights", "extracted_weights",
      "projected_curves")], nrow, integer(1))), "\n")
} else {
  old <- readRDS(args[[2L]])
  pack <- readRDS(args[[3L]])
  stopifnot(identical(old$terminal_names, pack$terminal_names),
    identical(old$fit_cells, pack$fit_cells),
    identical(old$fit_engine, pack$fit_engine),
    identical(old$fit_lineages, pack$fit_lineages),
    identical(old$effective_parameters, pack$effective_parameters))
  for (name in c("legacy_pseudotime", "fit_pseudotime", "fit_weights",
    "extracted_weights")) {
    a <- old[[name]]; b <- pack[[name]]
    stopifnot(identical(a$key, b$key), identical(a$cell_id, b$cell_id),
      isTRUE(all.equal(a$value, b$value, tolerance = 1e-8)))
  }
  for (name in c("legacy_curves", "projected_curves")) {
    a <- old[[name]]; b <- pack[[name]]
    stopifnot(identical(a$key, b$key), identical(a$point, b$point),
      isTRUE(all.equal(a$x, b$x, tolerance = 1e-8)),
      isTRUE(all.equal(a$y, b$y, tolerance = 1e-8)))
  }
  cat("PASS S3 old-core/new-pack terminal closure, ordered axes, values,",
    "effective parameters (numeric tolerance 1e-8)\n")
}
