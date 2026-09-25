test_trajectory_s2_ownership_catalog_contract <- function() {
  script <- tempfile("trajectory-s2-catalog-", fileext = ".R")
  log <- tempfile("trajectory-s2-catalog-", fileext = ".log")
  on.exit(unlink(c(script, log)), add = TRUE)
  lib_paths <- paste(capture.output(dput(.libPaths())), collapse = "")
  writeLines(c(
    paste0(".libPaths(", lib_paths, ")"),
    "library(scworkbench)",
    "legacy <- c('trajectory_de', 'trajectory_smooth',",
    "  'trajectory_signature_smooth', 'trajectory_geneset_dynamics',",
    "  'trajectory_vector_smooth', 'trajectory_signal_dynamics')",
    "explicit <- c('trajectory_geneset_dynamics_explicit',",
    "  'trajectory_signal_dynamics_explicit')",
    "owned <- c(legacy, explicit)",
    "before <- wb_list_ops()",
    "stopifnot(all(vapply(owned, function(name)",
    "  sum(before == name) == 0L &&",
    "  sum(before == paste0('scwbtrajectory::', name)) == 0L,",
    "  logical(1))))",
    "before_identity <- wb_operation_catalog_identity()$sha256",
    "manifest <- system.file('scworkbench', 'operation-pack.json',",
    "  package = 'scwbTrajectory')",
    "stopifnot(nzchar(manifest))",
    "wb_operation_pack_activate(manifest, lib.loc = .libPaths())",
    "after <- wb_list_ops()",
    "stopifnot(all(vapply(owned, function(name)",
    "  sum(after == name) == 0L &&",
    "  sum(after == paste0('scwbtrajectory::', name)) == 1L,",
    "  logical(1))),",
    "  !identical(before_identity, wb_operation_catalog_identity()$sha256))",
    "cat('RESULT no-pack=0 pack-qualified=1 names=8\\n')"
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    args = shQuote(script), stdout = log, stderr = log)
  if (!identical(as.integer(status), 0L))
    stop(paste(readLines(log, warn = FALSE), collapse = "\n"),
      call. = FALSE)
  stopifnot(any(readLines(log, warn = FALSE) ==
    "RESULT no-pack=0 pack-qualified=1 names=8"))
  TRUE
}
