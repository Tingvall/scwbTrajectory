test_dir <- if (file.exists("helper_trajectory_fixture.R")) "." else "tests"
source(file.path(test_dir, "helper_trajectory_fixture.R"),
  local = .GlobalEnv)
files <- sort(list.files(test_dir, pattern = "^test_.*\\.R$",
  full.names = TRUE))
for (file in files) source(file, local = .GlobalEnv)
tests <- sort(ls(.GlobalEnv, pattern = "^test_trajectory_"))
passed <- failed <- skipped <- 0L
for (name in tests) {
  value <- tryCatch(get(name, .GlobalEnv)(), error = identity)
  if (inherits(value, "error")) {
    failed <- failed + 1L
    cat("FAIL", name, conditionMessage(value), "\n")
  } else if (!is.null(attr(value, "skip"))) {
    skipped <- skipped + 1L
    cat("SKIP", name, attr(value, "skip"), "\n")
  } else {
    passed <- passed + 1L
    cat("PASS", name, "\n")
  }
}
cat(sprintf("RESULT discovered=%d executed=%d passed=%d skipped=%d failed=%d\n",
  length(tests), length(tests), passed, skipped, failed))
if (failed) quit(status = 1L)
