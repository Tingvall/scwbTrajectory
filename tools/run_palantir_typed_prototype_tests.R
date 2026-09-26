if (!nzchar(Sys.getenv("SCWORKBENCH_SOURCE", "")) ||
    !nzchar(Sys.getenv("SCWORKBENCH_PYBENCH_PYTHON", "")) ||
    !nzchar(Sys.getenv("PYTHONPATH", "")))
  stop("Set SCWORKBENCH_SOURCE, SCWORKBENCH_PYBENCH_PYTHON and PYTHONPATH.")
library(scworkbench)
library(scwbTrajectory)
source("tests/helper_palantir_typed_prototype.R")
source("tests/test_palantir_typed_prototype_contract.R")
tests <- sort(ls(.GlobalEnv,
  pattern = "^test_palantir_typed_prototype_.*_contract$"))
passed <- failed <- 0L
for (name in tests) {
  result <- tryCatch(get(name, .GlobalEnv)(), error = identity)
  if (inherits(result, "error")) {
    failed <- failed + 1L
    cat("FAIL", name, conditionMessage(result), "\n")
  } else {
    passed <- passed + 1L
    cat("PASS", name, "\n")
  }
}
cat(sprintf("RESULT discovered=%d executed=%d passed=%d skipped=0 failed=%d\n",
  length(tests), length(tests), passed, failed))
if (failed) quit(status = 1L)
