if (!nzchar(Sys.getenv("SCWORKBENCH_SOURCE", "")) ||
    !nzchar(Sys.getenv("SCWBTRAJECTORY_SOURCE", "")) ||
    !nzchar(Sys.getenv("SCWORKBENCH_PYBENCH_PYTHON", "")) ||
    !nzchar(Sys.getenv("PYTHONPATH", "")))
  stop("Set exact core, pack, Python interpreter and private PyBench source pins.")

library(scworkbench)
library(scwbTrajectory)
source("tests/helper_trajectory_fixture.R")
source("tests/test_s5_monocle3_contract.R")
source("tests/test_s5_palantir_registered_contract.R")

manifest <- jsonlite::read_json(system.file("scworkbench",
  "operation-pack.json", package = "scwbTrajectory"),
  simplifyVector = FALSE)
names <- vapply(manifest$operations, `[[`, character(1), "qualified_name")
stopifnot(length(names) == 25L,
  identical(tail(names, 4L), c("scwbtrajectory::monocle3_fit",
    "scwbtrajectory::monocle3_project", "scwbtrajectory::palantir_fit",
    "scwbtrajectory::palantir_project")))

tests <- c("test_trajectory_s5_monocle3_registered_oracle_reopen_contract",
  "test_trajectory_palantir_registered_fit_project_oracle_contract")
for (name in tests) {
  get(name, .GlobalEnv)()
  cat("PASS", name, "\n")
}
cat("S5 INSTALLED PAIR COMPOSITE OK discovered=2 executed=2 passed=2 skipped=0 failed=0\n")
