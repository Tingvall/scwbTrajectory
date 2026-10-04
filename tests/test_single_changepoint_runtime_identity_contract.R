# Runtime-identity contract tests for
# scwbtrajectory::trajectory_single_changepoint_regions.
#
# The accepted ruling is "probe behaviour plus recorded identity, no exact-version
# gate". These tests attack that claim from the direction the other tests cannot
# reach: they shadow the ruptures *module* on PYTHONPATH with a stub, so a run
# happens under a library that is demonstrably not the one the distribution
# metadata advertises. Four things are pinned. An implausible version is recorded
# and never rejected. The numbers that reach the terminals come from whatever
# library was imported. A dependency that cannot serve the API the core needs, or
# answers with something that is not a boundary, stops the run with the documented
# code instead of publishing anything. And the probe is a correctness gate, not a
# liveness check: a library that fails to recover the planted step is refused,
# which is the property that makes the recorded identity worth having.
#
# The stub shadows the module only. importlib.metadata still reports the installed
# distribution, so ruptures_version and ruptures_dunder_version disagree during
# these runs — that disagreement is precisely why both are published.
sti_params <- list(smooth_window = 3L, window_width = 8L, window_jump = 1L,
  window_min_size = 2L, direction_window = 3L, density_smooth_window = 1L,
  density_scale = 1, change_threshold = 2, chunk_size = 256L,
  degenerate_bkp_policy = "drop", non_finite_policy = "drop")

# The only real change point is the step at rank 16, so a boundary reported
# anywhere else is visibly the library's answer and not the data's.
sti_fixture <- function(n = 30L) {
  trajectory_activate_pack()
  bundle <- wb_init_bundle(tempfile("scwbTrajectory-sti-"), overwrite = TRUE)
  cells <- sprintf("c%02d", seq_len(n))
  values <- matrix(0, nrow = 1L, ncol = n)
  values[seq.int(16L, n)] <- 9
  dimnames(values) <- list("step", cells)
  wb_define_cells(bundle, cells)
  wb_attach_matrix(bundle, "RNA", Matrix::Matrix(values, sparse = TRUE),
    feature_type = "gene")
  list(bundle = bundle, cells = cells,
    pseudotime = wb_prim_cell_vector(setNames((seq_len(n) - 1L) / (n - 1L),
      cells), "continuous"))
}

sti_run <- function(x, params = list()) {
  outputs <- wb_session_run_op(wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_single_changepoint_regions",
    inputs = c(list(matrix = "RNA", pseudotime = x$pseudotime),
      modifyList(sti_params, params)))$run$outputs
  list(features = outputs$feature_single_changepoints$table,
    summary = outputs$run_summary$table)
}

sti_error <- function(x, params = list()) {
  caught <- tryCatch(sti_run(x, params), error = function(error)
    conditionMessage(error))
  if (!is.character(caught))
    stop("the stubbed run was expected to fail and did not.")
  paste(caught, collapse = "\n")
}

sti_with_stub <- function(source, expr) {
  dir <- tempfile("scwbTrajectory-stub-")
  if (!dir.create(dir, recursive = TRUE, mode = "0700"))
    stop("could not stage the stub ruptures directory.")
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(source, file.path(dir, "ruptures.py"))
  previous <- Sys.getenv("PYTHONPATH", unset = NA_character_)
  Sys.setenv(PYTHONPATH = dir)
  on.exit(if (is.na(previous)) Sys.unsetenv("PYTHONPATH") else
    Sys.setenv(PYTHONPATH = previous), add = TRUE)
  force(expr)
}

# The stub answers the probe (a 200-point signal with a step planted at 100) with
# the boundary the probe demands, so the run reaches the features phase; every
# other signal gets `on_data`, which is whatever boundary the test is about.
sti_stub <- function(version, on_probe, on_data) paste0(
  '"""Stub used only by scwbTrajectory runtime-identity contract tests."""\n',
  '__version__ = "', version, '"\n\n\n',
  'class Window:\n',
  '    def __init__(self, **kw):\n',
  '        pass\n\n',
  '    def fit(self, signal, **kw):\n',
  '        self.signal = signal\n',
  '        return self\n\n',
  '    def predict(self, n_bkps=1, **kw):\n',
  '        if len(self.signal) == 200:\n',
  '            return ', on_probe, '\n',
  '        return ', on_data, '\n')

sti_broken_import <- '"""Stub that fails the way a broken install does."""
raise ImportError("stub: ruptures is present but will not import")
'

test_trajectory_single_changepoint_no_version_gate_contract <- function() {
  x <- sti_fixture()
  on.exit(trajectory_close(x), add = TRUE)
  ran <- sti_with_stub(sti_stub("0.0.1-dev-not-a-real-release", "[100, 200]",
    "[12, len(self.signal)]"), sti_run(x))
  summary <- ran$summary
  stopifnot(
    # No version gate exists: a version no ruptures release ever carried runs to
    # completion and is published verbatim.
    identical(summary$ruptures_dunder_version, "0.0.1-dev-not-a-real-release"),
    # The distribution metadata is untouched by shadowing the module, so the two
    # published identity fields disagree and nothing reconciles them.
    is.character(summary$ruptures_version), nzchar(summary$ruptures_version),
    !identical(summary$ruptures_version, summary$ruptures_dunder_version),
    # The break that reached the terminal is the stub's, not the step's: the
    # operation trusts the library it imported and says which one it used.
    identical(nrow(ran$features), 1L),
    identical(ran$features$analysed, TRUE),
    identical(ran$features$boundary_rank, 12L),
    identical(ran$features$boundary_before_cell_id, "c12"),
    identical(ran$features$boundary_after_cell_id, "c13"))
}

test_trajectory_single_changepoint_probe_missed_step_contract <- function() {
  x <- sti_fixture()
  on.exit(trajectory_close(x), add = TRUE)
  before <- wb_art_list(x$bundle)
  # Structurally perfect — two integers, interior, terminal count correct — and
  # still refused, because the probe planted a step at 100 and this answer is not
  # it. The probe is a correctness gate, so the recorded identity means the run
  # finished under a dependency that demonstrably found the thing it planted.
  message <- sti_with_stub(sti_stub("1.1.10", "[7, 200]", "[12, 30]"),
    sti_error(x))
  stopifnot(grepl("RUNTIME_API_UNAVAILABLE", message, fixed = TRUE),
    grepl("probe recovered boundary 7 for the planted step at 100", message,
      fixed = TRUE),
    identical(before, wb_art_list(x$bundle)))
}

test_trajectory_single_changepoint_probe_sample_count_contract <- function() {
  x <- sti_fixture()
  on.exit(trajectory_close(x), add = TRUE)
  before <- wb_art_list(x$bundle)
  message <- sti_with_stub(sti_stub("1.1.10", "[100, 100]", "[12, 30]"),
    sti_error(x))
  stopifnot(grepl("RUNTIME_API_UNAVAILABLE", message, fixed = TRUE),
    grepl("probe terminal sample count is 100, expected 200", message,
      fixed = TRUE),
    identical(before, wb_art_list(x$bundle)))
}

test_trajectory_single_changepoint_window_api_rejection_contract <- function() {
  x <- sti_fixture()
  on.exit(trajectory_close(x), add = TRUE)
  message <- sti_with_stub('"""Stub used only by scwbTrajectory contract tests."""
__version__ = "1.1.10"


class Window:
    def __init__(self, model=None, width=None, jump=None, min_size=None, **kw):
        raise TypeError("no cost kernel named " + str(model))

    def fit(self, signal, **kw):
        return self

    def predict(self, n_bkps=1, **kw):
        return [100, 200]
', sti_error(x))
  stopifnot(grepl("RUNTIME_API_UNAVAILABLE", message, fixed = TRUE),
    grepl("could not be constructed", message, fixed = TRUE),
    grepl("model='l1'", message, fixed = TRUE))
}

test_trajectory_single_changepoint_malformed_predict_contract <- function() {
  x <- sti_fixture()
  on.exit(trajectory_close(x), add = TRUE)
  message <- sti_with_stub('"""Stub used only by scwbTrajectory contract tests."""
__version__ = "1.1.10"


class Window:
    def __init__(self, model=None, width=None, jump=None, min_size=None, **kw):
        pass

    def fit(self, signal, **kw):
        self.signal = signal
        return self

    def predict(self, n_bkps=1, **kw):
        return ["not", "boundaries"]
', sti_error(x))
  stopifnot(grepl("RUNTIME_API_UNAVAILABLE", message, fixed = TRUE),
    grepl("returned non-integer", message, fixed = TRUE))
}

test_trajectory_single_changepoint_import_failure_contract <- function() {
  x <- sti_fixture()
  on.exit(trajectory_close(x), add = TRUE)
  before <- wb_art_list(x$bundle)
  # A distribution that is installed but will not import is a different failure
  # from a missing interpreter, and the identity probe survives it (it only reads
  # metadata), so the stop has to come from the run itself.
  message <- sti_with_stub(sti_broken_import, sti_error(x))
  stopifnot(grepl("RUNTIME_IMPORT_UNAVAILABLE", message, fixed = TRUE),
    grepl("ruptures is not importable", message, fixed = TRUE),
    identical(before, wb_art_list(x$bundle)))
}

test_trajectory_single_changepoint_out_of_range_predict_contract <- function() {
  x <- sti_fixture()
  on.exit(trajectory_close(x), add = TRUE)
  stub <- sti_stub("1.1.10", "[100, 200]", "[999, 1000]")
  dropped <- sti_with_stub(stub, sti_run(x))
  stopifnot(identical(nrow(dropped$features), 1L),
    identical(dropped$features$analysed, FALSE),
    identical(dropped$features$reason, "degenerate"),
    identical(dropped$features$boundary_rank, NA_integer_))
  # The same answer under degenerate_bkp_policy = "error" stops the run instead.
  message <- sti_with_stub(stub, sti_error(x,
    list(degenerate_bkp_policy = "error")))
  stopifnot(grepl("no interior breakpoint", message, fixed = TRUE))
}
