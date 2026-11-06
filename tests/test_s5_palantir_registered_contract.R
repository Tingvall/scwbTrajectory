test_dir <- if (file.exists("helper_trajectory_fixture.R")) "." else "tests"
source(file.path(test_dir, "helper_palantir_typed_prototype.R"),
  local = .GlobalEnv)
source(file.path(test_dir, "test_palantir_typed_prototype_contract.R"),
  local = .GlobalEnv)

.trajectory_palantir_registered <- function() {
  trajectory_activate_pack()
  list(fit = "scwbtrajectory::palantir_fit",
    project = "scwbtrajectory::palantir_project")
}

test_trajectory_palantir_registered_manifest_contract <- function() {
  proto <- .trajectory_palantir_registered()
  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbTrajectory"),
    simplifyVector = FALSE)
  names <- vapply(manifest$operations, `[[`, character(1),
    "qualified_name")
  stopifnot(length(names) == 31L,
    identical(sum(names == proto$fit), 1L),
    identical(sum(names == proto$project), 1L),
    !any(c("palantir_fit", "palantir_project") %in%
      scworkbench::wb_list_ops()),
    all(c(proto$fit, proto$project) %in% scworkbench::wb_list_ops()))
  fit <- manifest$operations[[match(proto$fit, names)]]
  project <- manifest$operations[[match(proto$project, names)]]
  inherited <- manifest
  inherited$operations <- Filter(function(operation)
    !as.character(operation$name) %in% c(
      "trajectory_pseudotime_from_cell_vector",
      "trajectory_signal_regions"), manifest$operations)
  inherited$operations <- inherited$operations[seq_len(23L)]
  inherited_fit <- match("slingshot_fit", vapply(inherited$operations,
    `[[`, character(1), "name"))
  inherited$operations[[inherited_fit]]$params_schema$properties$omega <- NULL
  inherited$operations[[inherited_fit]]$params_schema$properties$omega_scale <-
    NULL
  # Freeze the current inherited prefix. Its identity changes intentionally when
  # the earlier gene-set contracts gain depth-aware inputs and outputs; the two
  # Palantir operation definitions remain independently frozen below.
  inherited_identity <- getFromNamespace(".wb_operation_pack_v3_id",
    "scworkbench")("scworkbench.operation_pack.manifest.v3",
      inherited)$sha256
  definition_identity <- getFromNamespace(".wb_operation_pack_v3_id",
    "scworkbench")
  closure <- sort(c("pseudotime", "branch_membership",
    "branch_summary", "entropy", "model"))
  stopifnot(identical(inherited_identity,
      "ea4a033caa66ce4f301dbae60539297d9ac46691642891deea612a7ffcd4a8c1"),
    identical(definition_identity(
      "scworkbench.operation_pack.definition.v3", fit)$sha256,
      "d2da4a19f9ff497c8894231bb0eae7d4b9dfc6e4266590eb0a6e7af3a92d81d5"),
    identical(definition_identity(
      "scworkbench.operation_pack.definition.v3", project)$sha256,
      "64c41c2ca94c2475cecf0da93b5308d5fe18b6d353b190c1445dbff511234846"),
    identical(sort(vapply(fit$outputs, `[[`, character(1),
      "name")), closure),
    all(vapply(fit$promotion_closure, function(value)
      identical(as.character(unlist(value)), closure), logical(1))),
    identical(fit$outputs[[1L]]$derived_fitted_axis_policy$mode,
      "source_relative_subset_v1"),
    identical(fit$outputs[[1L]]$
      contained_model_sibling_closure_policy$primary_path,
      "model.pkl"),
    identical(fit$dimension_policies[[2L]]$axis_source$source,
      "g2_derived_fitted_axis_v1"),
    identical(project$self_linked_contained_model_reads[[1L]]$
      allowed_models[[1L]], list(engine = "palantir",
      model_type = "palantir_fit")),
    identical(project$params_schema$required,
      list("pickle_execution_trust")),
    identical(project$params_schema$properties$
      pickle_execution_trust$enum,
      list("operator_verified_source_v1")),
    identical(project$promotion_closure$curves, list("curves")),
    identical(project$replay_equivalence$curves,
      "exact_canonical_value"))
  TRUE
}

test_trajectory_palantir_registered_fit_project_oracle_contract <-
    function() {
  test_palantir_typed_prototype_fit_project_oracle_contract(
    .trajectory_palantir_registered())
}

test_trajectory_palantir_registered_missing_pack_predispatch_contract <-
    function() {
  proto <- .trajectory_palantir_registered()
  x <- palantir_proto_fixture()
  on.exit(palantir_proto_close(x), add = TRUE)
  registry <- getFromNamespace(".wb_operation_pack_registry",
    "scworkbench")
  saved_pack <- registry$packs$scwbtrajectory
  saved_operations <- registry$operations[saved_pack$operation_names]
  on.exit({
    registry$packs$scwbtrajectory <- saved_pack
    registry$operations[names(saved_operations)] <- saved_operations
  }, add = TRUE)
  for (name in saved_pack$operation_names)
    registry$operations[[name]] <- NULL
  registry$packs$scwbtrajectory <- NULL
  before <- palantir_proto_state(x)
  for (name in c(proto$fit, proto$project)) {
    rejected <- try(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), x$bundle, name,
      inputs = palantir_proto_fit_inputs()), silent = TRUE)
    stopifnot(inherits(rejected, "try-error"),
      identical(before, palantir_proto_state(x)))
  }
  stopifnot(!any(c("palantir_fit", "palantir_project") %in%
    scworkbench::wb_list_ops()))
  TRUE
}

test_trajectory_palantir_registered_pack_free_historical_inspection_contract <-
    function() {
  proto <- .trajectory_palantir_registered()
  x <- palantir_proto_fixture()
  on.exit(palantir_proto_close(x), add = TRUE)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, proto$fit, inputs = palantir_proto_fit_inputs())
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("pseudotime", "branch_membership",
      "branch_summary", "entropy", "model"))
  scworkbench::wb_session_discard(fit$session)
  model_id <- as.character(promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]])
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    x$bundle)$ready))
  scworkbench::wb_close_bundle(x$bundle)
  x$bundle <- scworkbench::wb_open_bundle(x$root, read_only = TRUE)
  registry <- getFromNamespace(".wb_operation_pack_registry",
    "scworkbench")
  saved_pack <- registry$packs$scwbtrajectory
  saved_operations <- registry$operations[saved_pack$operation_names]
  on.exit({
    registry$packs$scwbtrajectory <- saved_pack
    registry$operations[names(saved_operations)] <- saved_operations
  }, add = TRUE)
  for (name in saved_pack$operation_names)
    registry$operations[[name]] <- NULL
  registry$packs$scwbtrajectory <- NULL
  before <- palantir_proto_state(x)
  loaded <- scworkbench::wb_art_load(x$bundle, model_id,
    expected_kind = "model_ref")
  dag <- scworkbench::wb_prov_bundle_dag(x$bundle)
  report <- scworkbench::wb_bundle_validation_report(x$root,
    level = "recipient")
  if (!isTRUE(report$ok))
    stop(paste(capture.output(print(report$diagnostics)),
      collapse = "\n"), call. = FALSE)
  stopifnot(identical(loaded$engine, "palantir"),
    is.list(dag),
    identical(before, palantir_proto_state(x)))
  TRUE
}

test_trajectory_palantir_registered_staged_g1_contract <- function() {
  proto <- .trajectory_palantir_registered()
  x <- palantir_proto_fixture()
  on.exit(palantir_proto_close(x), add = TRUE)
  before <- palantir_proto_state(x)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, proto$fit, inputs = palantir_proto_fit_inputs())
  staged_model <- scworkbench::wb_ref_run(fit$run$run_id, "model")
  inputs <- palantir_proto_project_inputs(staged_model)
  project <- scworkbench::wb_session_run_op(fit$session, x$bundle,
    proto$project, inputs = inputs)
  curves <- project$run$outputs$curves
  stopifnot(length(curves$key) == 40L,
    identical(unique(as.character(curves$key)),
      c("branch:a09", "branch:b09")),
    identical(as.integer(curves$point), rep(seq_len(20L), 2L)),
    all(is.finite(curves$x)), all(is.finite(curves$y)),
    identical(before, palantir_proto_state(x)))
  scworkbench::wb_session_discard(project$session)
  TRUE
}

test_trajectory_palantir_registered_pickle_trust_contract <- function() {
  proto <- .trajectory_palantir_registered()
  x <- palantir_proto_fixture()
  on.exit(palantir_proto_close(x), add = TRUE)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, proto$fit, inputs = palantir_proto_fit_inputs())
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("pseudotime", "branch_membership",
      "branch_summary", "entropy", "model"))
  scworkbench::wb_session_discard(fit$session)
  id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  model <- scworkbench::wb_art_load(x$bundle, id,
    expected_kind = "model_ref")
  entry <- Filter(function(value) identical(value$path, "model.pkl"),
    model$model_tree_manifest$entries)
  bytes <- readBin(file.path(model$model_tree_root, "model.pkl"), "raw",
    n = file.info(file.path(model$model_tree_root, "model.pkl"))$size)
  stopifnot(length(entry) == 1L,
    identical(as.character(entry[[1L]]$sha256),
      digest::digest(bytes, algo = "sha256", serialize = FALSE)))
  # This hash-valid bundle is still untrusted by default. Integrity is not
  # permission to execute its pickle; only an explicit operator attestation
  # can cross the private-computation trust boundary.
  inputs <- palantir_proto_project_inputs(
    scworkbench::wb_ref_artifact(id, "model"))
  inputs$pickle_execution_trust <- NULL
  before <- palantir_proto_state(x)
  prior <- getOption("scwbTrajectory.palantir_typed_python")
  on.exit(options(scwbTrajectory.palantir_typed_python = prior), add = TRUE)
  false_program <- unname(Sys.which("false"))
  stopifnot(nzchar(false_program))
  options(scwbTrajectory.palantir_typed_python = false_program)
  denied <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, proto$project,
    inputs = inputs), silent = TRUE)
  stopifnot(inherits(denied, "try-error"),
    !grepl("Private Palantir computation failed", as.character(denied),
      fixed = TRUE), identical(before, palantir_proto_state(x)))
  inputs$pickle_execution_trust <- "untrusted_bundle"
  denied <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, proto$project,
    inputs = inputs), silent = TRUE)
  stopifnot(inherits(denied, "try-error"),
    !grepl("Private Palantir computation failed", as.character(denied),
      fixed = TRUE), identical(before, palantir_proto_state(x)))
  inputs$pickle_execution_trust <- "operator_verified_source_v1"
  reached <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, proto$project,
    inputs = inputs), silent = TRUE)
  if (!grepl("Private Palantir computation failed",
      as.character(reached), fixed = TRUE))
    stop("Unexpected trusted-model rejection: ",
      paste(as.character(reached), collapse = " "), call. = FALSE)
  stopifnot(inherits(reached, "try-error"),
    grepl("Private Palantir computation failed", as.character(reached),
      fixed = TRUE), identical(before, palantir_proto_state(x)))
  options(scwbTrajectory.palantir_typed_python = prior)
  retry <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, proto$project, inputs = inputs)
  stopifnot(length(retry$run$outputs$curves$key) == 40L)
  scworkbench::wb_session_discard(retry$session)
  TRUE
}

test_trajectory_palantir_registered_private_failure_contract <-
    function() {
  test_palantir_typed_prototype_private_failure_zero_residue_contract(
    .trajectory_palantir_registered())
}

test_trajectory_palantir_registered_source_model_tamper_contract <-
    function() {
  test_palantir_typed_prototype_source_model_tamper_contract(
    .trajectory_palantir_registered())
}

test_trajectory_palantir_registered_engine_linkage_contract <- function() {
  proto <- .trajectory_palantir_registered()
  x <- palantir_proto_fixture()
  on.exit(palantir_proto_close(x), add = TRUE)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, proto$fit, inputs = palantir_proto_fit_inputs())
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("pseudotime", "branch_membership",
      "branch_summary", "entropy", "model"))
  scworkbench::wb_session_discard(fit$session)
  id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  model <- scworkbench::wb_art_load(x$bundle, id,
    expected_kind = "model_ref")
  before <- palantir_proto_state(x)
  previous <- getOption("scwbTrajectory.palantir_typed_python")
  on.exit(options(scwbTrajectory.palantir_typed_python = previous),
    add = TRUE)
  options(scwbTrajectory.palantir_typed_python =
    unname(Sys.which("false")))
  for (field in c("engine", "model_type", "source_linkage")) {
    forged <- unserialize(serialize(model, NULL))
    if (identical(field, "source_linkage")) {
      forged$model_source_linkage$source_space_id <- "foreign_space"
    } else {
      forged[[field]] <- "untrusted_engine"
    }
    inputs <- palantir_proto_project_inputs(forged)
    rejected <- try(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), x$bundle, proto$project,
      inputs = inputs), silent = TRUE)
    stopifnot(inherits(rejected, "try-error"),
      !grepl("Private Palantir computation failed", as.character(rejected),
        fixed = TRUE), identical(before, palantir_proto_state(x)))
  }
  TRUE
}
