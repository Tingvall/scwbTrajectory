source(file.path(test_dir, "helper_s2_fixture.R"), local = .GlobalEnv)
source(file.path(test_dir, "test_s2_model_execution_contract.R"),
  local = .GlobalEnv)
source(file.path(test_dir, "helper_s2_registered_fixture.R"),
  local = .GlobalEnv)

.trajectory_s2_registered_explicit_cases <- function(x) {
  list(
    trajectory_geneset_dynamics_explicit = list(
      inputs = list(matrix = "RNA", pseudotime = x$explicit,
        membership = x$membership, namespace = "trajectory",
        pseudotime_key = "L1", min_set_size = 1L, n_bins = 5L,
        min_cells_per_bin = 1L, min_bins_observed = 2L,
        bootstrap_n = 8L, bootstrap_seed = 17L,
        n_permutations = 8L, permutation_seed = 19L,
        chunk_size = 2L)),
    trajectory_signal_dynamics_explicit = list(
      inputs = list(signal = x$signal, pseudotime = x$explicit,
        membership = x$membership, pseudotime_key = "L1",
        n_bins = 5L, min_cells_per_bin = 1L,
        min_bins_observed = 2L, bootstrap_n = 8L,
        bootstrap_seed = 17L)))
}

test_trajectory_s2_registered_manifest_contract <- function() {
  trajectory_activate_pack()
  names <- paste0("scwbtrajectory::", c(
    "trajectory_de", "trajectory_smooth",
    "trajectory_signature_smooth", "trajectory_geneset_dynamics",
    "trajectory_vector_smooth", "trajectory_signal_dynamics",
    "trajectory_geneset_dynamics_explicit",
    "trajectory_signal_dynamics_explicit"))
  stopifnot(all(names %in% scworkbench::wb_list_ops()))
  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbTrajectory"),
    simplifyVector = FALSE)
  declared <- Filter(function(item)
    as.character(item$qualified_name) %in% names,
    manifest$operations)
  stopifnot(length(declared) == 8L)
  for (item in declared) {
    outputs <- vapply(item$outputs, `[[`, character(1), "name")
    stopifnot(all(vapply(item$replay_equivalence,
      identical, logical(1), "exact_canonical_value")),
      identical(names(item$replay_equivalence), outputs),
      all(vapply(item$promotion_closure, function(closure)
        identical(as.character(unlist(closure)), sort(outputs)),
        logical(1))))
    if (grepl("_explicit$", item$name)) {
      stopifnot(is.null(item$self_linked_contained_model_reads),
        "pseudotime" %in% vapply(item$inputs, `[[`, character(1), "name"))
    } else {
      stopifnot(length(item$self_linked_contained_model_reads) == 1L,
        identical(item$self_linked_contained_model_reads[[1L]]$reader,
          "read_self_linked_contained_model_payload"))
      models <- item$self_linked_contained_model_reads[[1L]]$allowed_models
      pairs <- sort(vapply(models, function(model) paste(
        model$engine, model$model_type, sep = "/"), character(1)))
      stopifnot(identical(pairs,
        c("monocle3/monocle3_cds", "slingshot/slingshot")))
    }
  }
  TRUE
}

test_trajectory_s2_registered_explicit_publication_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  cases <- .trajectory_s2_registered_explicit_cases(x)
  reference <- trajectory_s2_reference()
  artifacts <- character()
  for (name in names(cases)) {
    result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
      x$bundle, paste0("scwbtrajectory::", name),
      inputs = cases[[name]]$inputs)
    trajectory_s2_compare_outputs(x$bundle, name,
      list(outputs = result$run$outputs), reference)
    closure <- names(result$run$outputs)
    promoted <- scworkbench::wb_session_promote_current(x$bundle,
      result$session, outputs = closure,
      created_by = "s2-registered-publication")
    details <- scworkbench::wb_prov_run_details(x$bundle,
      result$run$run_id, promoted$provenance_snapshot_id)
    stopifnot(setequal(promoted$promoted$output, closure),
      all(promoted$promoted$status == "ok"),
      all(details$outputs$verification_status == "content_verified"),
      all(details$outputs$replay_equivalence_policy ==
        "exact_canonical_value"))
    artifacts <- c(artifacts, promoted$promoted$artifact_id)
  }
  stopifnot(length(artifacts) == 10L,
    all(vapply(artifacts, function(id)
      inherits(scworkbench::wb_art_load(x$bundle, id), "wb_primitive"),
      logical(1))),
    isTRUE(scworkbench::wb_finalize_recipient_bundle(x$bundle)$ready))
  scworkbench::wb_close_bundle(x$bundle); x$bundle <- NULL
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
  reopened <- scworkbench::wb_open_bundle(x$root, read_only = TRUE)
  x$bundle <- reopened
  stopifnot(all(vapply(artifacts, function(id)
    inherits(scworkbench::wb_art_load(reopened, id), "wb_primitive"),
    logical(1))),
    isTRUE(scworkbench::wb_bundle_validation_report(x$root,
      "recipient")$ok))
  TRUE
}

test_trajectory_s2_registered_explicit_rollback_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit({
    options(scworkbench.promotion_workspace_publish_hook = NULL)
    trajectory_close(x)
  }, add = TRUE)
  for (name in names(.trajectory_s2_registered_explicit_cases(x))) {
    inputs <- .trajectory_s2_registered_explicit_cases(x)[[name]]$inputs
    result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
      x$bundle, paste0("scwbtrajectory::", name), inputs = inputs)
    before_artifacts <- scworkbench::wb_art_list(x$bundle)
    before_state <- getFromNamespace(".wb_session_source_identity",
      "scworkbench")(x$bundle)$scientific_state_fingerprint
    options(scworkbench.promotion_workspace_publish_hook = function(...)
      stop("injected S2 publication failure"))
    rejected <- try(scworkbench::wb_session_promote_current(x$bundle,
      result$session, outputs = names(result$run$outputs)), silent = TRUE)
    options(scworkbench.promotion_workspace_publish_hook = NULL)
    stopifnot(inherits(rejected, "try-error"),
      identical(before_artifacts, scworkbench::wb_art_list(x$bundle)),
      identical(before_state, getFromNamespace(
        ".wb_session_source_identity", "scworkbench")(
          x$bundle)$scientific_state_fingerprint))
    retry <- scworkbench::wb_session_promote_current(x$bundle,
      result$session, outputs = names(result$run$outputs))
    stopifnot(setequal(retry$promoted$output,
      names(result$run$outputs)),
      all(retry$promoted$status == "ok"))
  }
  TRUE
}

test_trajectory_s2_registered_missing_pack_dispatch_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  registry <- getFromNamespace(".wb_operation_pack_registry",
    "scworkbench")
  saved_pack <- registry$packs$scwbtrajectory
  saved_operations <- registry$operations[saved_pack$operation_names]
  on.exit({
    registry$packs$scwbtrajectory <- saved_pack
    registry$operations[names(saved_operations)] <- saved_operations
  }, add = TRUE)
  before <- scworkbench::wb_art_list(x$bundle)
  for (name in saved_pack$operation_names)
    registry$operations[[name]] <- NULL
  registry$packs$scwbtrajectory <- NULL
  rejected <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_signal_dynamics_explicit",
    inputs = .trajectory_s2_registered_explicit_cases(x)$
      trajectory_signal_dynamics_explicit$inputs), silent = TRUE)
  stopifnot(inherits(rejected, "try-error"),
    identical(before, scworkbench::wb_art_list(x$bundle)))
  TRUE
}

test_trajectory_s2_registered_explicit_rejection_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  before <- scworkbench::wb_art_list(x$bundle)
  for (name in names(.trajectory_s2_registered_explicit_cases(x))) {
    inputs <- .trajectory_s2_registered_explicit_cases(x)[[name]]$inputs
    inputs$pseudotime_key <- "missing"
    bad_key <- try(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), x$bundle,
      paste0("scwbtrajectory::", name), inputs = inputs), silent = TRUE)
    inputs <- .trajectory_s2_registered_explicit_cases(x)[[name]]$inputs
    inputs$membership$cell_id[[1L]] <- "foreign_cell"
    bad_membership <- try(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), x$bundle,
      paste0("scwbtrajectory::", name), inputs = inputs), silent = TRUE)
    stopifnot(inherits(bad_key, "try-error"),
      inherits(bad_membership, "try-error"),
      identical(before, scworkbench::wb_art_list(x$bundle)))
  }
  TRUE
}

test_trajectory_s2_registered_model_g1_core_gap_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  producer <- trajectory_s2_g1_test_producer("slingshot", "slingshot",
    .trajectory_s2_synthetic_slingshot)
  on.exit(producer$cleanup(), add = TRUE)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbv3phase5::fit_derived_axis_model", inputs = list(
      source = "latent", dims = c(1L, 3L),
      fitted_cells = x$cells))
  before_staged <- scworkbench::wb_art_list(x$bundle)
  staged <- tryCatch(scworkbench::wb_session_run_op(fit$session,
    x$bundle, "scwbtrajectory::trajectory_de", inputs = list(
      matrix = "RNA", model = scworkbench::wb_ref_run(
        fit$run$run_id, "model"), lineage = "L1",
      features = x$features[1:3], chunk_size = 2L)), error = identity)
  stopifnot(inherits(staged, "error"),
    grepl("missing required property source_space_content_identity",
      conditionMessage(staged), fixed = TRUE),
    identical(before_staged, scworkbench::wb_art_list(x$bundle)))
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("embedding", "model", "weights"),
    created_by = "s2-g1-gap-reproducer")
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  row <- DBI::dbGetQuery(x$bundle$con, paste0(
    "SELECT model_path, model_primary_path FROM artifact_model_ref ",
    "WHERE artifact_id = ?"), params = list(model_id))
  stopifnot(nrow(row) == 1L)
  tree <- file.path(x$root, row$model_path[[1L]])
  stopifnot(dir.exists(tree), file.exists(file.path(tree,
    row$model_primary_path[[1L]])))
  load_error <- tryCatch(scworkbench::wb_art_load(x$bundle,
    model_id, expected_kind = "model_ref"), error = identity)
  before_durable <- scworkbench::wb_art_list(x$bundle)
  durable <- tryCatch(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_de", inputs = list(
      matrix = "RNA", model = scworkbench::wb_ref_artifact(
        model_id, "model"), lineage = "L1",
      features = x$features[1:3], chunk_size = 2L)), error = identity)
  stopifnot(inherits(load_error, "error"),
    inherits(durable, "error"),
    grepl("cannot open the connection", conditionMessage(load_error),
      fixed = TRUE),
    grepl("cannot open the connection", conditionMessage(durable),
      fixed = TRUE),
    identical(before_durable, scworkbench::wb_art_list(x$bundle)))
  TRUE
}
