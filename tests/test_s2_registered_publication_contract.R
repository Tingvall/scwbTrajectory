test_dir <- if (file.exists("helper_s2_fixture.R")) "." else "tests"
source(file.path(test_dir, "helper_s2_fixture.R"), local = .GlobalEnv)
source(file.path(test_dir, "test_s2_model_execution_contract.R"),
  local = .GlobalEnv)
source(file.path(test_dir, "helper_s2_registered_fixture.R"),
  local = .GlobalEnv)
source(file.path(test_dir, "test_s2_monocle3_execution_contract.R"),
  local = .GlobalEnv)

.trajectory_s2_registered_explicit_cases <- function(x) {
  list(
    trajectory_geneset_dynamics_explicit = list(
      inputs = list(matrix = "RNA", pseudotime = x$explicit,
        membership = x$membership, namespace = "trajectory",
        pseudotime_key = "L1", feature_sets = "selected",
        min_set_size = 1L, n_bins = 5L,
        adjust_depth = FALSE,
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
  stopifnot(length(artifacts) == 12L,
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

.trajectory_s2_registered_model_cases <- function(x, engine, model_ref,
                                                   lineage = NULL,
                                                   signal = x$signal) {
  if (is.null(lineage)) lineage <- if (identical(engine, "slingshot"))
    "L1" else "root:A::lineage:C"
  operations <- c("de", "smooth", "signature_smooth",
    "geneset_dynamics", "vector_smooth", "signal_dynamics")
  stats::setNames(lapply(operations, function(operation) {
    params <- .trajectory_s2_model_params(operation)
    params$lineage <- lineage
    input <- switch(operation,
      de = list(matrix = "RNA"), smooth = list(matrix = "RNA"),
      signature_smooth = list(matrix = "RNA"),
      geneset_dynamics = list(matrix = "RNA"),
      vector_smooth = list(signal = signal),
      signal_dynamics = list(signal = signal))
    c(input, list(model = model_ref), params)
  }), paste0("trajectory_", operations))
}

.trajectory_s2_registered_model_oracle_key <- function(operation, engine,
                                                        lineage) {
  if (identical(engine, "slingshot")) {
    if (operation %in% c("trajectory_geneset_dynamics",
      "trajectory_signal_dynamics")) return(paste0(operation, "_model"))
    return(operation)
  }
  paste(sub("^trajectory_", "", operation),
    if (identical(lineage, "root:A::lineage:C")) "C" else "D",
    "full", sep = "__")
}

.trajectory_s2_registered_model_reference <- function(engine) {
  if (identical(engine, "slingshot")) trajectory_s2_reference() else
    jsonlite::read_json(file.path(test_dir, "fixtures",
      "trajectory_s2_monocle3_reference.json"), simplifyVector = FALSE)
}

.trajectory_s2_registered_model_setup <- function(x, engine,
                                                   cells = x$cells) {
  factory <- if (identical(engine, "slingshot"))
    .trajectory_s2_synthetic_slingshot else .trajectory_s2_branching_cds
  model_type <- if (identical(engine, "slingshot")) "slingshot" else
    "monocle3_cds"
  producer <- trajectory_s2_g1_test_producer(engine, model_type, factory)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, producer$fit_op, inputs = list(
      source = "latent", dims = c(1L, 3L), fitted_cells = cells))
  list(producer = producer, fit = fit)
}

.trajectory_s2_registered_durable_snapshot <- function(x) {
  con <- x$bundle$con
  tables <- sort(DBI::dbListTables(con))
  counts <- stats::setNames(vapply(tables, function(table) {
    query <- paste0("SELECT COUNT(*) AS n FROM ",
      as.character(DBI::dbQuoteIdentifier(con, table)))
    as.integer(DBI::dbGetQuery(con, query)$n[[1L]])
  }, integer(1)), tables)
  files <- list.files(x$root, recursive = TRUE, all.files = TRUE,
    include.dirs = FALSE)
  db_file <- basename(DBI::dbGetInfo(con)$dbname)
  files <- sort(files[!basename(files) %in% c("manifest.json", db_file,
    paste0(db_file, c("-shm", "-wal", "-journal")))])
  payload <- stats::setNames(vapply(files, function(relative)
    digest::digest(file = file.path(x$root, relative), algo = "sha256"),
    character(1)), files)
  list(artifacts = scworkbench::wb_art_list(x$bundle),
    state = getFromNamespace(".wb_session_source_identity",
      "scworkbench")(x$bundle)$scientific_state_fingerprint,
    counts = counts, payload = payload)
}

.trajectory_s2_registered_without_pack <- function(code) {
  registry <- getFromNamespace(".wb_operation_pack_registry",
    "scworkbench")
  saved <- registry$packs$scwbtrajectory
  operations <- registry$operations[saved$operation_names]
  on.exit({
    registry$packs$scwbtrajectory <- saved
    registry$operations[names(operations)] <- operations
  }, add = TRUE)
  for (name in saved$operation_names) registry$operations[[name]] <- NULL
  registry$packs$scwbtrajectory <- NULL
  force(code)
}

test_trajectory_s2_registered_de_full_axis_stream_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  setup <- .trajectory_s2_registered_model_setup(x, "slingshot")
  on.exit(setup$producer$cleanup(), add = TRUE)
  reference <- jsonlite::read_json(file.path(test_dir, "fixtures",
    "trajectory_s2_de_streaming_reference.json"), simplifyVector = FALSE)
  model_ref <- scworkbench::wb_ref_run(setup$fit$run$run_id, "model")
  before <- .trajectory_s2_registered_durable_snapshot(x)
  for (size in c(2L, 3L)) {
    full <- scworkbench::wb_session_run_op(setup$fit$session,
      x$bundle, "scwbtrajectory::trajectory_de", inputs = list(
        matrix = "RNA", model = model_ref, chunk_size = size))
    trajectory_s2_compare_outputs(x$bundle,
      "all_lineages_all_features", list(outputs = full$run$outputs),
      reference)
    selected <- scworkbench::wb_session_run_op(setup$fit$session,
      x$bundle, "scwbtrajectory::trajectory_de", inputs = list(
        matrix = "RNA", model = model_ref, chunk_size = size,
        features = c("g5", "g1", "g3"), lineage = "L1"))
    trajectory_s2_compare_outputs(x$bundle, "selected_reordered",
      list(outputs = selected$run$outputs), reference)
    stopifnot(identical(before,
      .trajectory_s2_registered_durable_snapshot(x)))
  }
  TRUE
}

test_trajectory_s2_registered_model_staged_oracles_contract <-
    function() {
  for (engine in c("slingshot", "monocle3")) {
    x <- trajectory_s2_fixture()
    setup <- .trajectory_s2_registered_model_setup(x, engine)
    reference <- .trajectory_s2_registered_model_reference(engine)
    staged_ref <- scworkbench::wb_ref_run(setup$fit$run$run_id, "model")
    staged_cases <- .trajectory_s2_registered_model_cases(x, engine,
      staged_ref)
    before <- scworkbench::wb_art_list(x$bundle)
    for (operation in names(staged_cases)) {
      result <- scworkbench::wb_session_run_op(setup$fit$session,
        x$bundle, paste0("scwbtrajectory::", operation),
        inputs = staged_cases[[operation]])
      key <- .trajectory_s2_registered_model_oracle_key(operation,
        engine, staged_cases[[operation]]$lineage)
      trajectory_s2_compare_outputs(x$bundle, key,
        list(outputs = result$run$outputs), reference)
      stopifnot(identical(before, scworkbench::wb_art_list(x$bundle)))
    }
    if (identical(engine, "monocle3")) {
      alternate <- .trajectory_s2_registered_model_cases(x, engine,
        staged_ref, lineage = "root:A::lineage:D")
      for (operation in names(alternate)) {
        result <- scworkbench::wb_session_run_op(setup$fit$session,
          x$bundle, paste0("scwbtrajectory::", operation),
          inputs = alternate[[operation]])
        key <- .trajectory_s2_registered_model_oracle_key(operation,
          engine, alternate[[operation]]$lineage)
        trajectory_s2_compare_outputs(x$bundle, key,
          list(outputs = result$run$outputs), reference)
        stopifnot(identical(before, scworkbench::wb_art_list(x$bundle)))
      }
    }
    setup$producer$cleanup()
    trajectory_close(x)
  }
  TRUE
}

test_trajectory_s2_registered_model_promoted_oracles_contract <-
    function() {
  on.exit(options(scworkbench.promotion_workspace_publish_hook = NULL),
    add = TRUE)
  for (engine in c("slingshot", "monocle3")) {
    x <- trajectory_s2_fixture()
    setup <- .trajectory_s2_registered_model_setup(x, engine)
    promoted <- scworkbench::wb_session_promote_current(x$bundle,
      setup$fit$session, outputs = c("embedding", "model", "weights"),
      created_by = "s2-registered-model-oracle")
    model_id <- promoted$promoted$artifact_id[
      promoted$promoted$output == "model"][[1L]]
    signal_id <- promoted$promoted$artifact_id[
      promoted$promoted$output == "weights"][[1L]]
    model <- scworkbench::wb_art_load(x$bundle, model_id,
      expected_kind = "model_ref")
    stopifnot(inherits(model, "wb_primitive"),
      identical(model$engine, engine),
      !is.null(model$model_source_linkage),
      !is.null(model$model_tree_manifest),
      !is.null(model$model_obj))
    expected_axis <- getFromNamespace(
      ".wb_operation_pack_v2_identity_plain", "scworkbench")(
        "scworkbench.cell_order.v1", unname(as.list(x$cells)))
    stopifnot(identical(model$model_source_linkage$
      fitted_cell_order_identity, expected_axis))
    reference <- .trajectory_s2_registered_model_reference(engine)
    model_ref <- scworkbench::wb_ref_artifact(model_id, "model")
    signal_ref <- scworkbench::wb_ref_artifact(signal_id, "weights")
    cases <- .trajectory_s2_registered_model_cases(x, engine, model_ref,
      signal = signal_ref)
    if (identical(engine, "monocle3")) {
      cases <- c(cases, stats::setNames(
        .trajectory_s2_registered_model_cases(x, engine, model_ref,
          lineage = "root:A::lineage:D", signal = signal_ref),
        paste0(names(cases), "_alternate")))
    }
    terminal_ids <- character()
    for (case_name in names(cases)) {
      operation <- sub("_alternate$", "", case_name)
      result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
        x$bundle, paste0("scwbtrajectory::", operation),
        inputs = cases[[case_name]])
      key <- .trajectory_s2_registered_model_oracle_key(operation,
        engine, cases[[case_name]]$lineage)
      trajectory_s2_compare_outputs(x$bundle, key,
        list(outputs = result$run$outputs), reference)
      closure <- names(reference$cases[[key]]$outputs)
      if (identical(operation, "trajectory_geneset_dynamics"))
        closure <- c(closure, "adjusted_scores", "adjustment_diagnostics")
      stopifnot(setequal(names(result$run$outputs), closure))
      if (identical(engine, "slingshot")) {
        before <- .trajectory_s2_registered_durable_snapshot(x)
        options(scworkbench.promotion_workspace_publish_hook = function(...)
          stop("injected S2 model publication failure"))
        failure <- tryCatch(scworkbench::wb_session_promote_current(
          x$bundle, result$session, outputs = closure), error = identity)
        options(scworkbench.promotion_workspace_publish_hook = NULL)
        stopifnot(inherits(failure, "error"),
          identical(before, .trajectory_s2_registered_durable_snapshot(x)))
      }
      publication <- scworkbench::wb_session_promote_current(x$bundle,
        result$session, outputs = closure,
        created_by = "s2-registered-model-oracle")
      details <- scworkbench::wb_prov_run_details(x$bundle,
        result$run$run_id, publication$provenance_snapshot_id)
      stopifnot(setequal(publication$promoted$output, closure),
        all(publication$promoted$status == "ok"),
        all(details$outputs$verification_status == "content_verified"),
        all(details$outputs$replay_equivalence_policy ==
          "exact_canonical_value"))
      replay <- scworkbench::wb_prov_run_replay_inputs(x$bundle,
        result$run$run_id, publication$provenance_snapshot_id)
      stopifnot(inherits(replay, "wb_prov_run_replay_inputs"),
        is.list(replay$inputs),
        identical(replay$inputs$model$ref_type, "artifact"),
        identical(replay$inputs$model$artifact_id, model_id))
      terminal_ids <- c(terminal_ids, publication$promoted$artifact_id)
    }
    stopifnot(length(terminal_ids) == if (identical(engine, "slingshot"))
      20L else 40L,
      all(vapply(terminal_ids, function(id)
        inherits(scworkbench::wb_art_load(x$bundle, id), "wb_primitive"),
        logical(1))))
    before_missing <- .trajectory_s2_registered_durable_snapshot(x)
    missing <- .trajectory_s2_registered_without_pack(tryCatch(
      scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
        x$bundle, "scwbtrajectory::trajectory_de",
        inputs = cases$trajectory_de), error = identity))
    stopifnot(inherits(missing, "wb_operation_pack_unavailable"),
      identical(before_missing, .trajectory_s2_registered_durable_snapshot(x)))
    stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
      x$bundle)$ready))
    scworkbench::wb_close_bundle(x$bundle); x$bundle <- NULL
    relocated <- tempfile(paste0("s2-", engine, "-relocated-"))
    fs::dir_copy(x$root, relocated)
    reopened <- scworkbench::wb_open_bundle(relocated, read_only = TRUE)
    .trajectory_s2_registered_without_pack({
      stopifnot(inherits(scworkbench::wb_art_load(reopened, model_id,
        expected_kind = "model_ref"), "wb_primitive"),
        all(vapply(terminal_ids, function(id)
          inherits(scworkbench::wb_art_load(reopened, id), "wb_primitive"),
          logical(1))),
        isTRUE(scworkbench::wb_bundle_validation_report(relocated,
          "recipient")$ok))
    })
    scworkbench::wb_close_bundle(reopened)
    unlink(relocated, recursive = TRUE, force = TRUE)
    setup$producer$cleanup()
    trajectory_close(x)
  }
  TRUE
}

test_trajectory_s2_registered_model_rejection_residue_contract <- function() {
  x <- trajectory_s2_fixture(); on.exit(trajectory_close(x), add = TRUE)
  setup <- .trajectory_s2_registered_model_setup(x, "slingshot")
  on.exit(setup$producer$cleanup(), add = TRUE)
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    setup$fit$session, outputs = c("embedding", "model", "weights"))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  model <- scworkbench::wb_art_load(x$bundle, model_id,
    expected_kind = "model_ref")
  model_ref <- scworkbench::wb_ref_artifact(model_id, "model")
  cases <- .trajectory_s2_registered_model_cases(x, "slingshot",
    model_ref, lineage = "missing")
  before <- .trajectory_s2_registered_durable_snapshot(x)
  for (operation in names(cases)) {
    rejected <- tryCatch(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), x$bundle,
      paste0("scwbtrajectory::", operation), inputs = cases[[operation]]),
      error = identity)
    stopifnot(inherits(rejected, "error"),
      identical(before, .trajectory_s2_registered_durable_snapshot(x)))
  }
  valid <- .trajectory_s2_registered_model_cases(x, "slingshot",
    model_ref)$trajectory_de
  attack <- function(path) {
    clean <- .trajectory_s2_registered_durable_snapshot(x)
    original <- readBin(path, "raw", n = file.info(path)$size)
    mode <- file.info(path)$mode
    Sys.chmod(path, mode = "0644")
    restore <- function() {
      con <- file(path, "wb"); writeBin(original, con); close(con)
      Sys.chmod(path, mode = mode)
    }
    on.exit(restore(), add = TRUE)
    con <- file(path, "ab"); writeBin(charToRaw("altered"), con); close(con)
    tampered <- .trajectory_s2_registered_durable_snapshot(x)
    stopifnot(!identical(clean$payload, tampered$payload))
    rejected <- tryCatch(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), x$bundle,
      "scwbtrajectory::trajectory_de", inputs = valid),
      error = identity)
    stopifnot(inherits(rejected, "error"),
      identical(tampered, .trajectory_s2_registered_durable_snapshot(x)))
    restore()
    stopifnot(identical(clean, .trajectory_s2_registered_durable_snapshot(x)))
    on.exit(NULL, add = FALSE)
  }
  attack(file.path(model$model_tree_root, model$model_primary_path))
  attack(file.path(model$model_tree_root,
    model$model_source_linkage$selected_axis_evidence$relative_path))
  foreign <- trajectory_s2_g1_test_producer("slingshot", "foreign_model",
    .trajectory_s2_synthetic_slingshot)
  on.exit(foreign$cleanup(), add = TRUE)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, foreign$fit_op, inputs = list(source = "latent",
      dims = c(1L, 3L), fitted_cells = x$cells))
  publication <- scworkbench::wb_session_promote_current(x$bundle,
    fit$session, outputs = c("embedding", "model", "weights"))
  foreign_id <- publication$promoted$artifact_id[
    publication$promoted$output == "model"][[1L]]
  invalid_type <- .trajectory_s2_registered_model_cases(x, "slingshot",
    scworkbench::wb_ref_artifact(foreign_id, "model"))$trajectory_de
  before_type <- .trajectory_s2_registered_durable_snapshot(x)
  rejected_type <- tryCatch(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_de", inputs = invalid_type),
    error = identity)
  stopifnot(inherits(rejected_type, "error"),
    identical(before_type, .trajectory_s2_registered_durable_snapshot(x)))
  TRUE
}
