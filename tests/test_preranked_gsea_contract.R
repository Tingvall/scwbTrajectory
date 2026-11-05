.trajectory_gsea_fixture <- function(n_cells = 80L) {
  cells <- sprintf("gsea_cell_%03d", seq_len(n_cells))
  features <- sprintf("gene_%02d", seq_len(48L))
  time <- seq(0, 1, length.out = n_cells)
  values <- matrix(0, length(features), n_cells,
    dimnames = list(features, cells))
  wiggle <- outer(seq_along(features), seq_along(cells), function(i, j)
    0.03 * sin(i * 0.7 + j * 0.19))
  values[1:6, ] <- 8 * rep(time, each = 6L) + wiggle[1:6, ]
  values[7:12, ] <- 8 * rep(1 - time, each = 6L) + wiggle[7:12, ]
  peak <- exp(-((time - 0.42) / 0.07)^2)
  trough <- exp(-((time - 0.72) / 0.07)^2)
  values[13:18, ] <- 8 * rep(peak, each = 6L) + wiggle[13:18, ]
  values[19:24, ] <- 6 - 5 * rep(trough, each = 6L) + wiggle[19:24, ]
  values[25:45, ] <- 2 + wiggle[25:45, ]
  values[41, ] <- values[40, ]
  values[42, ] <- values[40, ]
  values[46, ] <- 0
  values[47, ] <- 3
  values[48, 5L] <- NA_real_
  set <- function(name, members) list(namespace = "synthetic",
    feature_set = name, memberships = lapply(members, function(feature_id)
      list(feature_id = feature_id)))
  registry <- list(matrix_id = "RNA",
    feature_universe = as.list(features), sets = list(
      set("increasing", features[1:6]),
      set("decreasing", features[7:12]),
      set("transient_peak", features[13:18]),
      set("transient_trough", features[19:24]),
      set("null", features[25:30]),
      set("too_small", features[1]),
      set("no_overlap", c("not_in_matrix_1", "not_in_matrix_2"))))
  list(cells = cells, features = features, values = values, time = time,
    pseudotime = scworkbench::wb_prim_cell_vector_family(
      rep("lineage", n_cells), cells, time, value_type = "continuous"),
    samples = scworkbench::wb_prim_cell_vector(stats::setNames(
      rep(paste0("sample_", 1:4), length.out = n_cells), cells),
      value_type = "categorical"), registry = registry)
}

.trajectory_gsea_call <- function(x, ...) {
  method <- getFromNamespace(".trajectory_preranked_gsea",
    "scwbTrajectory")
  args <- list(matrix = x$values, registry = x$registry,
    effective_cells = x$cells, pseudotime = x$pseudotime,
    sample_id = x$samples, namespace = "synthetic",
    pseudotime_key = "lineage", window_width_fraction = 0.20,
    window_step_fraction = 0.10, min_cells_per_window = 8L,
    min_samples_per_window = 2L, n_permutations = 39L, seed = 53L,
    min_set_size = 3L, max_set_size = 10L, sample_size = 31L,
    n_perm_simple = 100L, eps = 1e-8, global_q_cutoff = 0.25,
    local_q_cutoff = 1, merge_gap_fraction = 0.10)
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(method, args)
}

test_trajectory_preranked_gsea_optimized_kernel_parity_contract <- function() {
  monotonic <- getFromNamespace(".trajectory_gsea_monotonic",
    "scwbTrajectory")
  plan_fun <- getFromNamespace(".trajectory_gsea_permutation_plan",
    "scwbTrajectory")
  monotonic_batched <- getFromNamespace(
    ".trajectory_gsea_monotonic_batched", "scwbTrajectory")
  windows_fun <- getFromNamespace(".trajectory_gsea_windows",
    "scwbTrajectory")
  local_effects <- getFromNamespace(".trajectory_gsea_local_effects",
    "scwbTrajectory")
  contrast <- getFromNamespace(".trajectory_gsea_contrast",
    "scwbTrajectory")
  contrast_weights <- getFromNamespace(".trajectory_gsea_contrast_weights",
    "scwbTrajectory")
  window_weights <- getFromNamespace(".trajectory_gsea_window_weights",
    "scwbTrajectory")
  contrasts_batched <- getFromNamespace(
    ".trajectory_gsea_contrasts_batched", "scwbTrajectory")
  permute <- getFromNamespace(".trajectory_regions_permute",
    "scwbTrajectory")
  with_seed <- getFromNamespace(".trajectory_regions_with_seed",
    "scwbTrajectory")

  set.seed(4101L)
  expression <- matrix(stats::rnorm(37L * 80L), 37L, 80L)
  expression[1L, ] <- seq_len(80L)
  expression[2L, ] <- expression[1L, ]
  expression[3L, ] <- 1
  pseudotime <- sort(stats::runif(80L))
  blocks <- rep(sprintf("sample_%02d", 1:8), each = 10L)
  permutations <- 23L
  seed <- 991L
  for (one_blocks in list(NULL, blocks)) {
    one_plan <- plan_fun(length(pseudotime), one_blocks, permutations, seed)
    one_legacy <- with_seed(seed, lapply(seq_len(permutations),
      function(ii) permute(pseudotime, one_blocks)))
    stopifnot(all(vapply(seq_len(permutations), function(ii)
      identical(pseudotime[one_plan[, ii]], one_legacy[[ii]]), logical(1))))
  }
  plan <- plan_fun(length(pseudotime), blocks, permutations, seed)
  legacy_pseudotime <- with_seed(seed, lapply(seq_len(permutations),
    function(ii) permute(pseudotime, blocks)))

  legacy_sum <- legacy_sumsq <- numeric(nrow(expression))
  for (permuted in legacy_pseudotime) {
    value <- monotonic(expression, permuted)
    value[!is.finite(value)] <- 0
    legacy_sum <- legacy_sum + value
    legacy_sumsq <- legacy_sumsq + value^2
  }
  reference <- monotonic(expression, pseudotime)
  for (batch_size in c(1L, 4L, 32L)) for (chunk_size in c(1L, 7L, 128L)) {
    optimized <- monotonic_batched(expression, pseudotime, plan,
      batch_size, chunk_size)
    stopifnot(isTRUE(all.equal(optimized$observed, reference,
        tolerance = 1e-14)),
      isTRUE(all.equal(optimized$sum, legacy_sum, tolerance = 1e-12)),
      isTRUE(all.equal(optimized$sumsq, legacy_sumsq, tolerance = 1e-12)))
  }

  endpoint_weights <- contrast_weights(pseudotime, 0.2, 0.2)$weights
  endpoint_observed <- contrast(expression, pseudotime, 0.2, 0.2)$value
  endpoint_sum <- endpoint_sumsq <- numeric(nrow(expression))
  for (permuted in legacy_pseudotime) {
    value <- contrast(expression, permuted, 0.2, 0.2)$value
    value[!is.finite(value)] <- 0
    endpoint_sum <- endpoint_sum + value
    endpoint_sumsq <- endpoint_sumsq + value^2
  }
  endpoint <- contrasts_batched(expression, endpoint_weights, plan, 5L, 9L)
  stopifnot(isTRUE(all.equal(as.numeric(endpoint$observed),
      endpoint_observed, tolerance = 1e-12)),
    isTRUE(all.equal(as.numeric(endpoint$sum), endpoint_sum,
      tolerance = 1e-11)),
    isTRUE(all.equal(as.numeric(endpoint$sumsq), endpoint_sumsq,
      tolerance = 1e-11)))

  windows <- windows_fun(pseudotime, blocks, 0.2, 0.1, 5L, 1L)
  windows <- windows[vapply(windows, `[[`, logical(1), "supported")]
  weights <- window_weights(windows, length(pseudotime))
  legacy_observed <- local_effects(expression, windows)
  legacy_sum <- legacy_sumsq <- matrix(0, nrow(expression), length(windows))
  for (permuted in legacy_pseudotime) {
    permuted_windows <- windows_fun(permuted, blocks, 0.2, 0.1, 5L, 1L)
    ids <- vapply(windows, `[[`, character(1), "ranking_id")
    available <- vapply(permuted_windows, `[[`, character(1), "ranking_id")
    value <- local_effects(expression, permuted_windows[match(ids, available)])
    value[!is.finite(value)] <- 0
    legacy_sum <- legacy_sum + value
    legacy_sumsq <- legacy_sumsq + value^2
  }
  for (batch_size in c(1L, 3L, 8L)) for (chunk_size in c(1L, 9L, 128L)) {
    optimized <- contrasts_batched(expression, weights, plan,
      batch_size, chunk_size)
    stopifnot(isTRUE(all.equal(optimized$observed, legacy_observed,
        tolerance = 1e-12)),
      isTRUE(all.equal(optimized$sum, legacy_sum, tolerance = 1e-11)),
      isTRUE(all.equal(optimized$sumsq, legacy_sumsq, tolerance = 1e-11)))
  }
  TRUE
}

test_trajectory_preranked_gsea_streamed_chunk_replay_contract <- function() {
  x <- .trajectory_gsea_fixture()
  reference <- .trajectory_gsea_call(x, analysis_scope = "both",
    n_permutations = 11L, feature_chunk_size = 128L)
  for (chunk_size in c(1L, 7L, 19L)) {
    observed <- new.env(parent = emptyenv())
    observed$chunks <- list()
    reader <- function(features, cells) {
      observed$chunks[[length(observed$chunks) + 1L]] <- features
      list(values = x$values[features, cells, drop = FALSE])
    }
    streamed <- .trajectory_gsea_call(x, matrix = NULL,
      matrix_source = list(reader = reader, features = x$features,
        cells = x$cells), feature_chunk_size = chunk_size,
      analysis_scope = "both", n_permutations = 11L)
    stopifnot(identical(serialize(streamed$outputs, NULL, version = 3L),
        serialize(reference$outputs, NULL, version = 3L)),
      identical(unlist(observed$chunks, use.names = FALSE), x$features),
      all(lengths(observed$chunks) <= chunk_size))
  }
  TRUE
}

test_trajectory_preranked_gsea_manifest_contract <- function() {
  trajectory_activate_pack()
  manifest <- scworkbench::wb_operation_pack_inspect(
    "scwbtrajectory")$manifest
  operations <- stats::setNames(manifest$operations,
    vapply(manifest$operations, `[[`, character(1), "name"))
  operation <- operations$trajectory_preranked_gsea
  stopifnot(!is.null(operation),
    identical(operation$qualified_name,
      "scwbtrajectory::trajectory_preranked_gsea"),
    identical(operation$title, "Trajectory Preranked GSEA"),
    identical(operation$params_schema$required,
      list("namespace", "pseudotime_key")),
    identical(operation$params_schema$properties$analysis_scope$default,
      "both"),
    identical(vapply(operation$inputs, `[[`, character(1), "name"),
      c("matrix", "pseudotime", "membership", "sample_id", "covariates",
        "depth")),
    identical(vapply(operation$outputs, `[[`, character(1), "name"),
      c("global_enrichment", "window_enrichment",
        "enrichment_intervals", "leading_edge", "ranked_features",
        "fit_summary")),
    identical(operation$dependencies[[2L]]$name, "fgsea"),
    identical(operation$dependencies[[2L]]$version$specifier,
      "==1.34.2"),
    identical(length(operation$typed_slice_reads), 1L),
    identical(length(operation$registry_reads), 1L),
    identical(operation$params_schema$properties$adjust_depth$default,
      FALSE),
    identical(operation$params_schema$properties$depth_field$default,
      "nCount_RNA"),
    length(operation$provenance_reads) == 5L,
    all(vapply(operation$promotion_closure, function(value)
      setequal(as.character(unlist(value)),
        c("global_enrichment", "window_enrichment",
          "enrichment_intervals", "leading_edge", "ranked_features",
          "fit_summary")), logical(1))))
  TRUE
}

test_trajectory_preranked_gsea_scopes_and_signal_contract <- function() {
  x <- .trajectory_gsea_fixture()
  set.seed(712L)
  before <- .Random.seed
  both <- .trajectory_gsea_call(x, analysis_scope = "both")
  after <- .Random.seed
  repeat_both <- .trajectory_gsea_call(x, analysis_scope = "both")
  stopifnot(identical(before, after),
    identical(serialize(both$outputs, NULL, version = 3L),
      serialize(repeat_both$outputs, NULL, version = 3L)))
  global <- both$outputs$global_enrichment$table
  local <- both$outputs$window_enrichment$table
  intervals <- both$outputs$enrichment_intervals$table
  leading <- both$outputs$leading_edge$table
  increasing <- global[global$feature_set == "increasing", , drop = FALSE]
  decreasing <- global[global$feature_set == "decreasing", , drop = FALSE]
  peak <- local[local$feature_set == "transient_peak" &
    local$status == "ok", , drop = FALSE]
  trough <- local[local$feature_set == "transient_trough" &
    local$status == "ok", , drop = FALSE]
  stopifnot(nrow(increasing) == 1L, nrow(decreasing) == 1L,
    increasing$normalized_enrichment_score > 0,
    decreasing$normalized_enrichment_score < 0,
    peak$window_start_fraction[which.max(
      peak$normalized_enrichment_score)] <= 0.5,
    peak$window_end_fraction[which.max(
      peak$normalized_enrichment_score)] >= 0.35,
    trough$window_start_fraction[which.min(
      trough$normalized_enrichment_score)] <= 0.8,
    trough$window_end_fraction[which.min(
      trough$normalized_enrichment_score)] >= 0.65,
    nrow(intervals) > 0L,
    all(intervals$n_constituent_windows >= 1L),
    any(intervals$n_constituent_windows > 1L),
    all(nzchar(intervals$constituent_ranking_ids)),
    any(leading$analysis_scope == "local" & leading$peak_window),
    any(leading$analysis_scope == "local" &
      leading$consensus_leading_edge),
    all(nzchar(leading$interval_id[leading$analysis_scope == "local" &
      (leading$peak_window | leading$consensus_leading_edge)])),
    identical(global$q_global_family[global$status == "ok"],
      stats::p.adjust(global$p_value[global$status == "ok"], method = "BH")),
    identical(local$q_global_family[local$status == "ok"],
      stats::p.adjust(local$p_value[local$status == "ok"], method = "BH")))
  for (ranking_id in unique(local$ranking_id)) {
    rows <- local$ranking_id == ranking_id & local$status == "ok" &
      is.finite(local$p_value)
    stopifnot(identical(local$q_within_ranking[rows],
      stats::p.adjust(local$p_value[rows], method = "BH")))
  }

  global_only <- .trajectory_gsea_call(x, analysis_scope = "global",
    n_permutations = 11L)
  local_only <- .trajectory_gsea_call(x, analysis_scope = "local",
    n_permutations = 11L)
  stopifnot(nrow(global_only$outputs$window_enrichment$table) == 0L,
    identical(global_only$outputs$window_enrichment$meta$status,
      "not_requested"),
    nrow(local_only$outputs$global_enrichment$table) == 0L,
    identical(local_only$outputs$global_enrichment$meta$status,
      "not_requested"))
  TRUE
}

test_trajectory_preranked_gsea_leading_edge_and_parity_contract <- function() {
  x <- .trajectory_gsea_fixture()
  result <- .trajectory_gsea_call(x, analysis_scope = "global")
  enrichment <- result$outputs$global_enrichment$table
  leading <- result$outputs$leading_edge$table
  ranked <- result$outputs$ranked_features$table
  summary <- result$outputs$fit_summary$table
  ranking <- ranked[ranked$eligible, , drop = FALSE]
  ranking <- ranking[order(ranking$gene_rank), , drop = FALSE]
  statistic <- stats::setNames(ranking$standardized_ranking_statistic,
    ranking$feature_id)
  sets <- lapply(x$registry$sets, function(item) vapply(item$memberships,
    `[[`, character(1), "feature_id"))
  names(sets) <- vapply(x$registry$sets, `[[`, character(1), "feature_set")
  eligible_sets <- sets[vapply(sets, function(value) {
    n <- length(intersect(value, names(statistic)))
    n >= 3L && n <= 10L
  }, logical(1))]
  seed <- summary$fgsea_seed[summary$analysis_scope == "global"]
  with_seed <- getFromNamespace(".trajectory_regions_with_seed",
    "scwbTrajectory")
  direct <- with_seed(seed, fgsea::fgseaMultilevel(pathways = eligible_sets,
    stats = statistic, sampleSize = 31L, minSize = 3L, maxSize = 10L,
    eps = 1e-8, scoreType = "std", nproc = 1L,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE), gseaParam = 1,
    nPermSimple = 100L))
  direct <- as.data.frame(direct)
  observed <- enrichment[match(direct$pathway, enrichment$feature_set),
    , drop = FALSE]
  stopifnot(identical(observed$enrichment_score, as.numeric(direct$ES)),
    identical(observed$normalized_enrichment_score, as.numeric(direct$NES)),
    identical(observed$p_value, as.numeric(direct$pval)))
  for (ii in seq_len(nrow(direct))) {
    observed_genes <- leading$feature_id[leading$feature_set ==
      direct$pathway[ii] & leading$analysis_scope == "global"]
    stopifnot(identical(observed_genes,
      as.character(direct$leadingEdge[[ii]])))
  }
  stopifnot(all(c("peak_window", "consensus_leading_edge",
      "leading_edge_boundary_tied") %in% names(leading)))
  TRUE
}

test_trajectory_preranked_gsea_frozen_fgsea_parity_contract <- function() {
  directory <- if (exists("test_dir", inherits = TRUE))
    get("test_dir", inherits = TRUE) else "tests"
  fixture <- jsonlite::read_json(file.path(directory, "fixtures",
    "trajectory_preranked_gsea_frozen.json"), simplifyVector = TRUE)
  feature_id <- as.character(fixture$feature_id)
  statistic <- as.numeric(fixture$statistic)
  ranked <- data.frame(analysis_scope = "global", namespace = "frozen",
    lineage_id = "lineage", ranking_id = "frozen_global",
    feature_id = feature_id, feature_index = seq_along(feature_id),
    gene_rank = seq_along(feature_id), raw_ranking_statistic = statistic,
    standardized_ranking_statistic = statistic, eligible = TRUE,
    exclusion_reason = "", tie_group_id = seq_along(feature_id),
    tie_group_size = 1L, stringsAsFactors = FALSE)
  ranking <- list(table = ranked,
    stats = stats::setNames(statistic, feature_id))
  pathways <- lapply(fixture$pathways, as.character)
  sets <- lapply(names(pathways), function(name) list(feature_set = name,
    raw = pathways[[name]], overlap = pathways[[name]]))
  names(sets) <- names(pathways)
  engine <- getFromNamespace(".trajectory_gsea_fgsea", "scwbTrajectory")
  observed <- engine(ranking, sets, "global", "frozen", "lineage",
    list(n_cells = 20L, n_samples = NA_integer_), "std", 1, 31L, 100L,
    1e-8, 3L, 10L, "BH", as.integer(fixture$seed))
  with_seed <- getFromNamespace(".trajectory_regions_with_seed",
    "scwbTrajectory")
  direct <- with_seed(as.integer(fixture$seed), fgsea::fgseaMultilevel(
    pathways = pathways, stats = ranking$stats, sampleSize = 31L,
    minSize = 3L, maxSize = 10L, eps = 1e-8, scoreType = "std",
    nproc = 1L, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE),
    gseaParam = 1, nPermSimple = 100L))
  direct <- as.data.frame(direct)
  result <- observed$enrichment[match(direct$pathway,
    observed$enrichment$feature_set), , drop = FALSE]
  stopifnot(identical(result$enrichment_score, as.numeric(direct$ES)),
    identical(result$normalized_enrichment_score, as.numeric(direct$NES)),
    identical(result$p_value, as.numeric(direct$pval)),
    identical(result$q_within_ranking, as.numeric(direct$padj)))
  for (ii in seq_len(nrow(direct))) stopifnot(identical(
    observed$leading$feature_id[observed$leading$feature_set ==
      direct$pathway[ii]], as.character(direct$leadingEdge[[ii]])))
  tied <- ranking
  tied$stats[["f04"]] <- tied$stats[["f03"]]
  tied$table$standardized_ranking_statistic[tied$table$feature_id == "f04"] <-
    tied$stats[["f03"]]
  tied_set <- list(boundary = list(feature_set = "boundary",
    raw = c("f01", "f02", "f03"), overlap = c("f01", "f02", "f03")))
  tied_result <- suppressWarnings(engine(tied, tied_set, "global", "frozen",
    "lineage", list(n_cells = 20L, n_samples = NA_integer_), "std", 1,
    31L, 100L, 1e-8, 3L, 10L, "BH", as.integer(fixture$seed) + 1L))
  stopifnot(nrow(tied_result$leading) > 0L,
    all(tied_result$leading$leading_edge_boundary_tied))
  TRUE
}

test_trajectory_preranked_gsea_exclusion_inversion_and_rejection_contract <-
    function() {
  x <- .trajectory_gsea_fixture()
  forward <- .trajectory_gsea_call(x, analysis_scope = "global",
    n_permutations = 19L)
  inverse <- .trajectory_gsea_call(x, analysis_scope = "global",
    n_permutations = 19L, pseudotime_direction = "invert")
  f <- forward$outputs$global_enrichment$table
  i <- inverse$outputs$global_enrichment$table
  stopifnot(sign(f$normalized_enrichment_score[f$feature_set ==
      "increasing"]) == -sign(i$normalized_enrichment_score[
        i$feature_set == "increasing"]))
  ranked <- forward$outputs$ranked_features$table
  stopifnot(identical(ranked$exclusion_reason[ranked$feature_id ==
      "gene_46"], "undetected"),
    identical(ranked$exclusion_reason[ranked$feature_id == "gene_47"],
      "zero_information"),
    identical(ranked$exclusion_reason[ranked$feature_id == "gene_48"],
      "non_finite_expression"),
    any(ranked$tie_group_size > 1L, na.rm = TRUE))
  method <- getFromNamespace(".trajectory_preranked_gsea",
    "scwbTrajectory")
  duplicated <- x$values
  rownames(duplicated)[2L] <- rownames(duplicated)[1L]
  bad_duplicate <- try(.trajectory_gsea_call(
    within(x, values <- duplicated), analysis_scope = "global",
    n_permutations = 9L), silent = TRUE)
  missing <- x
  missing$pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("lineage", length(x$cells)), x$cells, rep(NA_real_, length(x$cells)),
    value_type = "continuous")
  bad_missing <- try(.trajectory_gsea_call(missing,
    analysis_scope = "global", n_permutations = 9L), silent = TRUE)
  stopifnot(inherits(bad_duplicate, "try-error"),
    inherits(bad_missing, "try-error"))
  TRUE
}

test_trajectory_preranked_gsea_support_permutation_and_adjustment_contract <-
    function() {
  x <- .trajectory_gsea_fixture()
  blocked <- .trajectory_gsea_call(x, analysis_scope = "global",
    global_rank_method = "early_late_contrast", permutation_mode = "blocked",
    n_permutations = 19L)
  unblocked <- .trajectory_gsea_call(x, analysis_scope = "global",
    sample_id = NULL, min_samples_per_window = 1L,
    permutation_mode = "unblocked", n_permutations = 19L)
  blocked_summary <- blocked$outputs$fit_summary$table
  unblocked_summary <- unblocked$outputs$fit_summary$table
  stopifnot(identical(blocked_summary$permutation_mode[1L], "blocked"),
    identical(unblocked_summary$permutation_mode[1L], "unblocked"),
    blocked_summary$early_fraction[1L] == 0.2,
    blocked_summary$late_fraction[1L] == 0.2,
    blocked_summary$early_n_cells[1L] >= 8L,
    blocked_summary$late_n_cells[1L] >= 8L,
    blocked_summary$early_n_samples[1L] >= 2L,
    blocked_summary$late_n_samples[1L] >= 2L)
  increasing <- blocked$outputs$global_enrichment$table
  increasing <- increasing[increasing$feature_set == "increasing", ]
  stopifnot(increasing$normalized_enrichment_score > 0)

  covariate <- scworkbench::wb_prim_cell_vector_family(
    rep("technical", length(x$cells)), x$cells,
    sin(seq(0, 4 * pi, length.out = length(x$cells))),
    value_type = "continuous")
  adjusted <- .trajectory_gsea_call(x, analysis_scope = "global",
    covariates = covariate, covariate_adjustment = "linear",
    n_permutations = 19L)
  stopifnot(identical(adjusted$outputs$fit_summary$table$
    covariate_adjustment[1L], "linear"))

  no_sample_support <- try(.trajectory_gsea_call(x,
    analysis_scope = "local", sample_id = NULL,
    min_samples_per_window = 2L, n_permutations = 9L), silent = TRUE)
  no_cell_support <- try(.trajectory_gsea_call(x,
    analysis_scope = "local", min_cells_per_window = length(x$cells),
    n_permutations = 9L), silent = TRUE)
  endpoint_support <- try(.trajectory_gsea_call(x,
    analysis_scope = "global", global_rank_method = "early_late_contrast",
    early_fraction = 0.05, late_fraction = 0.05,
    min_cells_per_window = 10L, n_permutations = 9L), silent = TRUE)
  stopifnot(inherits(no_sample_support, "try-error"),
    inherits(no_cell_support, "try-error"),
    inherits(endpoint_support, "try-error"))
  TRUE
}

test_trajectory_preranked_gsea_registered_lifecycle_contract <- function() {
  x <- trajectory_fixture(); on.exit(trajectory_close(x), add = TRUE)
  before <- list(artifacts = scworkbench::wb_art_list(x$bundle),
    fingerprint = getFromNamespace(".wb_session_source_identity",
      "scworkbench")(x$bundle)$scientific_state_fingerprint)
  missing_pseudotime <- scworkbench::wb_prim_cell_vector_family(
    rep("A", length(x$cells)), x$cells, rep(NA_real_, length(x$cells)),
    value_type = "continuous")
  rejected <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle,
    "scwbtrajectory::trajectory_preranked_gsea", inputs = list(
      matrix = "RNA", pseudotime = missing_pseudotime,
      namespace = "trajectory", pseudotime_key = "A",
      analysis_scope = "global", min_cells_per_window = 3L,
      n_permutations = 9L, min_set_size = 1L, max_set_size = 5L,
      sample_size = 31L, n_perm_simple = 100L, eps = 1e-8,
      chunk_size = 3L, min_cells = 1L)), silent = TRUE)
  after_rejection <- list(artifacts = scworkbench::wb_art_list(x$bundle),
    fingerprint = getFromNamespace(".wb_session_source_identity",
      "scworkbench")(x$bundle)$scientific_state_fingerprint)
  stopifnot(inherits(rejected, "try-error"),
    identical(before, after_rejection))
  result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbtrajectory::trajectory_preranked_gsea", inputs = list(
      matrix = "RNA", pseudotime = x$pseudotime_family,
      namespace = "trajectory", pseudotime_key = "A",
      analysis_scope = "global", min_cells_per_window = 3L,
      n_permutations = 19L, min_set_size = 1L, max_set_size = 5L,
      sample_size = 31L, n_perm_simple = 100L, eps = 1e-8,
      chunk_size = 3L, min_cells = 1L))
  after <- list(artifacts = scworkbench::wb_art_list(x$bundle),
    fingerprint = getFromNamespace(".wb_session_source_identity",
      "scworkbench")(x$bundle)$scientific_state_fingerprint)
  stopifnot(identical(before, after),
    identical(names(result$run$outputs), c("global_enrichment",
      "window_enrichment", "enrichment_intervals", "leading_edge",
      "ranked_features", "fit_summary")),
    nrow(result$run$outputs$window_enrichment$table) == 0L)
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    result$session, outputs = names(result$run$outputs),
    created_by = "trajectory-gsea-contract")
  finalized <- scworkbench::wb_finalize_recipient_bundle(x$bundle)
  scworkbench::wb_close_bundle(x$bundle)
  x$bundle <- NULL
  recipient <- scworkbench::wb_bundle_validation_report(x$root, "recipient")
  x$bundle <- scworkbench::wb_open_bundle(x$root, read_only = TRUE)
  stopifnot(all(promoted$promoted$status == "ok"), isTRUE(finalized$ready),
    isTRUE(recipient$ok))
  TRUE
}
