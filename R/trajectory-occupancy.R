.trajectory_occupancy_vector <- function(value, cells, name,
                                         value_type = "categorical") {
  .trajectory_kind(value, "cell_vector", paste0(" for ", name))
  .trajectory_assert(identical(as.character(value$value_type)[1L], value_type),
    paste0(name, " must be a ", value_type, " cell_vector."))
  ids <- as.character(value$cell_id)
  .trajectory_assert(length(ids) == length(value$value) &&
      !anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids) &&
      setequal(ids, cells),
    paste0(name, " must cover the complete pseudotime cell axis exactly."))
  value$value[match(cells, ids)]
}

.trajectory_occupancy_strata <- function(group, sample = NULL) {
  groups <- sort(unique(group), method = "radix")
  out <- lapply(groups, function(name) list(
    stratum = "group", group_id = name, sample_id = NA_character_,
    index = which(group == name)))
  if (!is.null(sample)) for (name in groups) {
    samples <- sort(unique(sample[group == name]), method = "radix")
    out <- c(out, lapply(samples, function(sample_name) list(
      stratum = "group_sample", group_id = name, sample_id = sample_name,
      index = which(group == name & sample == sample_name))))
  }
  out
}

.trajectory_occupancy_compute <- function(pseudotime, group_id, sample_id,
    mapping_acceptance, params) {
  .trajectory_kind(pseudotime, "cell_vector_family", " for pseudotime")
  .trajectory_assert(identical(as.character(pseudotime$value_type)[1L],
      "continuous"), "trajectory_occupancy requires continuous pseudotime.")
  pseudotime_key <- as.character(params$pseudotime_key %||% "")
  acceptance_key <- as.character(params$acceptance_output_key %||% "")
  analysis_mode <- as.character(params$analysis_mode %||% "descriptive")
  n_bins <- as.integer(params$n_bins %||% 10L)
  .trajectory_assert(length(pseudotime_key) == 1L && nzchar(pseudotime_key),
    "pseudotime_key must select one non-empty pseudotime key.")
  .trajectory_assert(length(acceptance_key) == 1L && nzchar(acceptance_key),
    "acceptance_output_key must select one non-empty acceptance key.")
  .trajectory_assert(identical(analysis_mode, "descriptive"),
    paste("trajectory_occupancy is descriptive-only; inferential comparisons",
      "require a separately qualified replicate-level method."))
  .trajectory_assert(length(n_bins) == 1L && !is.na(n_bins) &&
      n_bins >= 2L && n_bins <= 1000L,
    "n_bins must be one integer from 2 through 1000.")

  keys <- as.character(pseudotime$key)
  selected <- which(keys == pseudotime_key)
  .trajectory_assert(length(selected) > 0L,
    "pseudotime_key is absent from the pseudotime family.")
  cells <- as.character(pseudotime$cell_id[selected])
  .trajectory_assert(!anyNA(cells) && all(nzchar(cells)) &&
      !anyDuplicated(cells),
    "Selected pseudotime must contain one unique row per cell.")
  raw_pseudotime <- pseudotime$value[selected]
  values <- suppressWarnings(as.numeric(raw_pseudotime))
  .trajectory_assert(!any(is.na(values) & !is.na(raw_pseudotime)) &&
      all(is.na(values) | is.finite(values)),
    "Pseudotime values must be finite numeric values or missing.")

  group <- as.character(.trajectory_occupancy_vector(group_id, cells,
    "group_id"))
  group[is.na(group) | !nzchar(group)] <- "<missing>"
  sample <- NULL
  if (!is.null(sample_id)) {
    sample <- as.character(.trajectory_occupancy_vector(sample_id, cells,
      "sample_id"))
    sample[is.na(sample) | !nzchar(sample)] <- "<missing>"
  }

  .trajectory_kind(mapping_acceptance, "obs_table",
    " for mapping_acceptance")
  acceptance <- mapping_acceptance$table
  required <- c("query_cell_id", "output_key", "accepted",
    "acceptance_reason")
  .trajectory_assert(is.data.frame(acceptance) &&
      all(required %in% names(acceptance)),
    paste("mapping_acceptance must contain query_cell_id, output_key,",
      "accepted, and acceptance_reason."))
  acceptance <- acceptance[as.character(acceptance$output_key) ==
    acceptance_key, , drop = FALSE]
  acceptance_cells <- as.character(acceptance$query_cell_id)
  .trajectory_assert(nrow(acceptance) == length(cells) &&
      !anyNA(acceptance_cells) && all(nzchar(acceptance_cells)) &&
      !anyDuplicated(acceptance_cells) && setequal(acceptance_cells, cells),
    paste("Selected mapping acceptance must contain exactly one row for every",
      "pseudotime cell."))
  acceptance <- acceptance[match(cells, acceptance_cells), , drop = FALSE]
  accepted <- acceptance$accepted
  .trajectory_assert(is.logical(accepted) && !anyNA(accepted),
    "mapping_acceptance accepted values must be non-missing logical values.")
  reason <- as.character(acceptance$acceptance_reason)
  .trajectory_assert(!anyNA(reason) && all(nzchar(reason)) &&
      all((accepted & reason == "accepted") |
        (!accepted & reason != "accepted")),
    "mapping_acceptance accepted flags and reasons are inconsistent.")

  finite <- is.finite(values)
  binned <- accepted & finite
  .trajectory_assert(any(binned),
    "trajectory_occupancy requires at least one accepted finite pseudotime.")
  lower <- min(values[binned]); upper <- max(values[binned])
  .trajectory_assert(upper > lower,
    "Accepted finite pseudotime must contain at least two distinct values.")
  breaks <- seq(lower, upper, length.out = n_bins + 1L)
  bin <- rep(NA_integer_, length(cells))
  bin[binned] <- findInterval(values[binned], breaks,
    rightmost.closed = TRUE, all.inside = TRUE)
  .trajectory_assert(all(bin[binned] >= 1L & bin[binned] <= n_bins),
    "Internal pseudotime bin assignment failed.")

  strata <- .trajectory_occupancy_strata(group, sample)
  occupancy_rows <- list(); cursor <- 0L
  for (stratum in strata) {
    denominator <- sum(binned[stratum$index])
    for (i in seq_len(n_bins)) {
      count <- sum(bin[stratum$index] == i, na.rm = TRUE)
      cursor <- cursor + 1L
      occupancy_rows[[cursor]] <- data.frame(
        stratum = stratum$stratum, group_id = stratum$group_id,
        sample_id = stratum$sample_id, bin_index = as.integer(i),
        bin_lower = breaks[[i]], bin_upper = breaks[[i + 1L]],
        lower_inclusive = TRUE, upper_inclusive = identical(i, n_bins),
        n_cells = as.integer(count),
        denominator_accepted_finite = as.integer(denominator),
        fraction_within_stratum = if (denominator > 0L)
          count / denominator else NA_real_, stringsAsFactors = FALSE)
    }
  }
  occupancy <- do.call(rbind, occupancy_rows); rownames(occupancy) <- NULL

  accounting_status <- ifelse(!accepted, reason,
    ifelse(finite, "accepted_binned", "accepted_pseudotime_missing"))
  status_levels <- c("accepted_binned", "accepted_pseudotime_missing",
    sort(unique(reason[!accepted]), method = "radix"))
  accounting_rows <- list(); cursor <- 0L
  for (stratum in strata) for (status in status_levels) {
    count <- sum(accounting_status[stratum$index] == status)
    denominator <- length(stratum$index)
    cursor <- cursor + 1L
    accounting_rows[[cursor]] <- data.frame(
      stratum = stratum$stratum, group_id = stratum$group_id,
      sample_id = stratum$sample_id, accounting_status = status,
      n_cells = as.integer(count), n_cells_in_stratum = as.integer(denominator),
      fraction_within_stratum = count / denominator,
      stringsAsFactors = FALSE)
  }
  accounting <- do.call(rbind, accounting_rows); rownames(accounting) <- NULL

  aggregate_accounting <- accounting[accounting$stratum == "group", ,
    drop = FALSE]
  .trajectory_assert(sum(aggregate_accounting$n_cells) == length(cells),
    "Acceptance accounting does not close over the complete cell axis.")
  summary <- data.frame(
    pseudotime_key = pseudotime_key,
    acceptance_output_key = acceptance_key,
    analysis_mode = analysis_mode,
    n_bins = n_bins, bin_lower = lower, bin_upper = upper,
    n_cells_input = as.integer(length(cells)),
    n_cells_accepted = as.integer(sum(accepted)),
    n_cells_withheld = as.integer(sum(!accepted)),
    n_cells_accepted_finite = as.integer(sum(binned)),
    n_cells_accepted_pseudotime_missing = as.integer(sum(accepted & !finite)),
    n_groups = as.integer(length(unique(group))),
    n_samples = as.integer(if (is.null(sample)) 0L else length(unique(sample))),
    accounting_complete = TRUE, stringsAsFactors = FALSE)

  meta <- list(operation = "trajectory_occupancy",
    pseudotime_key = pseudotime_key,
    acceptance_output_key = acceptance_key,
    analysis_mode = analysis_mode,
    bin_rule = "global_equal_width_left_closed_final_right_closed")
  list(
    bin_occupancy = scworkbench::wb_prim_obs_table(occupancy,
      table_type = "trajectory_occupancy_bins",
      entity_kind = "trajectory_occupancy_bin", meta = meta),
    acceptance_accounting = scworkbench::wb_prim_obs_table(accounting,
      table_type = "trajectory_occupancy_acceptance_accounting",
      entity_kind = "trajectory_occupancy_status", meta = meta),
    run_summary = scworkbench::wb_prim_obs_table(summary,
      table_type = "trajectory_occupancy_summary",
      entity_kind = "trajectory_occupancy_run", meta = meta))
}

.run_trajectory_occupancy <- function(readers, writers, params, context) {
  sample <- if (is.function(readers$sample_id)) readers$sample_id() else NULL
  outputs <- .trajectory_occupancy_compute(readers$pseudotime(),
    readers$group_id(), sample, readers$mapping_acceptance(), params)
  for (name in names(outputs))
    writers$write_terminal_primitive(name, outputs[[name]])
  invisible(list(status = "ok"))
}
