# Interpret a single authenticated feature-set-registry snapshot. Core owns
# the snapshot identity and validation; these are method selection rules.
.trajectory_s2_registry_sets <- function(snapshot, namespace, feature_sets,
                                         universe, min_set_size,
                                         max_set_size) {
  .trajectory_assert(is.list(snapshot) && is.list(snapshot$sets),
    "trajectory_geneset_dynamics: feature-set snapshot is unavailable.")
  matches <- Filter(function(item) identical(
    as.character(item$namespace %||% "")[1L], namespace), snapshot$sets)
  .trajectory_assert(length(matches) > 0L,
    paste0("trajectory_geneset_dynamics: no feature sets found for namespace='",
      namespace, "'."))
  names_all <- unique(vapply(matches, function(item)
    as.character(item$feature_set)[1L], character(1)))
  names_all <- names_all[nzchar(names_all)]
  if (!is.null(feature_sets)) {
    requested <- unique(as.character(unlist(feature_sets,
      use.names = FALSE)))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    .trajectory_assert(length(requested) > 0L,
      "trajectory_geneset_dynamics: feature_sets is empty after cleaning.")
    .trajectory_assert(!length(setdiff(requested, names_all)),
      "trajectory_geneset_dynamics: requested feature_set not found.")
    names_all <- requested
  }
  .trajectory_assert(length(names_all) > 0L,
    "trajectory_geneset_dynamics: no feature sets remain after filtering.")
  members <- list()
  for (name in names_all) {
    item <- Filter(function(row) identical(
      as.character(row$feature_set)[1L], name), matches)
    .trajectory_assert(length(item) == 1L,
      "trajectory_geneset_dynamics: feature-set name is ambiguous.")
    features <- vapply(item[[1L]]$memberships %||% list(),
      function(row) as.character(row$feature_id)[1L], character(1))
    features <- unique(intersect(features, universe))
    features <- features[nzchar(features)]
    if (length(features) >= min_set_size &&
        length(features) <= max_set_size) members[[name]] <- features
  }
  .trajectory_assert(length(members) > 0L,
    "trajectory_geneset_dynamics: no feature sets pass size filters.")
  list(names = names_all, members = members)
}
