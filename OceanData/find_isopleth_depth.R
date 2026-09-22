# The purpose of this script is to:
# translate OceanData/find_isopleth_depth.m into R.

# Find the depth immediately above the first model level below a target isopleth.
#
# `property_3d_matrix` must have dimensions depth x latitude x longitude.
# The first level below `target_isopleth` is identified along each profile from
# the shallowest model level through the model depth nearest to `depth_max`.
# The returned value is the depth immediately above that first crossing in this
# shallow-to-deep search order. Monotonic depth vectors ordered either
# shallow-to-deep or deep-to-shallow are supported.
#
# `depth_res_interp` is accepted for API compatibility with the Matlab helper
# but is currently unused because the active Matlab logic works directly on the
# native model levels.
find_isopleth_depth <- function(property_3d_matrix,
                                target_isopleth,
                                depth_max,
                                depth,
                                depth_res_interp = NULL) {
  dims <- dim(property_3d_matrix)

  if (length(dims) != 3L) {
    stop("property_3d_matrix must be a 3D array with dimensions depth x latitude x longitude.")
  }

  if (length(depth) != dims[1]) {
    stop("depth must have the same length as the first dimension of property_3d_matrix.")
  }

  if (!is.numeric(depth) || !is.numeric(target_isopleth) || !is.numeric(depth_max)) {
    stop("depth, target_isopleth, and depth_max must be numeric.")
  }

  depth_max_loc <- which.min(abs(depth - depth_max))
  depth_steps <- diff(depth)

  if (all(depth_steps >= 0)) {
    search_idx <- seq_len(depth_max_loc)
  } else if (all(depth_steps <= 0)) {
    search_idx <- seq.int(length(depth), depth_max_loc, by = -1L)
  } else {
    stop("depth must be monotonic so the function can identify shallower and deeper levels.")
  }

  search_depth <- depth[search_idx]
  search_profiles <- property_3d_matrix[search_idx, , , drop = FALSE]
  search_matrix <- matrix(search_profiles, nrow = length(search_idx))
  full_matrix <- matrix(property_3d_matrix, nrow = dims[1])

  below_target <- search_matrix < target_isopleth
  below_target[is.na(below_target)] <- FALSE

  any_crossing <- colSums(below_target) > 0L
  all_missing <- colSums(!is.na(full_matrix)) == 0L

  isopleth_depth <- rep(depth_max, ncol(search_matrix))
  isopleth_depth[all_missing] <- NA_real_

  if (any(any_crossing)) {
    first_crossing <- max.col(t(below_target[, any_crossing, drop = FALSE]), ties.method = "first")
    crossing_depths <- rep(NA_real_, sum(any_crossing))
    valid_crossing <- first_crossing > 1L
    crossing_depths[valid_crossing] <- search_depth[first_crossing[valid_crossing] - 1L]
    isopleth_depth[any_crossing] <- crossing_depths
  }

  matrix(isopleth_depth, nrow = dims[2], ncol = dims[3])
}
