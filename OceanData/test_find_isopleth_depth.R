source("OceanData/find_isopleth_depth.R")

expect_equal_with_na <- function(object, expected, tolerance = sqrt(.Machine$double.eps)) {
  same_na <- identical(is.na(object), is.na(expected))
  same_values <- isTRUE(all.equal(object, expected, tolerance = tolerance, check.attributes = FALSE))

  if (!same_na || !same_values) {
    stop("Objects are not equal.", call. = FALSE)
  }
}

depth <- c(0, 50, 100, 150)

# Normal crossing
normal_crossing <- array(c(120, 110, 90, 80), dim = c(4, 1, 1))
expect_equal_with_na(find_isopleth_depth(normal_crossing, 100, 150, depth), matrix(50, nrow = 1, ncol = 1))

# All-NA profile
all_na <- array(c(NA, NA, NA, NA), dim = c(4, 1, 1))
expect_equal_with_na(find_isopleth_depth(all_na, 100, 150, depth), matrix(NA_real_, nrow = 1, ncol = 1))

# No crossing returns depth_max even when depth_max is between model levels
no_crossing <- array(c(120, 110, 105, 101), dim = c(4, 1, 1))
expect_equal_with_na(find_isopleth_depth(no_crossing, 100, 120, depth), matrix(120, nrow = 1, ncol = 1))

# Crossing at the first level has no depth above it
first_level_crossing <- array(c(90, 80, 70, 60), dim = c(4, 1, 1))
expect_equal_with_na(find_isopleth_depth(first_level_crossing, 100, 150, depth), matrix(NA_real_, nrow = 1, ncol = 1))

# Internal NAs before the first crossing are ignored, consistent with Matlab comparisons
na_before_crossing <- array(c(120, NA, 90, 80), dim = c(4, 1, 1))
expect_equal_with_na(find_isopleth_depth(na_before_crossing, 100, 150, depth), matrix(50, nrow = 1, ncol = 1))

# Descending depth vectors are searched from shallow to deep
depth_desc <- c(150, 100, 50, 0)
descending_profile <- array(c(80, 90, 130, 140), dim = c(4, 1, 1))
expect_equal_with_na(find_isopleth_depth(descending_profile, 100, 150, depth_desc), matrix(50, nrow = 1, ncol = 1))

# Multiple latitude/longitude cells
multi_cell <- array(
  c(
    120, 110, 90, 80,
    90, 80, 70, 60,
    120, 110, 105, 101,
    NA, NA, NA, NA
  ),
  dim = c(4, 2, 2)
)
expected_multi <- matrix(c(50, NA, 120, NA), nrow = 2, ncol = 2)
expect_equal_with_na(find_isopleth_depth(multi_cell, 100, 120, depth), expected_multi)

message("find_isopleth_depth.R validation passed")
