# Internal null-coalescing operator used throughout the package
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (isTRUE(all(is.na(x)))) return(y)
  x
}

# Internal: min-max normalize a numeric vector to [0,1], ignoring NAs.
# Returns all zeros (with a warning) when every finite value is identical,
# since (x - min) / (max - min) would otherwise divide by zero.
.minmax_normalize <- function(x, label = "values") {
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)

  if (max_x == min_x) {
    warning("All ", label, " are identical; returning zeros.")
    return(rep(0, length(x)))
  }

  (x - min_x) / (max_x - min_x)
}

#' Reverse pseudotime by min-max normalization
#'
#' This function computes a reversed version of a given pseudotime vector.
#' It first scales the pseudotime values to the range \[0, 1\] using min-max normalization,
#' then reverses the direction by subtracting the scaled values from 1.
#' This is useful when the biological trajectory may run in the opposite direction
#' and both orientations should be considered.
#'
#' @param pseudotime A numeric vector of pseudotime values.
#'
#' @return A numeric vector of reversed pseudotime values, scaled to the range \[0, 1\].
#'
#' @examples
#' pt = pseudotime
#' reverse_pseudotime(pt)
#'
#' @export
reverse_pseudotime <- function(pseudotime) {
  if (!is.numeric(pseudotime)) stop("pseudotime must be numeric.")
  if (all(is.na(pseudotime))) return(rep(NA_real_, length(pseudotime)))

  # Identical values are handled directly (rather than via 1 - .minmax_normalize())
  # so that this degenerate case still returns zeros, not ones.
  if (max(pseudotime, na.rm = TRUE) == min(pseudotime, na.rm = TRUE)) {
    warning("All pseudotime values are identical; returning zeros.")
    return(rep(0, length(pseudotime)))
  }

  1 - .minmax_normalize(pseudotime, label = "pseudotime values")
}
