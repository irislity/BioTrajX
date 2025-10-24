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

  min_pt <- min(pseudotime, na.rm = TRUE)
  max_pt <- max(pseudotime, na.rm = TRUE)

  if (max_pt == min_pt) {
    warning("All pseudotime values are identical; returning zeros.")
    return(rep(0, length(pseudotime)))
  }

  scaled <- (pseudotime - min_pt) / (max_pt - min_pt)
  reversed <- 1 - scaled
  return(reversed)
}
