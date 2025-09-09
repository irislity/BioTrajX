#' Compute Order consistency (O) via signed adjacent-step concordance
#'
#' @param expr numeric matrix (genes x cells), e.g. log/normalized expression.
#' @param naive_markers character vector of naive/root marker gene names.
#' @param terminal_markers character vector of terminal marker gene names.
#' @param pseudotime numeric vector of length ncol(expr); not necessarily scaled.
#' @param tol numeric tolerance to treat a change as a tie (default 1e-8).
#'
#' @return A list with:
#'   \item{O}{scalar order-consistency score in [0,1]}
#'   \item{O_f}{named per-gene scores in [0,1]}
#'   \item{genes_used}{character vector of genes actually used}
#'   \item{n_cells}{number of cells}
#'   \item{n_genes}{number of genes used}
#' @export
metrics_O <- function(expr,
                      naive_markers,
                      terminal_markers,
                      pseudotime,
                      tol = 1e-8) {

  stopifnot(is.matrix(expr) || inherits(expr, "Matrix"))
  stopifnot(is.numeric(pseudotime), length(pseudotime) == ncol(expr))

  # helper: map expected direction per gene (-1 for naive, +1 for terminal)
  make_delta <- function(naive_markers, terminal_markers) {
    c(
      stats::setNames(rep(-1, length(naive_markers)),   naive_markers),
      stats::setNames(rep(+1, length(terminal_markers)), terminal_markers)
    )
  }

  delta_all <- make_delta(naive_markers, terminal_markers)

  # Intersect with available genes and align
  genes <- intersect(names(delta_all), rownames(expr))
  if (length(genes) == 0L) {
    warning("No overlap between provided markers and rownames(expr). Returning NA.")
    return(list(
      O = NA_real_, O_f = setNames(numeric(0), character(0)),
      genes_used = character(0), n_cells = ncol(expr), n_genes = 0L
    ))
  }
  dlt <- delta_all[genes]

  # Sort cells by pseudotime
  ord   <- order(pseudotime)
  t_ord <- as.numeric(pseudotime[ord])
  N     <- length(t_ord)
  if (N < 2L) {
    warning("Need at least 2 cells to compute adjacent-step concordance. Returning NA.")
    return(list(
      O = NA_real_, O_f = setNames(numeric(0), character(0)),
      genes_used = genes, n_cells = N, n_genes = length(genes)
    ))
  }

  # Reorder expression by pseudotime (genes x cells)
  X_ord <- expr[genes, ord, drop = FALSE]

  # Adjacent differences along pseudotime for each gene (genes x (N-1))
  dx_adj <- X_ord[, -1, drop = FALSE] - X_ord[, -ncol(X_ord), drop = FALSE]

  # Apply expected sign per gene: positive means concordant with biological order
  dx_signed <- sweep(dx_adj, 1, dlt, `*`)

  # Concordant (> tol), ties (<= tol in magnitude)
  concord <- dx_signed > tol
  ties    <- abs(dx_adj) <= tol

  # Per-gene O: fraction of concordant steps + 0.5 * ties
  O_f <- (rowSums(concord) + 0.5 * rowSums(ties)) / (ncol(X_ord) - 1)
  names(O_f) <- rownames(X_ord)

  # Aggregate O across genes
  O <- mean(O_f, na.rm = TRUE)

  list(
    O = O,
    O_f = O_f,
    genes_used = genes,
    n_cells = N,
    n_genes = length(genes)
  )
}
