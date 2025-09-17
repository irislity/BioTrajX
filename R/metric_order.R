#' Compute Order consistency (O) via signed adjacent-step concordance
#'
#' @param expr numeric matrix (genes x cells), e.g. log/normalized expression.
#' @param naive_markers character vector of naive/root marker gene names.
#' @param terminal_markers character vector of terminal marker gene names.
#' @param pseudotime numeric vector of length ncol(expr); not necessarily scaled.
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
                      iqr_quantile = 0.5,
                      orientation_max = TRUE) {

  # --- basic checks
  stopifnot(is.matrix(expr) || inherits(expr, "Matrix"))
  stopifnot(is.numeric(pseudotime), length(pseudotime) == ncol(expr))
  if (!is.numeric(iqr_quantile) || length(iqr_quantile) != 1L || iqr_quantile < 0 || iqr_quantile > 1) {
    stop("iqr_quantile must be a number in [0,1].")
  }

  # expected sign per gene
  delta_all <- c(
    stats::setNames(rep(-1L, length(naive_markers)),    naive_markers),
    stats::setNames(rep(+1L, length(terminal_markers)), terminal_markers)
  )

  # de-dup rownames
  rn <- rownames(expr)
  if (anyDuplicated(rn)) {
    keep <- !duplicated(rn)
    expr <- expr[keep, , drop = FALSE]
  }

  genes_all <- intersect(names(delta_all), rownames(expr))
  if (!length(genes_all)) {
    warning("No overlap between markers and expr rownames.")
    return(list(O = NA_real_, O_f = setNames(numeric(0), character(0)),
                genes_used = character(0), n_cells = ncol(expr), n_genes = 0L))
  }

  # per-gene isotonic R^2 in expected direction
  O_once <- function(t) {
    ord <- order(t)
    X   <- as.matrix(expr[genes_all, ord, drop = FALSE])
    t   <- as.numeric(t[ord])
    N   <- ncol(X)
    if (N < 3L) {
      return(list(O = NA_real_, O_f = setNames(numeric(0), character(0)),
                  genes_used = character(0), n_cells = N, n_genes = 0L))
    }

    # IQR filter
    iqr_vals <- apply(X, 1L, stats::IQR)
    thr <- stats::quantile(iqr_vals, probs = iqr_quantile, na.rm = TRUE, names = FALSE)
    keep <- iqr_vals >= thr
    if (!any(keep)) {
      return(list(O = NA_real_, O_f = setNames(numeric(0), character(0)),
                  genes_used = character(0), n_cells = N, n_genes = 0L))
    }
    X   <- X[keep, , drop = FALSE]
    dlt <- delta_all[rownames(X)]  # scalar per row after subsetting

    # isotonic with direction constraint (expects scalar delta)
    iso_R2_dir <- function(y, x, delta_scalar) {
      # enforce expected direction
      yy  <- if (delta_scalar > 0) y else -y
      fit <- stats::isoreg(x, yy)
      # interpolate fitted values to original x
      yhat <- stats::approx(fit$x, fit$yf, xout = x, ties = "ordered")$y
      if (delta_scalar < 0) yhat <- -yhat
      ss_tot <- sum((y - mean(y))^2)
      if (ss_tot <= .Machine$double.eps) return(0)
      r2 <- 1 - sum((y - yhat)^2) / ss_tot
      max(0, min(1, r2))
    }

    # compute per-gene R^2, matching the scalar delta for each row
    Of <- vapply(seq_len(nrow(X)), function(i) {
      iso_R2_dir(y = X[i, ], x = t, delta_scalar = dlt[i])
    }, numeric(1))
    names(Of) <- rownames(X)

    list(O = mean(Of, na.rm = TRUE),
         O_f = Of,
         genes_used = names(Of),
         n_cells = N,
         n_genes = length(Of))
  }

  res_pos <- O_once(pseudotime)
  if (isTRUE(orientation_max)) {
    res_neg <- O_once(-pseudotime)
    if (is.na(res_pos$O) || (!is.na(res_neg$O) && res_neg$O > res_pos$O)) {
      res_neg$orientation <- "-"
      return(res_neg)
    }
  }
  res_pos$orientation <- "+"
  res_pos
}

