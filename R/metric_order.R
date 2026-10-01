# ── Null R² lookup table (10,000 simulations per K, seed 2024) ───────────────
# Expected isotonic R² for K i.i.d. N(0,1) points (non-decreasing fit).
.O_NULL_R2_TABLE <- c(
  "5"  = 0.31973,
  "6"  = 0.28868,
  "7"  = 0.26882,
  "8"  = 0.24526,
  "9"  = 0.22713,
  "10" = 0.21386,
  "12" = 0.19330,
  "15" = 0.16654,
  "20" = 0.13913,
  "25" = 0.11822,
  "30" = 0.10388,
  "40" = 0.08435,
  "50" = 0.07109
)

.O_get_null_r2 <- function(n_bins) {
  keys <- as.integer(names(.O_NULL_R2_TABLE))
  vals <- unname(.O_NULL_R2_TABLE)
  if (n_bins %in% keys) return(vals[match(n_bins, keys)])
  stats::approx(log(keys), vals, xout = log(n_bins), rule = 2)$y
}

#' Compute Order consistency (O) via pseudobulk isotonic R² (null-calibrated)
#'
#' Cells are binned into \code{n_bins} equal-width pseudotime intervals and
#' expression is averaged within each bin (pseudobulk).  Direction-constrained
#' isotonic regression is then fitted to the \code{n_bins} bin means, and the
#' resulting R² is rescaled by the theoretical null:
#'
#' \deqn{O = \max\!\left(0,\; \frac{R^2 - R^2_{\text{null}}}{1 - R^2_{\text{null}}}\right)}
#'
#' \eqn{R^2_{\text{null}}} is the expected isotonic R² for \code{n_bins}
#' i.i.d. data points (precomputed from 10,000 simulations; log-linearly
#' interpolated for off-table values).  This gives O \eqn{\approx} 0 for
#' flat/random data and O \eqn{\approx} 1 for a strong, clean monotone trend.
#'
#' @param expr numeric matrix (genes x cells), e.g. log-normalised expression.
#' @param early_markers character vector of early/root marker gene names.
#' @param terminal_markers character vector of terminal marker gene names.
#' @param pseudotime numeric vector of length \code{ncol(expr)}.
#' @param iqr_quantile numeric in \[0,1\]; genes whose IQR falls below this
#'   quantile are excluded before scoring (default 0.5).
#' @param orientation_max logical; if \code{TRUE} (default) both pseudotime
#'   orientations are evaluated and the higher O is returned.
#' @param n_bins integer \eqn{\ge 3}; number of pseudotime bins for pseudobulk
#'   averaging (default 10).  More cells per bin yields cleaner bin means and
#'   higher O for a true signal.
#'
#' @return A list with:
#'   \item{O}{scalar order-consistency score, null-calibrated to \[0,1\]}
#'   \item{O_f}{named per-gene calibrated scores}
#'   \item{genes_used}{character vector of genes used}
#'   \item{n_cells}{number of cells}
#'   \item{n_genes}{number of genes used}
#'   \item{n_bins}{number of bins actually used}
#'   \item{null_r2}{the \eqn{R^2_{\text{null}}} value used for rescaling}
#'   \item{orientation}{"+" or "-"}
#' @export

metrics_o <- function(expr,
                      early_markers,
                      terminal_markers,
                      pseudotime,
                      iqr_quantile    = 0.5,
                      orientation_max = TRUE,
                      n_bins          = 10L) {

  stopifnot(is.matrix(expr) || inherits(expr, "Matrix"))
  stopifnot(is.numeric(pseudotime), length(pseudotime) == ncol(expr))
  if (!is.numeric(iqr_quantile) || length(iqr_quantile) != 1L ||
      iqr_quantile < 0 || iqr_quantile > 1)
    stop("iqr_quantile must be a number in [0,1].")

  n_bins  <- max(3L, as.integer(n_bins))
  null_r2 <- .O_get_null_r2(n_bins)

  delta_all <- c(
    stats::setNames(rep(-1L, length(early_markers)),    early_markers),
    stats::setNames(rep(+1L, length(terminal_markers)), terminal_markers)
  )

  rn <- rownames(expr)
  if (anyDuplicated(rn))
    expr <- expr[!duplicated(rn), , drop = FALSE]

  genes_all <- intersect(names(delta_all), rownames(expr))
  if (!length(genes_all)) {
    warning("No overlap between markers and expr rownames.")
    return(list(O = NA_real_, O_f = setNames(numeric(0), character(0)),
                genes_used = character(0), n_cells = ncol(expr), n_genes = 0L,
                n_bins = n_bins, null_r2 = null_r2))
  }

  # ── bin cells by pseudotime, compute bin means ────────────────────────────
  pseudobulk <- function(t) {
    breaks  <- seq(min(t), max(t), length.out = n_bins + 1L)
    bins    <- cut(t, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    bin_ids <- sort(unique(bins))
    pb_mat  <- vapply(bin_ids, function(b) {
      idx <- which(bins == b)
      if (length(idx) == 1L) as.numeric(expr[genes_all, idx, drop = TRUE])
      else rowMeans(expr[genes_all, idx, drop = FALSE])
    }, numeric(length(genes_all)))
    rownames(pb_mat) <- genes_all
    bin_pt <- vapply(bin_ids, function(b) mean(t[which(bins == b)]), numeric(1))
    list(mat = pb_mat, pt = bin_pt, n_bins_used = length(bin_ids))
  }

  # ── direction-constrained isotonic R² then null-calibrate ─────────────────
  iso_R2_calibrated <- function(y, x, delta_scalar) {
    yy   <- if (delta_scalar > 0) y else -y
    fit  <- stats::isoreg(x, yy)
    yhat <- stats::approx(fit$x, fit$yf, xout = x, ties = "ordered")$y
    if (delta_scalar < 0) yhat <- -yhat
    ss_tot <- sum((y - mean(y))^2)
    if (ss_tot <= .Machine$double.eps) return(0)
    r2_raw <- max(0, min(1, 1 - sum((y - yhat)^2) / ss_tot))
    max(0, (r2_raw - null_r2) / (1 - null_r2))
  }

  O_once <- function(t) {
    N <- length(t)
    if (N < 3L)
      return(list(O = NA_real_, O_f = setNames(numeric(0), character(0)),
                  genes_used = character(0), n_cells = N, n_genes = 0L,
                  n_bins = n_bins, null_r2 = null_r2))

    pb <- pseudobulk(t)
    X  <- pb$mat
    tp <- pb$pt

    # IQR filter on original single-cell data
    iqr_vals <- apply(expr[genes_all, , drop = FALSE], 1L, stats::IQR)
    thr  <- stats::quantile(iqr_vals, probs = iqr_quantile, na.rm = TRUE, names = FALSE)
    keep <- iqr_vals >= thr
    if (!any(keep))
      return(list(O = NA_real_, O_f = setNames(numeric(0), character(0)),
                  genes_used = character(0), n_cells = N, n_genes = 0L,
                  n_bins = n_bins, null_r2 = null_r2))

    X   <- X[keep, , drop = FALSE]
    dlt <- delta_all[rownames(X)]

    Of <- vapply(seq_len(nrow(X)), function(i)
      iso_R2_calibrated(y = X[i, ], x = tp, delta_scalar = dlt[i]),
      numeric(1))
    names(Of) <- rownames(X)

    list(O          = mean(Of, na.rm = TRUE),
         O_f        = Of,
         genes_used = names(Of),
         n_cells    = N,
         n_genes    = length(Of),
         n_bins     = pb$n_bins_used,
         null_r2    = null_r2)
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

#' Plot Order Consistency (O): overview + isotonic fits
#'
#' @param expr numeric matrix (genes x cells), same used in metrics_o()
#' @param pseudotime numeric vector, same used in metrics_o()
#' @param O_res list returned by metrics_o()
#' @param early_markers character vector of early/root genes (expected to decrease)
#' @param terminal_markers character vector of terminal genes (expected to increase)
#' @param top_n integer, number of example genes to plot with isotonic fits
#' @param point_cex numeric, point size for per-cell dots
#' @param point_alpha numeric in \[0,1\], transparency for dots
#' @param line_lwd numeric, line width for isotonic fit
#' @param main title across the whole figure
#'
#' @details
#' Left panel: histogram of per-gene O_f with vertical line at mean O.
#' Right panel: for top_n genes (by O_f), scatter of expression vs (oriented) pseudotime
#' with their direction-constrained isotonic regression fit overlaid.
#'
#' @export
plot_metrics_o <- function(expr,
                           pseudotime,
                           O_res,
                           early_markers,
                           terminal_markers,
                           top_n       = 6,
                           point_cex   = 0.45,
                           point_alpha = 0.35,
                           line_lwd    = 2,
                           main        = "Order Consistency (O) score") {
  stopifnot(is.matrix(expr) || inherits(expr, "Matrix"))
  stopifnot(is.numeric(pseudotime), length(pseudotime) == ncol(expr))
  stopifnot(is.list(O_res), !is.null(O_res$O_f))

  # orientation from metrics_o()
  t_oriented <- if (identical(O_res$orientation, "-")) -pseudotime else pseudotime

  # delta map: -1 for early (expected down), +1 for terminal (expected up)
  delta_map <- c(
    stats::setNames(rep(-1L, length(early_markers)),    early_markers),
    stats::setNames(rep(+1L, length(terminal_markers)), terminal_markers)
  )

  # restrict to genes actually used by metrics_o()
  genes_used <- intersect(names(O_res$O_f), intersect(rownames(expr), names(delta_map)))
  if (!length(genes_used)) {
    warning("No overlap between O_res$O_f, expr rownames, and marker lists.")
    return(invisible(NULL))
  }

  Of <- O_res$O_f[genes_used]
  # pick top_n genes by O_f (highest concordance)
  top_n <- min(top_n, length(Of))
  top_genes <- names(sort(Of, decreasing = TRUE))[seq_len(top_n)]

  # helper: semi-transparent color
  ac <- function(col, alpha) grDevices::adjustcolor(col, alpha.f = alpha)

  # isotonic fit respecting expected direction (same logic as metrics_o)
  iso_fit_vals <- function(y, x, delta_scalar) {
    yy  <- if (delta_scalar > 0) y else -y
    fit <- stats::isoreg(x, yy)
    yhat <- stats::approx(fit$x, fit$yf, xout = x, ties = "ordered")$y
    if (delta_scalar < 0) yhat <- -yhat
    list(x = x, yhat = yhat)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(oma = c(0,0,2,0))


  ## top genes — expression vs oriented pseudotime + isotonic line
  # Choose a distinct palette (base R) for up to ~8 genes; recycle if needed
  base_cols <- c("#1f77b4","#d62728","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22")
  cols <- rep(base_cols, length.out = top_n)

  # set up empty plot
  rng_x <- range(t_oriented, na.rm = TRUE)
  # compute y-range across selected genes
  y_min <- Inf; y_max <- -Inf
  for (g in top_genes) {
    yy <- expr[g, , drop = TRUE]
    y_min <- min(y_min, min(yy, na.rm = TRUE))
    y_max <- max(y_max, max(yy, na.rm = TRUE))
  }
  plot(NA, NA, xlim = rng_x, ylim = c(y_min, y_max),
       xlab = "Oriented pseudotime", ylab = "Expression")

  # add points + fits gene-by-gene
  ord <- order(t_oriented)
  t_o <- t_oriented[ord]
  i <- 0
  for (g in top_genes) {
    i <- i + 1
    y   <- as.numeric(expr[g, ])
    y_o <- y[ord]
    dlt <- delta_map[[g]]
    fit <- iso_fit_vals(y_o, t_o, dlt)

    points(t_o, y_o, pch = 16, cex = point_cex, col = ac(cols[i], point_alpha))
    # draw the piecewise-constant isotonic fit as connected segments
    # (approx already returned yhat at each x; drawing as lines is fine)
    lines(fit$x, fit$yhat, col = cols[i], lwd = line_lwd)
  }
  legend("topleft",
         legend = paste0(seq_along(top_genes), ". ", top_genes,
                         " (O_f=", sprintf("%.2f", Of[top_genes]), ")"),
         col = cols[seq_along(top_genes)], lwd = 2, cex = 0.8, bty = "n")

  legend("topright",
         legend = sprintf("O = %.3f  (orientation %s)", O_res$O, ifelse(identical(O_res$orientation, "-"), "−", "+")),
         bty = "n")

  ## ---- Global title
  mtext(main, outer = TRUE, cex = 1.2, line = 0)
}


