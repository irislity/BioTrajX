#' Compute Directionality (D) metric
#'
#' @param expr numeric matrix (genes x cells), e.g. log-normalized expression
#' @param early_markers character vector of early/root marker genes
#' @param terminal_markers character vector of terminal marker genes
#' @param pseudotime numeric vector of pseudotime values (length = ncol(expr))
#'
#' @return A list with:
#'   \item{D_term}{terminal directionality score in [0,1]}
#'   \item{D_early}{early directionality score in [0,1]}
#'   \item{s_i_early}{per-cell early program scores}
#'   \item{s_i_term}{per-cell terminal program scores}
#' @export
metrics_d <- function(expr, early_markers, terminal_markers, pseudotime) {

  # helper: mean expression of marker set per cell
  per_cell_score <- function(expr, markers) {
    genes <- intersect(markers, rownames(expr))
    if (length(genes) == 0) {
      warning("No markers found in expression matrix.")
      return(rep(NA_real_, ncol(expr)))
    }
    mark_expr <- expr[genes, , drop = FALSE]
    colSums(mark_expr) / length(genes)
  }

  # compute per-cell scores
  s_i_early <- per_cell_score(expr, early_markers)
  s_i_term  <- per_cell_score(expr, terminal_markers)

  # compute Spearman correlations with pseudotime
  rho_early <- suppressWarnings(cor(pseudotime, s_i_early, method = "spearman", use = "pairwise.complete.obs"))
  rho_term  <- suppressWarnings(cor(pseudotime, s_i_term,  method = "spearman", use = "pairwise.complete.obs"))

  # ReLU: reward only the correct direction, clip wrong-direction to 0
  D_early <- max(0, -rho_early)
  D_term  <- max(0,  rho_term)

  list(
    D_term      = D_term,
    D_early     = D_early,
    s_i_early   = s_i_early,
    s_i_term    = s_i_term
  )
}



#' Regression plot for Directionality (D)
#'
#' @param pseudotime numeric vector
#' @param D_res list returned by metrics_d()
#' @param col_early,col_term colors
#' @param point_alpha transparency of points
#' @param point_cex point size
#' @export
plot_metrics_d <- function(pseudotime,
                               D_res,
                               col_early   = "#1f77b4",
                               col_term    = "#d62728",
                               point_alpha = 0.3,
                               point_cex   = 0.6,
                               main        = "Directionality (D) scores") {

  if (all(is.na(D_res$s_i_early)) || all(is.na(D_res$s_i_term))) {
    warning("Cannot plot D metrics: all scores are NA (no marker overlap).")
    return(invisible(NULL))
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))  # outer margin for title

  # alpha helper
  ac <- function(col, alpha) grDevices::adjustcolor(col, alpha.f = alpha)

  ## Early markers
  rho_early <- suppressWarnings(
    cor(pseudotime, D_res$s_i_early, method = "spearman", use = "pairwise.complete.obs"))
  plot(pseudotime, D_res$s_i_early,
       pch = 16, cex = point_cex, col = ac(col_early, point_alpha),
       xlab = "Pseudotime", ylab = "Early marker score",
       main = bquote(paste("Early  (", rho == .(round(rho_early, 2)),
                           ",  D = ", .(round(D_res$D_early, 2)), ")")))
  abline(lm(D_res$s_i_early ~ pseudotime), col = col_early, lwd = 2)

  ## Terminal markers
  rho_term <- suppressWarnings(
    cor(pseudotime, D_res$s_i_term, method = "spearman", use = "pairwise.complete.obs"))
  plot(pseudotime, D_res$s_i_term,
       pch = 16, cex = point_cex, col = ac(col_term, point_alpha),
       xlab = "Pseudotime", ylab = "Terminal marker score",
       main = bquote(paste("Terminal  (", rho == .(round(rho_term, 2)),
                           ",  D = ", .(round(D_res$D_term, 2)), ")")))
  abline(lm(D_res$s_i_term ~ pseudotime), col = col_term, lwd = 2)

  ## Global title
  mtext(main, outer = TRUE, cex = 1.2, line = 0)
}
