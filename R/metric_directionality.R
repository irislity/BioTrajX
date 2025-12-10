#' Compute Directionality (D) metric
#'
#' @param expr numeric matrix (genes x cells), e.g. log-normalized expression
#' @param naive_markers character vector of root/naïve marker genes
#' @param terminal_markers character vector of terminal marker genes
#' @param pseudotime numeric vector of pseudotime values (length = ncol(expr))
#'
#' @return A list with:
#'   \item{D_term}{terminal directionality score in [0,1]}
#'   \item{D_naive}{naïve directionality score in [0,1]}
#'   \item{s_i_naive}{per-cell naïve program scores}
#'   \item{s_i_term}{per-cell terminal program scores}
#' @export
metrics_D <- function(expr, naive_markers, terminal_markers, pseudotime) {

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
  s_i_naive <- per_cell_score(expr, naive_markers)
  s_i_term  <- per_cell_score(expr, terminal_markers)

  # compute Spearman correlations with pseudotime
  rho_naive <- suppressWarnings(cor(pseudotime, s_i_naive, method = "spearman", use = "pairwise.complete.obs"))
  rho_term  <- suppressWarnings(cor(pseudotime, s_i_term,  method = "spearman", use = "pairwise.complete.obs"))

  # rescale correlations from [-1,1] to [0,1]
  D_naive <- (1 + rho_naive) / 2
  D_term  <- (1 + rho_term)  / 2

  list(
    D_term      = D_term,
    D_naive     = D_naive,
    s_i_naive   = s_i_naive,
    s_i_term    = s_i_term
  )
}



#' Regression plot for Directionality (D)
#'
#' @param pseudotime numeric vector
#' @param D_res list returned by metrics_D()
#' @param col_naive,col_term colors
#' @param point_alpha transparency of points
#' @param point_cex point size
#' @export
plot_metrics_D <- function(pseudotime,
                               D_res,
                               col_naive   = "#1f77b4",
                               col_term    = "#d62728",
                               point_alpha = 0.3,
                               point_cex   = 0.6,
                               main        = "Directionality (D) scores") {

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))  # outer margin for title

  # alpha helper
  ac <- function(col, alpha) grDevices::adjustcolor(col, alpha.f = alpha)

  ## Naïve markers
  plot(pseudotime, D_res$s_i_naive,
       pch = 16, cex = point_cex, col = ac(col_naive, point_alpha),
       xlab = "Pseudotime", ylab = "Naïve marker score",
       main = sprintf("Naïve (ρ = %.2f, D = %.2f)",
                      suppressWarnings(cor(pseudotime, D_res$s_i_naive, method="spearman", use="pairwise.complete.obs")),
                      D_res$D_naive))
  abline(lm(D_res$s_i_naive ~ pseudotime), col = col_naive, lwd = 2)

  ## Terminal markers
  plot(pseudotime, D_res$s_i_term,
       pch = 16, cex = point_cex, col = ac(col_term, point_alpha),
       xlab = "Pseudotime", ylab = "Terminal marker score",
       main = sprintf("Terminal (ρ = %.2f, D = %.2f)",
                      suppressWarnings(cor(pseudotime, D_res$s_i_term, method="spearman", use="pairwise.complete.obs")),
                      D_res$D_term))
  abline(lm(D_res$s_i_term ~ pseudotime), col = col_term, lwd = 2)

  ## Global title
  mtext(main, outer = TRUE, cex = 1.2, line = 0)
}
