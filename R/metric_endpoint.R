#' Endpoints Validity (E) Metrics
#'
#' Calculate endpoints validity metrics for pseudotime analysis
#'
#' @param pseudotime Numeric vector of pseudotime values
#' @param cluster_labels Character vector of cluster labels (optional)
#' @param naive_marker_scores Numeric vector of naive marker scores (optional)
#' @param terminal_marker_scores Numeric vector of terminal marker scores (optional)
#' @param naive_clusters Character vector of cluster names considered as naive (for Case 2)
#' @param terminal_clusters Character vector of cluster names considered as terminal (for Case 2)
#' @param method Character, one of "gmm" (Case 1), "clusters" (Case 2), or "combined" (Case 3)
#' @param plot Logical, whether to generate density plots for GMM
#'
#' @return An S3 object of class "endpoints_validity" containing all results
#'
#' @examples
#' # Case 1: Using GMM on marker scores
#' result <- calculate_endpoints_validity(
#'   pseudotime = pseudotime_vector,
#'   naive_marker_scores = naive_scores,
#'   terminal_marker_scores = terminal_scores,
#'   method = "gmm"
#' )
#'
#' # Case 2: Using cluster labels
#' result <- calculate_endpoints_validity(
#'   pseudotime = pseudotime_vector,
#'   cluster_labels = cluster_vector,
#'   naive_clusters = c("CD8.NaiveLike"),
#'   terminal_clusters = c("CD8.TEX", "CD8.TPEX"),
#'   method = "clusters"
#' )

metrics_E <- function(pseudotime,
                                         cluster_labels = NULL,
                                         naive_marker_scores = NULL,
                                         terminal_marker_scores = NULL,
                                         naive_clusters = NULL,
                                         terminal_clusters = NULL,
                                         method = c("gmm", "clusters", "combined"),
                                         plot = TRUE) {

  # Input validation
  method <- match.arg(method)
  n_cells <- length(pseudotime)

  if (any(is.na(pseudotime))) {
    stop("Pseudotime vector contains NA values")
  }

  # Initialize result object
  result <- list(
    method = method,
    n_cells = n_cells,
    pseudotime = pseudotime,
    naive_labels = NULL,
    terminal_labels = NULL,
    E_naive = NA,
    E_term = NA,
    E_comp = NA,
    gmm_results = NULL,
    cluster_info = NULL,
    plots = NULL
  )

  # Helper functions
  precision_at_k <- function(order_idx, y, k) {
    if (k <= 0L || sum(!is.na(order_idx)) < k) return(NA_real_)
    mean(y[order_idx[seq_len(k)]] == 1, na.rm = TRUE)
  }

  E_naive_calc <- function(t, y_naive) {
    ord <- order(t, decreasing = FALSE)  # early first
    k <- sum(y_naive == 1, na.rm = TRUE)
    precision_at_k(ord, y_naive, k)
  }

  E_term_calc <- function(t, y_term) {
    ord <- order(t, decreasing = TRUE)   # late first
    k <- sum(y_term == 1, na.rm = TRUE)
    precision_at_k(ord, y_term, k)
  }

  E_norm <- function(E, p) {
    if (is.na(E) || is.na(p) || p == 0) return(NA_real_)
    # Normalize by expected precision under random ordering
    return(E / p)
  }

  create_gmm_labels <- function(scores, label_name) {
    if (!requireNamespace("mclust", quietly = TRUE)) {
      stop("Package 'mclust' is required for GMM method")
    }

    fit <- try(mclust::Mclust(scores, G = 2), silent = TRUE)

    if (inherits(fit, "try-error")) {
      warning(paste("GMM fitting failed for", label_name))
      return(list(labels = rep(NA, length(scores)), fit = NULL, plot_data = NULL))
    }

    # Extract parameters
    mu <- as.numeric(fit$parameters$mean)
    sd <- sqrt(as.numeric(fit$parameters$variance$sigmasq))
    pi <- as.numeric(fit$parameters$pro)
    ord <- order(mu)
    mu <- mu[ord]
    sd <- sd[ord]
    pi <- pi[ord]  # enforce low=1, high=2

    # Get posterior probabilities for high component
    post_high <- predict(fit)$z[, ord[2]]  # P(high component | score)
    labels <- as.integer(post_high >= 0.5)

    # Prepare plot data
    plot_data <- NULL
    if (plot) {
      s <- scores
      ok <- is.finite(s)
      x <- s[ok]
      lab <- labels[ok]

      plot_data <- list(
        x = x,
        labels = lab,
        title = paste("Score density with GMM labels -", label_name),
        fit_params = list(mu = mu, sd = sd, pi = pi)
      )
    }

    return(list(labels = labels, fit = fit, plot_data = plot_data))
  }

  create_density_plot <- function(plot_data) {
    if (is.null(plot_data)) return(NULL)

    x <- plot_data$x
    lab <- plot_data$labels
    d_all <- density(x)
    ylim <- c(0, max(d_all$y) * 2.5)

    plot(density(x), main = plot_data$title, xlab = "Module score",
         ylab = "Density", ylim = ylim)

    if (sum(lab == 1, na.rm = TRUE) >= 2) {
      lines(density(x[lab == 1]), col = "tomato", lwd = 2)
    }
    if (sum(lab == 0, na.rm = TRUE) >= 2) {
      lines(density(x[lab == 0]), col = "steelblue", lwd = 2, lty = 2)
    }

    rug(x[lab == 1], col = "tomato", ticksize = 0.03)
    rug(x[lab == 0], col = "steelblue", ticksize = 0.02)

    legend("topright",
           c("All", "Label=1 (high comp)", "Label=0 (low comp)"),
           lwd = c(1, 2, 2), lty = c(1, 1, 2),
           col = c("black", "tomato", "steelblue"), bty = "n")
  }

  # Method-specific processing
  if (method == "gmm") {
    # Case 1: GMM on marker scores
    if (is.null(naive_marker_scores) || is.null(terminal_marker_scores)) {
      stop("For GMM method, both naive_marker_scores and terminal_marker_scores are required")
    }

    if (length(naive_marker_scores) != n_cells || length(terminal_marker_scores) != n_cells) {
      stop("Marker scores must have the same length as pseudotime")
    }

    # Fit GMM for naive markers
    naive_gmm <- create_gmm_labels(naive_marker_scores, "Naive markers")
    terminal_gmm <- create_gmm_labels(terminal_marker_scores, "Terminal markers")

    result$naive_labels <- naive_gmm$labels
    result$terminal_labels <- terminal_gmm$labels
    result$gmm_results <- list(naive = naive_gmm$fit, terminal = terminal_gmm$fit)

    if (plot) {
      par(mfrow = c(1, 2))
      create_density_plot(naive_gmm$plot_data)
      create_density_plot(terminal_gmm$plot_data)
      par(mfrow = c(1, 1))
    }

  } else if (method == "clusters") {
    # Case 2: Using cluster labels
    if (is.null(cluster_labels)) {
      stop("For clusters method, cluster_labels is required")
    }

    if (is.null(naive_clusters) || is.null(terminal_clusters)) {
      stop("For clusters method, both naive_clusters and terminal_clusters are required")
    }

    if (length(cluster_labels) != n_cells) {
      stop("Cluster labels must have the same length as pseudotime")
    }

    # Create binary labels
    result$naive_labels <- as.integer(cluster_labels %in% naive_clusters)
    result$terminal_labels <- as.integer(cluster_labels %in% terminal_clusters)

    result$cluster_info <- list(
      unique_clusters = unique(cluster_labels),
      naive_clusters = naive_clusters,
      terminal_clusters = terminal_clusters,
      naive_count = sum(result$naive_labels),
      terminal_count = sum(result$terminal_labels)
    )

  } else if (method == "combined") {
    # Case 3: Combined approach
    stop("Combined method not yet implemented. Please use 'gmm' or 'clusters' method.")
  }

  # Calculate metrics
  if (!is.null(result$naive_labels) && !is.null(result$terminal_labels)) {
    result$E_naive <- E_naive_calc(pseudotime, result$naive_labels)
    result$E_term <- E_term_calc(pseudotime, result$terminal_labels)

    # Calculate composite score
    if (!is.na(result$E_naive) && !is.na(result$E_term)) {
      pr <- mean(result$naive_labels == 1, na.rm = TRUE)
      pt <- mean(result$terminal_labels == 1, na.rm = TRUE)

      Ern <- E_norm(result$E_naive, pr)
      Etn <- E_norm(result$E_term, pt)

      if (!any(is.na(c(Ern, Etn))) && !any(c(Ern, Etn) == 0)) {
        result$E_comp <- 2 / (1/Ern + 1/Etn)  # harmonic mean
      }
    }
  }

  # Add summary statistics
  result$summary <- list(
    n_naive = sum(result$naive_labels == 1, na.rm = TRUE),
    n_terminal = sum(result$terminal_labels == 1, na.rm = TRUE),
    prop_naive = mean(result$naive_labels == 1, na.rm = TRUE),
    prop_terminal = mean(result$terminal_labels == 1, na.rm = TRUE),
    E_naive = result$E_naive,
    E_term = result$E_term,
    E_comp = result$E_comp
  )

  # Set class
  class(result) <- "endpoints_validity"

  return(result)
}

#' Print method for endpoints_validity objects
#'
#' @param x An endpoints_validity object
#' @param ... Additional arguments (not used)
#'
#' @export
print.endpoints_validity <- function(x, ...) {
  cat("Endpoints Validity Analysis\n")
  cat("==========================\n\n")

  cat("Method:", x$method, "\n")
  cat("Number of cells:", x$n_cells, "\n\n")

  cat("Cell counts:\n")
  cat("  Naive cells:", x$summary$n_naive, "(", sprintf("%.1f%%", x$summary$prop_naive * 100), ")\n")
  cat("  Terminal cells:", x$summary$n_terminal, "(", sprintf("%.1f%%", x$summary$prop_terminal * 100), ")\n\n")

  cat("Metrics:\n")
  cat("  E_naive: ", sprintf("%.4f", x$E_naive), "\n")
  cat("  E_term:  ", sprintf("%.4f", x$E_term), "\n")
  cat("  E_comp:  ", sprintf("%.4f", x$E_comp), "\n")

  if (x$method == "clusters" && !is.null(x$cluster_info)) {
    cat("\nCluster information:\n")
    cat("  Naive clusters:", paste(x$cluster_info$naive_clusters, collapse = ", "), "\n")
    cat("  Terminal clusters:", paste(x$cluster_info$terminal_clusters, collapse = ", "), "\n")
  }
}

#' Summary method for endpoints_validity objects
#'
#' @param object An endpoints_validity object
#' @param ... Additional arguments (not used)
#'
#' @export
summary.endpoints_validity <- function(object, ...) {
  print(object)

  if (!is.null(object$gmm_results)) {
    cat("\nGMM fitting results:\n")
    if (!is.null(object$gmm_results$naive)) {
      cat("  Naive markers: Successfully fitted 2-component GMM\n")
    }
    if (!is.null(object$gmm_results$terminal)) {
      cat("  Terminal markers: Successfully fitted 2-component GMM\n")
    }
  }

  invisible(object)
}

# Example usage:
#
# # Case 1: Using GMM
# result_gmm <- calculate_endpoints_validity(
#   pseudotime = your_pseudotime,
#   naive_marker_scores = naive_scores,
#   terminal_marker_scores = terminal_scores,
#   method = "gmm",
#   plot = TRUE
# )
#
# # Case 2: Using clusters
# result_clusters <- calculate_endpoints_validity(
#   pseudotime = your_pseudotime,
#   cluster_labels = your_clusters,
#   naive_clusters = c("CD8.NaiveLike"),
#   terminal_clusters = c("CD8.TEX", "CD8.TPEX"),
#   method = "clusters"
# )
#
# print(result_gmm)
# summary(result_clusters)
