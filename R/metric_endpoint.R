#' Endpoints Validity (E) Metrics with Robust GMM Handling
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
#' @param min_cells_per_component Minimum cells required per GMM component (default: 10)
#' @param gmm_retry_methods Vector of fallback methods if GMM fails
#' @param verbose Logical, whether to print diagnostic messages
#'
#' @return An S3 object of class "endpoints_validity" containing all results

metrics_E <- function(pseudotime,
                      cluster_labels = NULL,
                      naive_marker_scores = NULL,
                      terminal_marker_scores = NULL,
                      naive_clusters = NULL,
                      terminal_clusters = NULL,
                      method = c("gmm", "clusters", "combined"),
                      plot = TRUE,
                      min_cells_per_component = 10,
                      gmm_retry_methods = c("quantile", "kmeans", "manual"),
                      verbose = FALSE) {

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
    plots = NULL,
    warnings = character(0),
    fallback_used = NULL
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

  # Robust GMM with multiple fallback strategies
  create_robust_gmm_labels <- function(scores, label_name) {
    if (length(scores) < min_cells_per_component * 2) {
      warning(paste("Too few cells for GMM fitting:", label_name, "- need at least", min_cells_per_component * 2))
      return(list(labels = rep(NA, length(scores)), fit = NULL, plot_data = NULL, method_used = "insufficient_data"))
    }

    # Remove NAs and check for variation
    finite_scores <- scores[is.finite(scores)]
    if (length(finite_scores) < min_cells_per_component * 2) {
      warning(paste("Too few finite scores for GMM fitting:", label_name))
      return(list(labels = rep(NA, length(scores)), fit = NULL, plot_data = NULL, method_used = "insufficient_finite_data"))
    }

    if (var(finite_scores) < 1e-10) {
      warning(paste("Insufficient variation in scores for GMM fitting:", label_name))
      return(list(labels = rep(NA, length(scores)), fit = NULL, plot_data = NULL, method_used = "no_variation"))
    }

    # Try standard GMM first
    gmm_result <- try_standard_gmm(scores, label_name)
    if (!is.null(gmm_result$labels) && !all(is.na(gmm_result$labels))) {
      return(c(gmm_result, list(method_used = "mclust_gmm")))
    }

    # Try fallback methods
    for (fallback_method in gmm_retry_methods) {
      if (verbose) message(paste("Trying fallback method:", fallback_method, "for", label_name))

      fallback_result <- switch(fallback_method,
                                "quantile" = try_quantile_split(scores, label_name),
                                "kmeans" = try_kmeans_split(scores, label_name),
                                "manual" = try_manual_split(scores, label_name),
                                NULL
      )

      if (!is.null(fallback_result$labels) && !all(is.na(fallback_result$labels))) {
        if (verbose) message(paste("Successfully used fallback method:", fallback_method, "for", label_name))
        return(c(fallback_result, list(method_used = fallback_method)))
      }
    }

    warning(paste("All methods failed for", label_name))
    return(list(labels = rep(NA, length(scores)), fit = NULL, plot_data = NULL, method_used = "all_failed"))
  }

  # Standard mclust GMM
  try_standard_gmm <- function(scores, label_name) {
    if (!requireNamespace("mclust", quietly = TRUE)) {
      warning("Package 'mclust' is required for GMM method")
      return(list(labels = rep(NA, length(scores)), fit = NULL, plot_data = NULL))
    }

    # Try different model types if default fails
    model_types <- c("E", "V", "EII", "VII")

    for (model_type in model_types) {
      fit <- try({
        mclust::Mclust(scores, G = 2, modelNames = model_type, verbose = FALSE)
      }, silent = TRUE)

      if (!inherits(fit, "try-error") && !is.null(fit)) {
        if (verbose && model_type != "E") {
          message(paste("GMM succeeded with model type:", model_type, "for", label_name))
        }

        # Validate the fit
        if (is.null(fit$parameters) || is.null(fit$classification)) {
          next
        }

        # Check component sizes
        component_sizes <- table(fit$classification)
        if (any(component_sizes < min_cells_per_component)) {
          if (verbose) message(paste("Components too small for", label_name, "- trying next model"))
          next
        }

        # Extract parameters safely
        mu <- tryCatch({
          as.numeric(fit$parameters$mean)
        }, error = function(e) NULL)

        if (is.null(mu) || length(mu) != 2) {
          next
        }

        # Order components by mean
        ord <- order(mu)

        # Get posterior probabilities for high component
        post_probs <- tryCatch({
          predict(fit)$z
        }, error = function(e) {
          # Fallback to classification
          matrix(0, nrow = length(scores), ncol = 2)
        })

        if (ncol(post_probs) >= 2) {
          post_high <- post_probs[, ord[2]]  # P(high component | score)
          labels <- as.integer(post_high >= 0.5)
        } else {
          # Use classification directly
          labels <- as.integer(fit$classification == ord[2])
        }

        # Prepare plot data
        plot_data <- NULL
        if (plot) {
          s <- scores
          ok <- is.finite(s)
          x <- s[ok]
          lab <- labels[ok]

          if (length(mu) == 2) {
            # Extract variance safely
            variance <- tryCatch({
              if (is.matrix(fit$parameters$variance$sigmasq)) {
                diag(fit$parameters$variance$sigmasq)
              } else {
                as.numeric(fit$parameters$variance$sigmasq)
              }
            }, error = function(e) rep(1, 2))

            sd <- sqrt(variance)
            pi <- as.numeric(fit$parameters$pro)

            if (length(sd) == 2 && length(pi) == 2) {
              plot_data <- list(
                x = x,
                labels = lab,
                title = paste("Score density with GMM labels -", label_name),
                fit_params = list(mu = mu[ord], sd = sd[ord], pi = pi[ord])
              )
            }
          }
        }

        return(list(labels = labels, fit = fit, plot_data = plot_data))
      }
    }

    return(list(labels = rep(NA, length(scores)), fit = NULL, plot_data = NULL))
  }

  # Fallback: Quantile-based split
  try_quantile_split <- function(scores, label_name) {
    finite_scores <- scores[is.finite(scores)]
    if (length(finite_scores) < 10) return(NULL)

    # Use 75th percentile as threshold (or median if too few high scorers)
    threshold <- quantile(finite_scores, 0.75, na.rm = TRUE)
    high_count <- sum(finite_scores > threshold)

    # Adjust threshold if too few cells in high group
    if (high_count < min_cells_per_component) {
      threshold <- quantile(finite_scores, 0.5, na.rm = TRUE)
    }

    labels <- as.integer(scores > threshold)
    labels[!is.finite(scores)] <- NA

    plot_data <- if (plot) {
      list(x = finite_scores, labels = labels[is.finite(scores)],
           title = paste("Quantile split -", label_name))
    } else NULL

    return(list(labels = labels, fit = NULL, plot_data = plot_data))
  }

  # Fallback: K-means clustering
  try_kmeans_split <- function(scores, label_name) {
    finite_scores <- scores[is.finite(scores)]
    if (length(finite_scores) < 10) return(NULL)

    kmeans_result <- try({
      kmeans(finite_scores, centers = 2, nstart = 10)
    }, silent = TRUE)

    if (inherits(kmeans_result, "try-error")) return(NULL)

    # Order clusters by center value
    center_order <- order(kmeans_result$centers)
    labels <- rep(NA_integer_, length(scores))
    labels[is.finite(scores)] <- as.integer(kmeans_result$cluster == center_order[2])

    plot_data <- if (plot) {
      list(x = finite_scores, labels = labels[is.finite(scores)],
           title = paste("K-means split -", label_name))
    } else NULL

    return(list(labels = labels, fit = NULL, plot_data = plot_data))
  }

  # Fallback: Manual threshold based on distribution
  try_manual_split <- function(scores, label_name) {
    finite_scores <- scores[is.finite(scores)]
    if (length(finite_scores) < 10) return(NULL)

    # Use mean + 0.5*sd as threshold
    threshold <- mean(finite_scores) + 0.5 * sd(finite_scores)
    labels <- as.integer(scores > threshold)
    labels[!is.finite(scores)] <- NA

    plot_data <- if (plot) {
      list(x = finite_scores, labels = labels[is.finite(scores)],
           title = paste("Manual threshold split -", label_name))
    } else NULL

    return(list(labels = labels, fit = NULL, plot_data = plot_data))
  }

  # Enhanced plotting function
  create_density_plot <- function(plot_data, method_used = "GMM") {
    if (is.null(plot_data)) return(NULL)

    x <- plot_data$x
    lab <- plot_data$labels

    if (length(x) < 2 || all(is.na(lab))) return(NULL)

    d_all <- density(x, na.rm = TRUE)
    ylim <- c(0, max(d_all$y, na.rm = TRUE) * 1.2)

    plot(d_all, main = paste(plot_data$title, "(", method_used, ")"),
         xlab = "Module score", ylab = "Density", ylim = ylim)

    # Plot class-specific densities if enough data
    if (sum(lab == 1, na.rm = TRUE) >= 2) {
      lines(density(x[lab == 1 & !is.na(lab)], na.rm = TRUE), col = "tomato", lwd = 2)
    }
    if (sum(lab == 0, na.rm = TRUE) >= 2) {
      lines(density(x[lab == 0 & !is.na(lab)], na.rm = TRUE), col = "steelblue", lwd = 2, lty = 2)
    }

    # Add rug plots
    rug(x[lab == 1 & !is.na(lab)], col = "tomato", ticksize = 0.03)
    rug(x[lab == 0 & !is.na(lab)], col = "steelblue", ticksize = 0.02)

    legend("topright",
           c("All", "Label=1 (high)", "Label=0 (low)"),
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

    # Fit robust GMM for naive markers
    naive_gmm <- create_robust_gmm_labels(naive_marker_scores, "Naive markers")
    terminal_gmm <- create_robust_gmm_labels(terminal_marker_scores, "Terminal markers")

    result$naive_labels <- naive_gmm$labels
    result$terminal_labels <- terminal_gmm$labels
    result$gmm_results <- list(naive = naive_gmm$fit, terminal = terminal_gmm$fit)
    result$fallback_used <- list(naive = naive_gmm$method_used, terminal = terminal_gmm$method_used)

    if (plot) {
      par(mfrow = c(1, 2))
      create_density_plot(naive_gmm$plot_data, naive_gmm$method_used)
      create_density_plot(terminal_gmm$plot_data, terminal_gmm$method_used)
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
    # Case 3: Combined approach - try GMM first, fall back to clusters
    if (verbose) message("Using combined method: trying GMM first, then clusters as fallback")

    gmm_success <- FALSE

    if (!is.null(naive_marker_scores) && !is.null(terminal_marker_scores)) {
      naive_gmm <- create_robust_gmm_labels(naive_marker_scores, "Naive markers")
      terminal_gmm <- create_robust_gmm_labels(terminal_marker_scores, "Terminal markers")

      if (!all(is.na(naive_gmm$labels)) && !all(is.na(terminal_gmm$labels))) {
        result$naive_labels <- naive_gmm$labels
        result$terminal_labels <- terminal_gmm$labels
        result$gmm_results <- list(naive = naive_gmm$fit, terminal = terminal_gmm$fit)
        result$fallback_used <- list(naive = naive_gmm$method_used, terminal = terminal_gmm$method_used)
        gmm_success <- TRUE

        if (plot) {
          par(mfrow = c(1, 2))
          create_density_plot(naive_gmm$plot_data, naive_gmm$method_used)
          create_density_plot(terminal_gmm$plot_data, terminal_gmm$method_used)
          par(mfrow = c(1, 1))
        }
      }
    }

    # Fall back to clusters if GMM failed
    if (!gmm_success) {
      if (verbose) message("GMM failed, falling back to cluster method")

      if (!is.null(cluster_labels) && !is.null(naive_clusters) && !is.null(terminal_clusters)) {
        result$naive_labels <- as.integer(cluster_labels %in% naive_clusters)
        result$terminal_labels <- as.integer(cluster_labels %in% terminal_clusters)
        result$fallback_used <- list(method = "clusters")

        result$cluster_info <- list(
          unique_clusters = unique(cluster_labels),
          naive_clusters = naive_clusters,
          terminal_clusters = terminal_clusters,
          naive_count = sum(result$naive_labels),
          terminal_count = sum(result$terminal_labels)
        )
      } else {
        stop("Combined method failed: GMM unsuccessful and cluster information incomplete")
      }
    }
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
print.endpoints_validity <- function(x, ...) {
  cat("Endpoints Validity Analysis\n")
  cat("==========================\n\n")

  cat("Method:", x$method, "\n")
  cat("Number of cells:", x$n_cells, "\n")

  if (!is.null(x$fallback_used)) {
    if (is.list(x$fallback_used) && !is.null(x$fallback_used$naive)) {
      cat("Naive method used:", x$fallback_used$naive, "\n")
      cat("Terminal method used:", x$fallback_used$terminal, "\n")
    } else if (!is.null(x$fallback_used$method)) {
      cat("Fallback method used:", x$fallback_used$method, "\n")
    }
  }
  cat("\n")

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
summary.endpoints_validity <- function(object, ...) {
  print(object)

  if (!is.null(object$gmm_results)) {
    cat("\nGMM fitting results:\n")
    if (!is.null(object$gmm_results$naive)) {
      cat("  Naive markers: Successfully fitted 2-component GMM\n")
    } else if (!is.null(object$fallback_used$naive)) {
      cat("  Naive markers: Used fallback method -", object$fallback_used$naive, "\n")
    }
    if (!is.null(object$gmm_results$terminal)) {
      cat("  Terminal markers: Successfully fitted 2-component GMM\n")
    } else if (!is.null(object$fallback_used$terminal)) {
      cat("  Terminal markers: Used fallback method -", object$fallback_used$terminal, "\n")
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



#' Plot precision@k curves for Endpoints Validity (E)
#'
#' @param E_obj An object returned by metrics_E() (class "endpoints_validity")
#' @param overlay Logical; if TRUE, plot Naive & Terminal on one axis; else 2 panels
#' @param main Global title
#' @param col_naive,col_term Colors for curves/points
#' @param lwd Line width for curves
#' @param point_cex Size for the marker at k = #positives
#' @param point_alpha Transparency for the scatter rug (optional visual cue)
#' @param show_rug Logical; draw a small rug showing positive positions along the ranked list
#'
#' @details
#' For each label vector y in {0,1}, ordered by the appropriate direction of pseudotime,
#' we compute precision@k = (# positives in top k) / k for k = 1..N.
#' The plotted dot highlights precision at k = number of positives (which is what E uses).
#'
#' @export
plot_metrics_E <- function(E_obj,
                           overlay      = TRUE,
                           main         = "Endpoints Validity: precision",
                           col_naive    = "#1f77b4",
                           col_term     = "#d62728",
                           lwd          = 2,
                           point_cex    = 1.1,
                           point_alpha  = 0.7,
                           show_rug     = TRUE) {
  stopifnot(inherits(E_obj, "endpoints_validity"))

  t <- E_obj$pseudotime
  yN <- E_obj$naive_labels
  yT <- E_obj$terminal_labels

  # helper: compute precision@k curve given ordering direction
  precision_curve <- function(t, y, decreasing = FALSE) {
    ok <- is.finite(t) & !is.na(y)
    if (!any(ok)) return(list(k = integer(0), prec = numeric(0),
                              k_star = NA_integer_, prec_star = NA_real_,
                              p_base = NA_real_, ranks_pos = integer(0)))
    t  <- t[ok]; y <- as.integer(y[ok])
    ord <- order(t, decreasing = decreasing)
    y_o <- y[ord]
    cpos <- cumsum(y_o == 1L)
    kvec <- seq_along(y_o)
    prec <- cpos / kvec

    k_star <- sum(y_o == 1L)                      # number of positives
    prec_star <- if (k_star > 0L) mean(y_o[seq_len(k_star)] == 1L) else NA_real_
    p_base <- mean(y_o == 1L)                     # random baseline
    ranks_pos <- which(y_o == 1L)                 # positions of positives

    list(k = kvec, prec = prec,
         k_star = k_star, prec_star = prec_star,
         p_base = p_base, ranks_pos = ranks_pos)
  }

  # compute curves
  curN <- if (!is.null(yN)) precision_curve(t, yN, decreasing = FALSE) else NULL
  curT <- if (!is.null(yT)) precision_curve(t, yT, decreasing = TRUE)  else NULL

  # alpha helper (for rug color)
  ac <- function(col, alpha) grDevices::adjustcolor(col, alpha.f = alpha)

  old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
  if (overlay) {
    par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))

    # set up empty plot bounds
    xmax <- max(c(curN$k, curT$k), na.rm = TRUE)
    ymax <- max(c(curN$prec, curT$prec, curN$p_base, curT$p_base), na.rm = TRUE)
    if (!is.finite(xmax)) { plot.new(); mtext(main, outer = TRUE); return(invisible(NULL)) }

    plot(NA, NA, xlim = c(1, xmax), ylim = c(0, ymax),
         xlab = "k (top-k cells by direction)", ylab = "precision@k",
         main = "Precision@k (overlay)")

    # baselines
    if (!is.null(curN) && is.finite(curN$p_base)) abline(h = curN$p_base, col = ac(col_naive, 0.5), lty = 3)
    if (!is.null(curT) && is.finite(curT$p_base)) abline(h = curT$p_base, col = ac(col_term, 0.5),  lty = 3)

    # curves
    if (!is.null(curN) && length(curN$k)) {
      lines(curN$k, curN$prec, col = col_naive, lwd = lwd)
      if (isTRUE(show_rug) && length(curN$ranks_pos)) rug(curN$ranks_pos, col = ac(col_naive, point_alpha))
      if (is.finite(curN$k_star) && curN$k_star > 0L) {
        points(curN$k_star, curN$prec_star, pch = 19, cex = point_cex, col = col_naive)
      }
    }
    if (!is.null(curT) && length(curT$k)) {
      lines(curT$k, curT$prec, col = col_term, lwd = lwd)
      if (isTRUE(show_rug) && length(curT$ranks_pos)) rug(curT$ranks_pos, col = ac(col_term, point_alpha), side = 1)
      if (is.finite(curT$k_star) && curT$k_star > 0L) {
        points(curT$k_star, curT$prec_star, pch = 19, cex = point_cex, col = col_term)
      }
    }

    legend("bottomright",
           legend = c(
             sprintf("Naïve:  E_naive = %.3f  (p=%.2f)", E_obj$E_naive, ifelse(is.null(curN), NA, curN$p_base)),
             sprintf("Terminal: E_term  = %.3f  (p=%.2f)", E_obj$E_term, ifelse(is.null(curT), NA, curT$p_base)),
             if (!is.na(E_obj$E_comp)) sprintf("E_comp = %.3f", E_obj$E_comp) else NULL
           ),
           col = c(col_naive, col_term, "black")[seq_len(2 + !is.na(E_obj$E_comp))],
           lwd = c(lwd, lwd, NA)[seq_len(2 + !is.na(E_obj$E_comp))],
           pch = c(19, 19, NA)[seq_len(2 + !is.na(E_obj$E_comp))],
           bty = "n", cex = 0.9)

    mtext(main, outer = TRUE, cex = 1.2, line = 0)

  } else {
    par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))

    # panel helper
    panel_plot <- function(cur, title, col_line) {
      if (is.null(cur) || !length(cur$k)) {
        plot.new(); title(main = paste0(title, " (no labels)")); return(invisible(NULL))
      }
      plot(cur$k, cur$prec, type = "l", col = col_line, lwd = lwd,
           xlab = "k (top-k cells)", ylab = "precision@k", main = title,
           ylim = c(0, max(cur$prec, cur$p_base, na.rm = TRUE)))
      abline(h = cur$p_base, col = ac(col_line, 0.5), lty = 3)
      if (isTRUE(show_rug) && length(cur$ranks_pos)) rug(cur$ranks_pos, col = ac(col_line, point_alpha))
      if (is.finite(cur$k_star) && cur$k_star > 0L) {
        points(cur$k_star, cur$prec_star, pch = 19, cex = point_cex, col = col_line)
      }
    }

    panel_plot(curN, sprintf("Naïve (E=%.3f)", E_obj$E_naive), col_naive)
    panel_plot(curT, sprintf("Terminal (E=%.3f)", E_obj$E_term), col_term)

    mtext(main, outer = TRUE, cex = 1.2, line = 0)
  }

  invisible(NULL)
}

