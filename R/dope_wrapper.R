#' DOPE: Comprehensive pseudotime trajectory quality assessment for multiple trajectories
#'
#' @description
#' Compute all four DOPE metrics (Directionality, Order, Program coherence, Endpoints validity)
#' for multiple pseudotime trajectories in a single function call. This function loops through
#' each provided pseudotime and computes DOPE metrics, allowing comparison across different
#' trajectory inference methods or parameters.
#'
#' @param expr_or_seurat Matrix (genes x cells) **or** a Seurat object containing expression data.
#' @param pseudotime_list Named list of numeric vectors, where each vector contains pseudotime values
#'                       for one trajectory method. Names will be used as trajectory identifiers.
#'                       Each vector should have length = #cells and can be named by cell for alignment.
#' @param naive_markers Character vector of naive/root marker gene names.
#' @param terminal_markers Character vector of terminal marker gene names.
#' @param cluster_labels Character vector of cluster labels (optional, for E metric cluster method).
#' @param naive_clusters Character vector of cluster names considered as naive (for E metric).
#' @param terminal_clusters Character vector of cluster names considered as terminal (for E metric).
#' @param pathways Character vector of KEGG pathway queries for P metric (e.g., "hsa04660" or name fragments).
#'                 Ignored if \code{gene_sets} is provided.
#' @param gene_sets Optional named list of gene vectors for P metric (ID space given by \code{id_type}).
#'                  If supplied, skips KEGG lookup.
#' @param E_method Character, one of "gmm" (Case 1), "clusters" (Case 2), or "combined" (Case 3) for E metric.
#' @param id_type One of "SYMBOL","ENSEMBL","ENTREZID" for P metric gene mapping. Default "SYMBOL".
#' @param species Species label for msigdbr in P metric (default "Homo sapiens").
#' @param assay,slot If a Seurat object is provided, which assay/slot to pull for P metric (default: current assay, slot "data").
#' @param min_remaining,min_fraction,min_genes_per_module Parameters for P metric pathway filtering.
#' @param tol Numeric tolerance for O metric tied changes (default 1e-8).
#' @param plot_E Logical, whether to generate density plots for E metric GMM.
#' @param verbose Logical; print progress messages.
#' @param parallel Logical; whether to use parallel processing for multiple trajectories (requires parallel package).
#' @param n_cores Integer; number of cores to use if parallel=TRUE (default: detectCores()-1).
#'
#' @return A list of class "multi_dope_results" containing:
#' \item{results}{Named list where each element contains DOPE results for one trajectory}
#' \item{comparison_summary}{Data frame comparing all trajectories across metrics}
#' \item{best_trajectory}{Name of trajectory with highest composite DOPE score}
#' \item{n_trajectories}{Number of trajectories analyzed}
#' \item{method_info}{Analysis parameters used}
#'
#' @examples
#' # Create multiple pseudotime vectors
#' pseudotime_methods <- list(
#'   "monocle3" = monocle3_pseudotime,
#'   "slingshot" = slingshot_pseudotime,
#'   "palantir" = palantir_pseudotime
#' )
#'
#' # Compare multiple trajectory methods
#' multi_results <- compute_multi_DOPE(
#'   expr = expression_matrix,
#'   pseudotime_list = pseudotime_methods,
#'   naive_markers = c("TCF7", "LEF1", "CCR7"),
#'   terminal_markers = c("GZMB", "PRF1", "IFNG"),
#'   pathways = c("hsa04660", "hsa04658"),
#'   E_method = "gmm",
#'   parallel = TRUE
#' )
#'
#' @export
compute_multi_DOPE <- function(expr_or_seurat,
                               pseudotime_list,
                               naive_markers,
                               terminal_markers,
                               cluster_labels = NULL,
                               naive_clusters = NULL,
                               terminal_clusters = NULL,
                               pathways = NULL,
                               gene_sets = NULL,
                               E_method = c("gmm", "clusters", "combined"),
                               id_type = c("SYMBOL", "ENSEMBL", "ENTREZID"),
                               species = "Homo sapiens",
                               assay = NULL,
                               slot = "data",
                               min_remaining = 10,
                               min_fraction = 0.20,
                               min_genes_per_module = 3,
                               tol = 1e-8,
                               plot_E = TRUE,
                               verbose = TRUE,
                               parallel = FALSE,
                               n_cores = NULL) {

  E_method <- match.arg(E_method)
  id_type <- match.arg(id_type)

  # Input validation
  if (!is.list(pseudotime_list)) {
    stop("pseudotime_list must be a list of pseudotime vectors")
  }

  if (is.null(names(pseudotime_list))) {
    names(pseudotime_list) <- paste0("trajectory_", seq_along(pseudotime_list))
  }

  n_trajectories <- length(pseudotime_list)

  if (verbose) {
    cat("Computing DOPE metrics for", n_trajectories, "trajectories...\n")
    cat("==========================================================\n")
  }

  # Set up parallel processing if requested
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not available, running sequentially")
      parallel <- FALSE
    } else {
      if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
      }
      if (verbose) {
        cat("Using", n_cores, "cores for parallel processing\n")
      }
    }
  }

  # Function to compute DOPE for a single trajectory
  compute_single_dope <- function(trajectory_name, pseudotime_vector) {
    if (verbose) {
      cat("\n--- Processing trajectory:", trajectory_name, "---\n")
    }

    tryCatch({
      result <- compute_DOPE_single(
        expr_or_seurat = expr_or_seurat,
        pseudotime = pseudotime_vector,
        naive_markers = naive_markers,
        terminal_markers = terminal_markers,
        cluster_labels = cluster_labels,
        naive_clusters = naive_clusters,
        terminal_clusters = terminal_clusters,
        pathways = pathways,
        gene_sets = gene_sets,
        E_method = E_method,
        id_type = id_type,
        species = species,
        assay = assay,
        slot = slot,
        min_remaining = min_remaining,
        min_fraction = min_fraction,
        min_genes_per_module = min_genes_per_module,
        tol = tol,
        plot_E = plot_E,
        verbose = verbose
      )

      # Add trajectory name to result
      result$trajectory_name <- trajectory_name
      return(result)

    }, error = function(e) {
      warning(paste("Error processing trajectory", trajectory_name, ":", e$message))
      return(list(
        trajectory_name = trajectory_name,
        D = list(D_naive = NA, D_term = NA),
        O = list(O = NA),
        P = list(P = NA),
        E = list(E_naive = NA, E_term = NA, E_comp = NA),
        DOPE_score = NA,
        error = e$message
      ))
    })
  }

  # Process trajectories
  if (parallel && n_trajectories > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export necessary objects to cluster
    parallel::clusterExport(cl, c("compute_DOPE_single", "metrics_D", "metrics_O",
                                  "metrics_P", "metrics_E"), envir = environment())

    trajectory_results <- parallel::Map(
      compute_single_dope,
      names(pseudotime_list),
      pseudotime_list,
      mc.cores = n_cores
    )
  } else {
    # Sequential processing
    trajectory_results <- Map(
      compute_single_dope,
      names(pseudotime_list),
      pseudotime_list
    )
  }

  names(trajectory_results) <- names(pseudotime_list)

  # Create comparison summary
  comparison_summary <- create_comparison_summary(trajectory_results)

  # Find best trajectory
  valid_scores <- comparison_summary$DOPE_score[!is.na(comparison_summary$DOPE_score)]
  if (length(valid_scores) > 0) {
    best_idx <- which.max(comparison_summary$DOPE_score)
    best_trajectory <- comparison_summary$trajectory[best_idx]
  } else {
    best_trajectory <- NA
  }

  # Create final results object
  multi_results <- list(
    results = trajectory_results,
    comparison_summary = comparison_summary,
    best_trajectory = best_trajectory,
    n_trajectories = n_trajectories,
    method_info = list(
      E_method = E_method,
      id_type = id_type,
      species = species,
      parallel = parallel,
      n_cores = if (parallel) n_cores else NA
    )
  )

  class(multi_results) <- "multi_dope_results"

  if (verbose) {
    cat("\n==========================================================\n")
    cat("Multi-trajectory DOPE analysis complete!\n")
    if (!is.na(best_trajectory)) {
      best_score <- comparison_summary$DOPE_score[comparison_summary$trajectory == best_trajectory]
      cat(sprintf("Best trajectory: %s (DOPE score: %.3f)\n", best_trajectory, best_score))
    }
  }

  return(multi_results)
}

#' Helper function to compute DOPE for a single trajectory
#' (This would be the original compute_DOPE function with minimal modifications)
.get_expr <- function(x, assay = NULL, slot = "data") {
  if (inherits(x, "Seurat")) {
    if (is.null(assay)) assay <- Seurat::DefaultAssay(x)
    return(as.matrix(Seurat::GetAssayData(x, assay = assay, slot = slot)))
  } else if (is.matrix(x)) {
    return(x)
  } else {
    stop("expr_or_seurat must be a matrix or a Seurat object")
  }
}

# function for getting naive or terminal gene scores
.per_cell_score <- function(expr, markers) {
  genes <- intersect(markers, rownames(expr))
  if (length(genes) == 0) {
    warning("No markers found in expression matrix.")
    return(rep(NA_real_, ncol(expr)))
  }
  mark_expr <- expr[genes, , drop = FALSE]
  colSums(mark_expr) / length(genes)
}



# Helper function for null coalescing (moved to top for availability)
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

#' Helper function to extract expression data from Seurat object or matrix
.get_expr <- function(x, assay = NULL, slot = "data") {
  if (inherits(x, "Seurat")) {
    if (is.null(assay)) assay <- Seurat::DefaultAssay(x)
    return(as.matrix(Seurat::GetAssayData(x, assay = assay, slot = slot)))
  } else if (is.matrix(x)) {
    return(x)
  } else {
    stop("expr_or_seurat must be a matrix or a Seurat object")
  }
}

#' Function for getting naive or terminal gene scores
.per_cell_score <- function(expr, markers) {
  # Ensure expr is a matrix
  if (!is.matrix(expr)) {
    expr <- .get_expr(expr)
  }

  genes <- intersect(markers, rownames(expr))
  if (length(genes) == 0) {
    warning("No markers found in expression matrix.")
    return(rep(NA_real_, ncol(expr)))
  }
  mark_expr <- expr[genes, , drop = FALSE]
  colSums(mark_expr) / length(genes)
}

#' Compute DOPE for a single trajectory
#'
#' @param expr_or_seurat Expression matrix or Seurat object
#' @param pseudotime Vector of pseudotime values
#' @param naive_markers Vector of naive marker genes
#' @param terminal_markers Vector of terminal marker genes
#' @param cluster_labels Vector of cluster labels (optional)
#' @param naive_clusters Vector of naive cluster identities (optional)
#' @param terminal_clusters Vector of terminal cluster identities (optional)
#' @param pathways List of pathway gene sets (optional)
#' @param gene_sets Alternative gene sets (optional)
#' @param E_method Method for computing E metric
#' @param id_type Gene identifier type
#' @param species Species for pathway analysis
#' @param assay Assay name for Seurat objects
#' @param slot Slot name for Seurat objects
#' @param min_remaining Minimum remaining genes for pathway analysis
#' @param min_fraction Minimum fraction for pathway analysis
#' @param min_genes_per_module Minimum genes per pathway module
#' @param tol Tolerance for numerical computations
#' @param plot_E Whether to plot E metric results
#' @param verbose Whether to print verbose output
#' @return List containing DOPE results
compute_DOPE_single <- function(expr_or_seurat,
                                pseudotime,
                                naive_markers,
                                terminal_markers,
                                cluster_labels = NULL,
                                naive_clusters = NULL,
                                terminal_clusters = NULL,
                                pathways = NULL,
                                gene_sets = NULL,
                                E_method = c("gmm", "clusters", "combined"),
                                id_type = c("SYMBOL", "ENSEMBL", "ENTREZID"),
                                species = "Homo sapiens",
                                assay = NULL,
                                slot = "data",
                                min_remaining = 10,
                                min_fraction = 0.20,
                                min_genes_per_module = 3,
                                tol = 1e-8,
                                plot_E = TRUE,
                                verbose = FALSE) {

  E_method <- match.arg(E_method)
  id_type  <- match.arg(id_type)

  # Extract expression matrix
  expr <- .get_expr(expr_or_seurat, assay = assay, slot = slot)

  # Validate inputs
  if (length(pseudotime) != ncol(expr)) {
    stop("Length of pseudotime must match number of cells in expression data")
  }

  if (!is.null(cluster_labels) && length(cluster_labels) != ncol(expr)) {
    stop("Length of cluster_labels must match number of cells in expression data")
  }

  # Compute D metrics
  D_res <- tryCatch({
    metrics_D(expr, naive_markers, terminal_markers, pseudotime)
  }, error = function(e) {
    if (verbose) message("D metric failed: ", e$message)
    list(D_naive = NA_real_, D_term = NA_real_, error = e$message)
  })

  # Compute per-cell scores for markers
  naive_scores <- .per_cell_score(expr_or_seurat, naive_markers)
  term_scores  <- .per_cell_score(expr_or_seurat, terminal_markers)

  # Compute O metrics
  O_res <- tryCatch({
    metrics_O(expr,
              naive_markers,
              terminal_markers,
              pseudotime,
              tol = tol)
  }, error = function(e) {
    if (verbose) message("O metric failed: ", e$message)
    list(O = NA_real_, error = e$message)
  })

  # Compute P metrics
  P_res <- tryCatch({
    metrics_P(expr,
              pseudotime,
              pathways = pathways,
              gene_sets = gene_sets,
              naive_markers = naive_markers,
              terminal_markers = terminal_markers,
              id_type = id_type,
              species = species,
              assay = assay,
              slot = slot,
              min_remaining = min_remaining,
              min_fraction = min_fraction,
              min_genes_per_module = min_genes_per_module,
              verbose = verbose)
  }, error = function(e) {
    if (verbose) message("P metric failed: ", e$message)
    list(P = NA_real_, error = e$message)
  })

  # Compute E metrics
  E_res <- tryCatch({
    metrics_E(pseudotime,
              cluster_labels = cluster_labels,
              naive_marker_scores = naive_scores,
              terminal_marker_scores = term_scores,
              naive_clusters = naive_clusters,
              terminal_clusters = terminal_clusters,
              method = E_method,
              plot = plot_E)
  }, error = function(e) {
    if (verbose) message("E metric failed: ", e$message)
    list(E_naive = NA_real_, E_term = NA_real_, E_comp = NA_real_, error = e$message)
  })

  # Calculate composite DOPE score
  comp_vec <- c(
    D_res$D_naive %||% NA_real_,
    D_res$D_term %||% NA_real_,
    O_res$O %||% NA_real_,
    P_res$P %||% NA_real_,
    E_res$E_comp %||% NA_real_
  )

  DOPE_score <- if (all(is.na(comp_vec))) {
    NA_real_
  } else {
    mean(comp_vec, na.rm = TRUE)
  }

  # Compile results
  out <- list(
    D = list(
      D_naive = D_res$D_naive %||% NA_real_,
      D_term = D_res$D_term %||% NA_real_
    ),
    O = list(O = O_res$O %||% NA_real_),
    P = list(P = P_res$P %||% NA_real_),
    E = list(
      E_naive = E_res$E_naive %||% NA_real_,
      E_term = E_res$E_term %||% NA_real_,
      E_comp = E_res$E_comp %||% NA_real_
    ),
    DOPE_score = DOPE_score,
    # Include any errors that occurred
    errors = list(
      D_error = D_res$error,
      O_error = O_res$error,
      P_error = P_res$error,
      E_error = E_res$error
    )
  )

  class(out) <- "dope_results"
  return(out)
}

#' Create comparison summary across trajectories
#'
#' @param trajectory_results List of trajectory results from compute_DOPE_single
#' @return Data frame with comparison summary and rankings
create_comparison_summary <- function(trajectory_results) {

  if (length(trajectory_results) == 0) {
    stop("trajectory_results cannot be empty")
  }

  # Initialize summary data frame
  summary_data <- data.frame(
    trajectory = character(),
    D_naive = numeric(),
    D_term = numeric(),
    O = numeric(),
    P = numeric(),
    E_naive = numeric(),
    E_term = numeric(),
    E_comp = numeric(),
    DOPE_score = numeric(),
    has_error = logical(),
    stringsAsFactors = FALSE
  )

  # Extract data for each trajectory
  for (i in seq_along(trajectory_results)) {
    result <- trajectory_results[[i]]
    trajectory_name <- names(trajectory_results)[i] %||% paste0("Trajectory_", i)

    # Check if any errors occurred
    has_error <- any(sapply(result$errors, function(x) !is.null(x)))

    # Create row for this trajectory
    row_data <- data.frame(
      trajectory = trajectory_name,
      D_naive = result$D$D_naive %||% NA_real_,
      D_term = result$D$D_term %||% NA_real_,
      O = result$O$O %||% NA_real_,
      P = result$P$P %||% NA_real_,
      E_naive = result$E$E_naive %||% NA_real_,
      E_term = result$E$E_term %||% NA_real_,
      E_comp = result$E$E_comp %||% NA_real_,
      DOPE_score = result$DOPE_score %||% NA_real_,
      has_error = has_error,
      stringsAsFactors = FALSE
    )

    summary_data <- rbind(summary_data, row_data)
  }

  # Add rankings for each metric (higher scores = better ranks)
  metrics <- c("D_naive", "D_term", "O", "P", "E_naive", "E_term", "E_comp", "DOPE_score")
  for (metric in metrics) {
    rank_col <- paste0(metric, "_rank")
    summary_data[[rank_col]] <- rank(-summary_data[[metric]],
                                     na.last = "keep",
                                     ties.method = "min")
  }

  # Sort by DOPE score (best first)
  summary_data <- summary_data[order(summary_data$DOPE_score,
                                     decreasing = TRUE,
                                     na.last = TRUE), ]

  # Reset row names
  rownames(summary_data) <- NULL

  return(summary_data)
}
#' Print method for multi-DOPE results
#'
#' @param x A multi_dope_results object
#' @param ... Additional arguments (not used)
#'
#' @export
print.multi_dope_results <- function(x, ...) {
  cat("Multi-Trajectory DOPE Analysis\n")
  cat("==============================\n\n")

  cat("Analysis info:\n")
  cat(sprintf("  Trajectories analyzed: %d\n", x$n_trajectories))
  cat(sprintf("  E method: %s\n", x$method_info$E_method))
  if (x$method_info$parallel) {
    cat(sprintf("  Parallel processing: %d cores\n", x$method_info$n_cores))
  }
  cat("\n")

  cat("Trajectory Comparison:\n")
  cat("---------------------\n")

  # Print top trajectories
  top_n <- min(5, nrow(x$comparison_summary))
  top_trajectories <- x$comparison_summary[1:top_n, ]

  for (i in seq_len(nrow(top_trajectories))) {
    traj <- top_trajectories[i, ]
    rank_suffix <- switch(as.character(i), "1" = "st", "2" = "nd", "3" = "rd", "th")
    dope_score <- if (is.na(traj$DOPE_score)) "NA" else sprintf("%.3f", traj$DOPE_score)

    cat(sprintf("%d%s: %s (DOPE: %s)\n", i, rank_suffix, traj$trajectory, dope_score))
  }

  if (!is.na(x$best_trajectory)) {
    cat(sprintf("\nBest trajectory: %s\n", x$best_trajectory))
  }

  cat("\nUse summary() for detailed comparison table\n")
}

#' Summary method for multi-DOPE results
#'
#' @param object A multi_dope_results object
#' @param ... Additional arguments (not used)
#'
#' @export
summary.multi_dope_results <- function(object, ...) {
  print(object)

  cat("\nDetailed Comparison Table:\n")
  cat("=========================\n")

  # Create a nicely formatted table
  display_cols <- c("trajectory", "DOPE_score", "DOPE_score_rank",
                    "D_term", "O", "P", "E_comp")

  display_data <- object$comparison_summary[, display_cols]

  # Format numeric columns
  numeric_cols <- sapply(display_data, is.numeric)
  display_data[numeric_cols] <- lapply(display_data[numeric_cols], function(x) {
    ifelse(is.na(x), "NA", sprintf("%.3f", x))
  })

  # Print the table
  print(display_data, row.names = FALSE)

  # Show detailed results for best trajectory if available
  if (!is.na(object$best_trajectory)) {
    cat(sprintf("\nDetailed results for best trajectory (%s):\n", object$best_trajectory))
    cat("==========================================\n")
    best_result <- object$results[[object$best_trajectory]]
    if (inherits(best_result, "dope_results")) {
      print(best_result)
    }
  }

  invisible(object)
}

#' Plot comparison of trajectories across DOPE metrics
#'
#' @param multi_dope_results A multi_dope_results object
#' @param metrics Character vector of metrics to plot (default: main DOPE metrics)
#' @param type Plot type: "bar"or "radar"
#'
#' @export
plot.multi_dope_results <- function(multi_dope_results,
                                    metrics = c("D_naive", "D_term", "O", "P", "E_naive", "E_term", "DOPE_score"),
                                    type = "bar") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }

  # Extract data for plotting
  plot_data <- multi_dope_results$comparison_summary[, c("trajectory", metrics)]

  if (type == "bar") {
    # Reshape data for ggplot
    if (!requireNamespace("reshape2", quietly = TRUE)) {
      stop("reshape2 package required for bar plots")
    }

    plot_data_long <- reshape2::melt(plot_data, id.vars = "trajectory",
                                     variable.name = "metric", value.name = "score")

    p <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = trajectory, y = score, fill = metric)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "DOPE Metrics Comparison Across Trajectories",
                    x = "Trajectory Method", y = "Score") +
      ggplot2::ylim(0, 1)

      return(p)
      }
      else if (type == "radar") {
      # Check for fmsb package
      if (!requireNamespace("fmsb", quietly = TRUE)) {
        stop("fmsb package required for radar plots. Install with: install.packages('fmsb')")
      }

      # Prepare data for radar plot
      # fmsb::radarchart expects data with max/min values in first two rows
      radar_data <- plot_data[, metrics, drop = FALSE]

      # Handle missing values - replace with 0 for radar plot
      radar_data[is.na(radar_data)] <- 0

      # Ensure all metric values are between 0 and 1 for radar plot
      for (metric in metrics) {
        if (metric %in% names(radar_data)) {
          radar_data[[metric]] <- pmax(0, pmin(1, radar_data[[metric]]))
        }
      }

      # Add max and min rows (required by fmsb)
      radar_data <- rbind(rep(1, length(metrics)), rep(0, length(metrics)), radar_data)
      rownames(radar_data) <- c("Max", "Min", plot_data$trajectory)

      # Set up plotting parameters
      n_trajectories <- nrow(plot_data)

      # Add transparency for fill and line
      colors <- rainbow(n_trajectories, alpha = 0.3)
      line_colors <- rainbow(n_trajectories, alpha = 0.8)

      # Create the radar chart
      # Note: fmsb::radarchart returns NULL and plots directly
      # We'll capture this as a base R plot and convert to a function that can be called

      plot_radar <- function() {
        # Set up plot layout for legend
        layout(matrix(c(1, 2), ncol = 2), widths = c(3, 1))

        # Create radar chart
        fmsb::radarchart(
          radar_data,
          axistype = 1,
          pcol = line_colors,
          pfcol = colors,
          plwd = 2,
          plty = 1,
          cglcol = "grey",
          cglty = 1,
          axislabcol = "grey",
          caxislabels = rep("", 5),
          cglwd = 0.5,
          vlcex = 0,
          title = "DOPE Metrics Radar Chart"
        )

        # Add legend
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend("center",
               legend = plot_data$trajectory,
               col = line_colors,
               lty = 1,
               lwd = 2,
               cex = 0.8,
               bty = "n")

        # Reset layout
        layout(1)
      }

      # For immediate plotting
      plot_radar()

      # Return a function that can recreate the plot
      return(invisible(plot_radar))

      }
      else if (type == "heatmap") {
      # Prepare data for heatmap
      if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop("reshape2 package required for heatmap plots")
      }

      # Reshape data for heatmap
      heatmap_data <- reshape2::melt(plot_data, id.vars = "trajectory",
                                     variable.name = "metric", value.name = "score")

      # Create heatmap
      p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = metric, y = trajectory, fill = score)) +
        ggplot2::geom_tile(color = "white", size = 0.5) +
        ggplot2::scale_fill_gradient2(
          low = "red",
          mid = "yellow",
          high = "green",
          midpoint = 0.5,
          name = "Score",
          na.value = "grey90"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::labs(
          title = "DOPE Metrics Heatmap",
          x = "Metric",
          y = "Trajectory Method"
        ) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", score)),
                           color = "black", size = 3)

      return(p)

    }

    else if (type == "heatmap") {
    # Prepare data for heatmap
    if (!requireNamespace("reshape2", quietly = TRUE)) {
      stop("reshape2 package required for heatmap plots")
    }

    # Reshape data for heatmap
    heatmap_data <- reshape2::melt(plot_data, id.vars = "trajectory",
                                   variable.name = "metric", value.name = "score")

    # Create heatmap
    p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = metric, y = trajectory, fill = score)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::scale_fill_gradient2(
        low = "red",
        mid = "yellow",
        high = "green",
        midpoint = 0.5,
        name = "Score",
        na.value = "grey90"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5)
      ) +
      ggplot2::labs(
        title = "DOPE Metrics Heatmap",
        x = "Metric",
        y = "Trajectory Method"
      ) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", score)),
                         color = "black", size = 3)

    return(p)

  } else {
    stop("Plot type must be one of: 'bar', 'radar', or 'heatmap'")
  }
}
