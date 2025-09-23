# ===============================
# Multi-trajectory wrapper
# ===============================
#' Compute DOPE metrics for multiple linear trajectories
#'
#' @description
#' Runs the full DOPE pipeline (Directionality **D**, Order consistency **O**,
#' Program coherence **P**, Endpoints validity **E**, and the combined
#' `DOPE_score`) across **multiple linear** pseudotime vectors sharing the
#' same expression object. Results are collected per trajectory and summarized
#' in a comparison table with the best trajectory highlighted.
#'
#' @details
#' This function is a convenience wrapper that repeatedly calls
#' [compute_single_DOPE_linear()] for each linear trajectory provided in
#' `pseudotime_list`. It supports optional parallelization via the
#' **parallel** package (fork/PSOCK, depending on platform).
#'
#' If `pseudotime_list` is unnamed, trajectories are auto-named as
#' `"trajectory_1"`, `"trajectory_2"`, etc.
#'
#' @param expr_or_seurat A gene expression matrix (genes × cells),
#'   a `dgCMatrix`, a `data.frame` (coerced to matrix), or a Seurat object.
#'   If Seurat, `assay` and `slot` control which assay/slot is used.
#' @param pseudotime_list A named or unnamed **list** of numeric pseudotime
#'   vectors (one per trajectory). Each vector must have length equal to the
#'   number of columns/cells of `expr_or_seurat` (after extraction).
#' @param naive_markers Character vector of gene symbols/IDs for naïve/early
#'   programs (used by D/O/E/P as relevant).
#' @param terminal_markers Character vector of gene symbols/IDs for terminal/late
#'   programs.
#' @param cluster_labels Optional character/factor vector of cluster labels
#'   (one per cell) used by E (when `E_method = "clusters"` or `"combined"`).
#' @param naive_clusters Optional character vector of cluster names considered
#'   naïve (E Case 2 / cluster-based).
#' @param terminal_clusters Optional character vector of cluster names considered
#'   terminal (E Case 2 / cluster-based).
#' @param pathways Optional list of pathway name → gene vector (used by P).
#' @param gene_sets Optional list of additional gene sets name → gene vector
#'   (also used by P).
#' @param E_method Endpoints validity method. One of `"gmm"`, `"clusters"`,
#'   or `"combined"`. See [compute_single_DOPE_linear()] for details.
#' @param id_type Gene identifier type for markers/pathways. One of
#'   `"SYMBOL"`, `"ENSEMBL"`, or `"ENTREZID"`.
#' @param species Species name used for ID mapping (e.g., `"Homo sapiens"`).
#' @param assay If `expr_or_seurat` is a Seurat object, the assay to use.
#'   If `NULL`, defaults to `Seurat::DefaultAssay()`.
#' @param slot If `expr_or_seurat` is a Seurat object, the assay slot to
#'   extract (e.g., `"data"`, `"counts"`). Default `"data"`.
#' @param min_remaining Minimum number of genes to retain after filtering
#'   within sub-steps (safety check).
#' @param min_fraction Minimum retained fraction of genes after filtering
#'   (safety check).
#' @param min_genes_per_module Minimum number of genes required per program/
#'   module to compute stable scores.
#' @param plot_E Logical; if `TRUE`, produce diagnostic density plots when
#'   `E_method` uses GMM.
#' @param verbose Logical; print progress.
#' @param parallel Logical; if `TRUE` and multiple trajectories are provided,
#'   use **parallel** workers.
#' @param n_cores Integer number of cores. If `NULL`, uses `detectCores()-1`.
#' @param drop_unused_levels Logical; drop unused factor levels (where relevant).
#' @param tol Numerical tolerance for internal numerical checks.
#'
#' @return
#' An object of class `"multi_dope_results"` with components:
#' \itemize{
#'   \item \code{results}: named list of per-trajectory DOPE results
#'         (each as returned by \code{compute_single_DOPE_linear()}).
#'   \item \code{comparison_summary}: \code{data.frame} summarizing
#'         D/O/P/E and \code{DOPE_score} per trajectory.
#'   \item \code{best_trajectory}: character scalar with the name of the
#'         highest-scoring trajectory (or \code{NA} if none valid).
#'   \item \code{n_trajectories}: integer count.
#'   \item \code{method_info}: list of key settings used.
#' }
#'
#' @seealso [compute_single_DOPE_linear()], plotting helpers like
#'   \code{plot.multi_dope_results()}.
#'
#' @examples
#' \dontrun{
#' # expr: genes x cells matrix; pt1, pt2: numeric pseudotime vectors (length = ncol(expr))
#' res <- compute_multi_DOPE_linear(
#'   expr_or_seurat = expr,
#'   pseudotime_list = list(linear_a = pt1, linear_b = pt2),
#'   naive_markers = c("TCF7","LEF1"),
#'   terminal_markers = c("GZMB","PRF1"),
#'   E_method = "combined",
#'   parallel = TRUE
#' )
#' res$comparison_summary
#' res$best_trajectory
#' }
#'
#' @export
compute_multi_DOPE_linear <- function(expr_or_seurat,
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
                                      plot_E = TRUE,
                                      verbose = TRUE,
                                      parallel = FALSE,
                                      n_cores = NULL,
                                      drop_unused_levels = TRUE,
                                      tol = 1e-8) {

  E_method <- match.arg(E_method)
  id_type  <- match.arg(id_type)

  if (!is.list(pseudotime_list)) stop("pseudotime_list must be a list of pseudotime vectors")
  if (is.null(names(pseudotime_list))) names(pseudotime_list) <- paste0("trajectory_", seq_along(pseudotime_list))
  n_trajectories <- length(pseudotime_list)

  if (verbose) {
    cat("Computing DOPE metrics for", n_trajectories, "linear trajectories...\n")
    cat("==========================================================\n")
  }

  # Set up parallel
  use_parallel <- FALSE
  if (parallel && n_trajectories > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not available, running sequentially")
    } else {
      if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
      if (verbose) cat("Using", n_cores, "cores for parallel processing\n")
      use_parallel <- TRUE
    }
  }

  compute_single_DOPE_linear <- function(trajectory_name, pseudotime_vector) {
    if (verbose) cat("\n--- Processing trajectory:", trajectory_name, "---\n")
    tryCatch({
      res <- compute_single_DOPE_linear(
        expr_or_seurat   = expr_or_seurat,
        pseudotime       = pseudotime_vector,
        naive_markers    = naive_markers,
        terminal_markers = terminal_markers,
        cluster_labels   = cluster_labels,
        naive_clusters   = naive_clusters,
        terminal_clusters = terminal_clusters,
        pathways         = pathways,
        gene_sets        = gene_sets,
        E_method         = E_method,
        id_type          = id_type,
        species          = species,
        assay            = assay,
        slot             = slot,
        min_remaining    = min_remaining,
        min_fraction     = min_fraction,
        min_genes_per_module = min_genes_per_module,
        plot_E           = plot_E,
        verbose          = verbose,
        tol              = tol
      )
      res$trajectory_name <- trajectory_name
      res
    }, error = function(e) {
      warning(paste("Error processing trajectory", trajectory_name, ":", e$message))
      list(
        trajectory_name = trajectory_name,
        D = list(D_naive = NA, D_term = NA),
        O = list(O = NA),
        P = list(P = NA),
        E = list(E_naive = NA, E_term = NA, E_comp = NA),
        DOPE_score = NA,
        errors = list(D_error = e$message, O_error = e$message, P_error = e$message, E_error = e$message)
      )
    })
  }

  if (use_parallel) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("compute_single_DOPE_linear",
                  "metrics_D","metrics_O","metrics_P","metrics_E",
                  ".get_expr",".per_cell_score",".align_to_cells","subset_by_clusters","%||%", ".has_formal"),
      envir = environment()
    )
    # If you need packages on workers, load them here:
    # parallel::clusterEvalQ(cl, { library(Seurat) })

    idxs <- seq_along(pseudotime_list)
    trajectory_results <- parallel::parLapply(
      cl, idxs,
      function(i) {
        nm <- names(pseudotime_list)[i]
        compute_single_dope(nm, pseudotime_list[[i]])
      }
    )
    names(trajectory_results) <- names(pseudotime_list)
  } else {
    trajectory_results <- Map(compute_single_dope, names(pseudotime_list), pseudotime_list)
    names(trajectory_results) <- names(pseudotime_list)
  }

  comparison_summary <- create_comparison_summary(trajectory_results)

  valid_scores <- comparison_summary$DOPE_score[!is.na(comparison_summary$DOPE_score)]
  best_trajectory <- if (length(valid_scores) > 0) {
    comparison_summary$trajectory[which.max(comparison_summary$DOPE_score)]
  } else NA_character_

  multi_results <- list(
    results = trajectory_results,
    comparison_summary = comparison_summary,
    best_trajectory = best_trajectory,
    n_trajectories = n_trajectories,
    method_info = list(
      E_method = E_method,
      id_type = id_type,
      species = species,
      parallel = use_parallel,
      n_cores = if (use_parallel) n_cores else NA_integer_
    )
  )
  class(multi_results) <- "multi_dope_results"

  if (verbose) {
    cat("\n==========================================================\n")
    cat("Multi-trajectory DOPE analysis (linear) complete!\n")
    if (!is.na(best_trajectory)) {
      best_score <- comparison_summary$DOPE_score[comparison_summary$trajectory == best_trajectory]
      cat(sprintf("Best trajectory: %s (DOPE score: %.3f)\n", best_trajectory, best_score))
    }
  }
}
