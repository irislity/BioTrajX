#' DOPE: Comprehensive pseudotime trajectory quality assessment for multiple trajectories
#'
#' @description
#' Compute all four DOPE metrics (Directionality, Order, Program coherence, Endpoints validity)
#' for multiple pseudotime trajectories in a single function call. This function loops through
#' each provided pseudotime and computes DOPE metrics, allowing comparison across different
#' trajectory inference methods or parameters. Supports linear or branched trajectories; for
#' branched analyses, cells can be subset using cluster labels before metric computation.
#'
#' @param expr_or_seurat Matrix (genes x cells) **or** a Seurat object containing expression data.
#' @param pseudotime_list Named list of numeric vectors, where each vector contains pseudotime values
#'                       for one trajectory method. Names will be used as trajectory identifiers.
#'                       Each vector should have length = #cells and can be named by cell for alignment.
#' @param naive_markers Character vector of naive/root marker gene names.
#' @param terminal_markers Character vector of terminal marker gene names.
#' @param cluster_labels Character scalar (Seurat meta column name) **or** per-cell vector of cluster labels.
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
#' @param plot_E Logical, whether to generate density plots for E metric GMM.
#' @param verbose Logical; print progress messages.
#' @param parallel Logical; whether to use parallel processing for multiple trajectories.
#' @param n_cores Integer; number of cores to use if parallel=TRUE (default: detectCores()-1).
#' @param trajectory One of c("linear","branched"). If "branched", will subset cells by clusters before metrics.
#' @param branch_include Character vector of cluster names to keep (optional; branched mode).
#' @param branch_exclude Character vector of cluster names to drop (optional; branched mode).
#' @param branch_min_cells Integer; require at least this many cells after subsetting (default 10).
#' @param drop_unused_levels Logical; drop unused factor levels in returned labels (default TRUE).
#'
#' @return A list of class "multi_dope_results" containing:
#' \item{results}{Named list where each element contains DOPE results for one trajectory}
#' \item{comparison_summary}{Data frame comparing all trajectories across metrics}
#' \item{best_trajectory}{Name of trajectory with highest composite DOPE score}
#' \item{n_trajectories}{Number of trajectories analyzed}
#' \item{method_info}{Analysis parameters used}
#'
#' @examples
#' # pseudotime_methods <- list("monocle3" = mono_pt, "slingshot" = sling_pt)
#' # multi_results <- compute_multi_DOPE(
#' #   expr_or_seurat = seu,
#' #   pseudotime_list = pseudotime_methods,
#' #   naive_markers = c("TCF7","LEF1","CCR7"),
#' #   terminal_markers = c("GZMB","PRF1","IFNG"),
#' #   E_method = "gmm",
#' #   trajectory = "branched",
#' #   cluster_labels = "Phenotype",
#' #   branch_include = c("Stem_Progenitors","Monocyte_progenitors","Monocytes")
#' # )
#'
#' @export
# ===============================
# DOPE: Multi-trajectory wrapper
# (linear or branched) — single file
# ===============================

# ---------- Small utilities ----------
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

# Detect whether a function supports an argument by name
.has_formal <- function(fun, arg) {
  is.function(fun) && !is.null(names(formals(fun))) && arg %in% names(formals(fun))
}

# ---------- Data access & alignment ----------
.get_expr <- function(x, assay = NULL, slot = "data") {
  if (inherits(x, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Seurat object provided but Seurat package is not available.")
    if (is.null(assay)) assay <- Seurat::DefaultAssay(x)
    return(as.matrix(Seurat::GetAssayData(x, assay = assay, slot = slot)))
  } else if (is.matrix(x) || inherits(x, "dgCMatrix")) {
    return(as.matrix(x))
  } else if (is.data.frame(x)) {
    return(as.matrix(x))
  } else {
    stop("expr_or_seurat must be a matrix/data.frame/dgCMatrix or a Seurat object")
  }
}

.align_to_cells <- function(x, target_cells, what = "vector") {
  if (!is.null(names(x))) {
    miss <- setdiff(target_cells, names(x))
    if (length(miss))
      stop("Named ", what, " is missing cells: ", paste(head(miss, 6), collapse = ", "),
           if (length(miss) > 6) " ...")
    x[target_cells]
  } else {
    if (length(x) != length(target_cells))
      stop("Length of ", what, " (", length(x), ") must match number of cells (", length(target_cells), "). ",
           "Provide names(x) to align by cell names.")
    x
  }
}

.per_cell_score <- function(expr_or_seurat, markers) {
  expr <- if (is.matrix(expr_or_seurat) || inherits(expr_or_seurat, "dgCMatrix") || is.data.frame(expr_or_seurat)) {
    as.matrix(expr_or_seurat)
  } else {
    .get_expr(expr_or_seurat)
  }
  genes <- intersect(markers, rownames(expr))
  if (length(genes) == 0) {
    warning("No markers found in expression matrix.")
    return(rep(NA_real_, ncol(expr)))
  }
  colSums(expr[genes, , drop = FALSE]) / length(genes)
}

# ---------- Simple branch subsetting helper ----------
# Works for Seurat or vector cluster labels; supports include/exclude, min_cells, level dropping
subset_by_clusters <- function(expr_or_seurat,
                               clusters,
                               include = NULL,
                               exclude = NULL,
                               min_cells = 10,
                               drop_unused_levels = TRUE,
                               assay = NULL) {
  if (inherits(expr_or_seurat, "Seurat")) {
    all_cells <- Seurat::Cells(expr_or_seurat)
  } else {
    x <- .get_expr(expr_or_seurat, assay = assay, slot = "data")
    all_cells <- colnames(x)
  }
  clusters <- .align_to_cells(clusters, all_cells, what = "cluster_labels")

  keep <- rep(TRUE, length(all_cells))
  if (!is.null(include)) keep <- clusters %in% include
  if (!is.null(exclude)) keep <- keep & !(clusters %in% exclude)
  idx <- which(keep)

  if (length(idx) < min_cells)
    stop("Branch subsetting left ", length(idx), " cells (< ", min_cells, ").")

  kept_cells <- all_cells[idx]
  if (inherits(expr_or_seurat, "Seurat")) {
    sub_obj <- subset(expr_or_seurat, cells = kept_cells)
    new_clusters <- sub_obj[[colnames(expr_or_seurat[[]])[match(TRUE, colnames(expr_or_seurat[[]]) == colnames(expr_or_seurat[[]])[1])]]] # dummy to avoid NSE
    # Safer: recompute from original clusters
    new_clusters <- clusters[kept_cells]
    if (drop_unused_levels) new_clusters <- as.character(factor(new_clusters))
    return(list(obj = sub_obj, cluster_labels = new_clusters))
  } else {
    mat <- .get_expr(expr_or_seurat, assay = assay, slot = "data")
    sub_mat <- mat[, kept_cells, drop = FALSE]
    new_clusters <- clusters[kept_cells]
    if (drop_unused_levels) new_clusters <- as.character(factor(new_clusters))
    return(list(obj = sub_mat, cluster_labels = new_clusters))
  }
}

# ===============================
# Core single-trajectory wrapper
# ===============================
#' Compute DOPE for a single trajectory (linear or branched)
#' @export
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
                                plot_E = TRUE,
                                verbose = FALSE,
                                trajectory = c("linear","branched"),
                                branch_include = NULL,
                                branch_exclude = NULL,
                                branch_min_cells = 10,
                                drop_unused_levels = TRUE,
                                tol = 1e-8) {
  E_method   <- match.arg(E_method)
  id_type    <- match.arg(id_type)
  trajectory <- match.arg(trajectory)

  expr <- .get_expr(expr_or_seurat, assay = assay, slot = slot)
  cells_all <- if (inherits(expr_or_seurat, "Seurat")) Seurat::Cells(expr_or_seurat) else colnames(expr)

  # align provided vectors to cell order
  pseudotime <- .align_to_cells(pseudotime, cells_all, what = "pseudotime")
  if (!is.null(cluster_labels)) {
    if (!(inherits(expr_or_seurat, "Seurat") && is.character(cluster_labels) && length(cluster_labels) == 1)) {
      cluster_labels <- .align_to_cells(cluster_labels, cells_all, what = "cluster_labels")
    }
  }

  # branched: subset first
  if (trajectory == "branched") {
    if (is.null(cluster_labels)) {
      stop("For trajectory='branched', provide 'cluster_labels' (Seurat meta column or per-cell vector).")
    }
    sub <- subset_by_clusters(
      expr_or_seurat = expr_or_seurat,
      clusters = if (inherits(expr_or_seurat, "Seurat") && is.character(cluster_labels) && length(cluster_labels) == 1)
        expr_or_seurat[[cluster_labels]][,1,drop=TRUE] else cluster_labels,
      include = branch_include,
      exclude = branch_exclude,
      min_cells = branch_min_cells,
      drop_unused_levels = drop_unused_levels,
      assay = assay
    )
    expr <- .get_expr(sub$obj, assay = assay, slot = slot)
    kept_cells <- if (inherits(sub$obj, "Seurat")) Seurat::Cells(sub$obj) else colnames(expr)
    pseudotime <- .align_to_cells(pseudotime, kept_cells, what = "pseudotime")
    cluster_labels <- sub$cluster_labels
  }

  if (length(pseudotime) != ncol(expr)) {
    stop("Length of pseudotime must match number of cells in expression data after alignment/subsetting.")
  }
  if (!is.null(cluster_labels) && length(cluster_labels) != ncol(expr)) {
    stop("Length of cluster_labels must match number of cells in expression data after alignment/subsetting.")
  }

  # ---- D metric
  D_res <- tryCatch({
    metrics_D(expr, naive_markers, terminal_markers, pseudotime)
  }, error = function(e) {
    if (verbose) message("D metric failed: ", e$message)
    list(D_naive = NA_real_, D_term = NA_real_, error = e$message)
  })

  # Precompute marker scores (E metric)
  naive_scores <- .per_cell_score(expr, naive_markers)
  term_scores  <- .per_cell_score(expr, terminal_markers)

  # ---- O metric (conditionally pass tol)
  O_res <- tryCatch({
    if (.has_formal(metrics_O, "tol")) {
      metrics_O(expr, naive_markers, terminal_markers, pseudotime, tol = tol)
    } else {
      metrics_O(expr, naive_markers, terminal_markers, pseudotime)
    }
  }, error = function(e) {
    if (verbose) message("O metric failed: ", e$message)
    list(O = NA_real_, error = e$message)
  })

  # ---- P metric
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

  # ---- E metric
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

  comp_vec <- c(
    D_res$D_naive %||% NA_real_,
    D_res$D_term  %||% NA_real_,
    O_res$O       %||% NA_real_,
    P_res$P       %||% NA_real_,
    E_res$E_comp  %||% NA_real_
  )
  DOPE_score <- if (all(is.na(comp_vec))) NA_real_ else mean(comp_vec, na.rm = TRUE)

  out <- list(
    D = list(D_naive = D_res$D_naive %||% NA_real_, D_term = D_res$D_term %||% NA_real_),
    O = list(O = O_res$O %||% NA_real_),
    P = list(P = P_res$P %||% NA_real_),
    E = list(E_naive = E_res$E_naive %||% NA_real_, E_term = E_res$E_term %||% NA_real_, E_comp = E_res$E_comp %||% NA_real_),
    DOPE_score = DOPE_score,
    errors = list(D_error = D_res$error, O_error = O_res$error, P_error = P_res$error, E_error = E_res$error)
  )
  class(out) <- "dope_results"
  out
}

# ===============================
# Multi-trajectory wrapper
# ===============================
#' DOPE: Comprehensive pseudotime trajectory quality assessment for multiple trajectories
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
                               plot_E = TRUE,
                               verbose = TRUE,
                               parallel = FALSE,
                               n_cores = NULL,
                               trajectory = c("linear","branched"),
                               branch_include = NULL,
                               branch_exclude = NULL,
                               branch_min_cells = 10,
                               drop_unused_levels = TRUE,
                               tol = 1e-8) {

  E_method   <- match.arg(E_method)
  id_type    <- match.arg(id_type)
  trajectory <- match.arg(trajectory)

  if (!is.list(pseudotime_list)) stop("pseudotime_list must be a list of pseudotime vectors")
  if (is.null(names(pseudotime_list))) names(pseudotime_list) <- paste0("trajectory_", seq_along(pseudotime_list))
  n_trajectories <- length(pseudotime_list)

  if (verbose) {
    cat("Computing DOPE metrics for", n_trajectories, "trajectories...\n")
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

  compute_single_dope <- function(trajectory_name, pseudotime_vector) {
    if (verbose) cat("\n--- Processing trajectory:", trajectory_name, "---\n")
    tryCatch({
      res <- compute_DOPE_single(
        expr_or_seurat = expr_or_seurat,
        pseudotime     = pseudotime_vector,
        naive_markers  = naive_markers,
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
        plot_E = plot_E,
        verbose = verbose,
        trajectory = trajectory,
        branch_include = branch_include,
        branch_exclude = branch_exclude,
        branch_min_cells = branch_min_cells,
        drop_unused_levels = drop_unused_levels,
        tol = tol
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
      varlist = c("compute_DOPE_single",
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
      n_cores = if (use_parallel) n_cores else NA_integer_,
      trajectory = trajectory,
      branch_include = branch_include,
      branch_exclude = branch_exclude,
      branch_min_cells = branch_min_cells
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
  multi_results
}

# ===============================
# Comparison summary / printing / plotting
# ===============================
#' Create comparison summary across trajectories
#' @export
create_comparison_summary <- function(trajectory_results) {
  if (length(trajectory_results) == 0) stop("trajectory_results cannot be empty")

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

  for (i in seq_along(trajectory_results)) {
    result <- trajectory_results[[i]]
    trajectory_name <- names(trajectory_results)[i] %||% paste0("Trajectory_", i)
    has_error <- if (!is.null(result$errors)) any(vapply(result$errors, function(x) !is.null(x), logical(1))) else FALSE

    row_data <- data.frame(
      trajectory = trajectory_name,
      D_naive = result$D$D_naive %||% NA_real_,
      D_term  = result$D$D_term  %||% NA_real_,
      O       = result$O$O       %||% NA_real_,
      P       = result$P$P       %||% NA_real_,
      E_naive = result$E$E_naive %||% NA_real_,
      E_term  = result$E$E_term  %||% NA_real_,
      E_comp  = result$E$E_comp  %||% NA_real_,
      DOPE_score = result$DOPE_score %||% NA_real_,
      has_error = has_error,
      stringsAsFactors = FALSE
    )
    summary_data <- rbind(summary_data, row_data)
  }

  metrics <- c("D_naive","D_term","O","P","E_naive","E_term","E_comp","DOPE_score")
  for (metric in metrics) {
    rank_col <- paste0(metric, "_rank")
    summary_data[[rank_col]] <- rank(-summary_data[[metric]], na.last = "keep", ties.method = "min")
  }

  summary_data <- summary_data[order(summary_data$DOPE_score, decreasing = TRUE, na.last = TRUE), ]
  rownames(summary_data) <- NULL
  summary_data
}

#' @export
print.multi_dope_results <- function(x, ...) {
  cat("Multi-Trajectory DOPE Analysis\n")
  cat("==============================\n\n")

  cat("Analysis info:\n")
  cat(sprintf("  Trajectories analyzed: %d\n", x$n_trajectories))
  cat(sprintf("  E method: %s\n", x$method_info$E_method))
  if (isTRUE(x$method_info$parallel)) cat(sprintf("  Parallel processing: %d cores\n", x$method_info$n_cores))
  cat(sprintf("  Trajectory mode: %s\n\n", x$method_info$trajectory))

  cat("Trajectory Comparison:\n")
  cat("---------------------\n")
  top_n <- min(5, nrow(x$comparison_summary))
  top_trajectories <- x$comparison_summary[1:top_n, ]

  for (i in seq_len(nrow(top_trajectories))) {
    traj <- top_trajectories[i, ]
    rank_suffix <- switch(as.character(i), "1"="st","2"="nd","3"="rd","th")
    dope_score <- if (is.na(traj$DOPE_score)) "NA" else sprintf("%.3f", traj$DOPE_score)
    cat(sprintf("%d%s: %s (DOPE: %s)\n", i, rank_suffix, traj$trajectory, dope_score))
  }
  if (!is.na(x$best_trajectory)) cat(sprintf("\nBest trajectory: %s\n", x$best_trajectory))
  cat("\nUse summary() for detailed comparison table\n")
}

#' @export
summary.multi_dope_results <- function(object, ...) {
  print(object)
  cat("\nDetailed Comparison Table:\n")
  cat("=========================\n")
  display_cols <- c("trajectory","DOPE_score","DOPE_score_rank","D_term","O","P","E_comp")
  display_cols <- display_cols[display_cols %in% names(object$comparison_summary)]
  display_data <- object$comparison_summary[, display_cols, drop = FALSE]
  numeric_cols <- vapply(display_data, is.numeric, logical(1))
  display_data[numeric_cols] <- lapply(display_data[numeric_cols], function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x)))
  print(display_data, row.names = FALSE)

  if (!is.na(object$best_trajectory)) {
    cat(sprintf("\nDetailed results for best trajectory (%s):\n", object$best_trajectory))
    cat("==========================================\n")
    best_result <- object$results[[object$best_trajectory]]
    if (inherits(best_result, "dope_results")) print(best_result)
  }
  invisible(object)
}

#' @export
plot.multi_dope_results <- function(multi_dope_results,
                                    metrics = c("D_naive","D_term","O","P","E_naive","E_term","DOPE_score"),
                                    type = "bar") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package required for plotting")
  plot_data <- multi_dope_results$comparison_summary[, c("trajectory", metrics), drop = FALSE]

  if (type == "bar") {
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("reshape2 package required for bar plots")
    plot_data_long <- reshape2::melt(plot_data, id.vars = "trajectory", variable.name = "metric", value.name = "score")
    return(
      ggplot2::ggplot(plot_data_long, ggplot2::aes(x = trajectory, y = score, fill = metric)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "DOPE Metrics Comparison Across Trajectories", x = "Trajectory Method", y = "Score") +
        ggplot2::ylim(0, 1)
    )
  } else if (type == "radar") {
    if (!requireNamespace("fmsb", quietly = TRUE)) stop("fmsb package required for radar plots. Install with: install.packages('fmsb')")
    radar_data <- plot_data[, metrics, drop = FALSE]
    radar_data[is.na(radar_data)] <- 0
    for (metric in metrics) if (metric %in% names(radar_data)) radar_data[[metric]] <- pmax(0, pmin(1, radar_data[[metric]]))
    radar_data <- rbind(rep(1, length(metrics)), rep(0, length(metrics)), radar_data)
    rownames(radar_data) <- c("Max","Min", plot_data$trajectory)
    n_trajectories <- nrow(plot_data)
    colors <- rainbow(n_trajectories, alpha = 0.3)
    line_colors <- rainbow(n_trajectories, alpha = 0.8)
    plot_radar <- function() {
      layout(matrix(c(1,2), ncol = 2), widths = c(3,1))
      fmsb::radarchart(radar_data, axistype = 1, pcol = line_colors, pfcol = colors, plwd = 2, plty = 1,
                       cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = rep("",5),
                       cglwd = 0.5, vlcex = 0, title = "DOPE Metrics Radar Chart")
      par(mar = c(0,0,0,0)); plot.new()
      legend("center", legend = plot_data$trajectory, col = line_colors, lty = 1, lwd = 2, cex = 0.8, bty = "n")
      layout(1)
    }
    plot_radar(); return(invisible(plot_radar))
  } else if (type == "heatmap") {
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("reshape2 package required for heatmap plots")
    heatmap_data <- reshape2::melt(plot_data, id.vars = "trajectory", variable.name = "metric", value.name = "score")
    return(
      ggplot2::ggplot(heatmap_data, ggplot2::aes(x = metric, y = trajectory, fill = score)) +
        ggplot2::geom_tile(color = "white", size = 0.5) +
        ggplot2::scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                                      midpoint = 0.5, name = "Score", na.value = "grey90") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title = "DOPE Metrics Heatmap", x = "Metric", y = "Trajectory Method") +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(score), "NA", sprintf("%.2f", score))),
                           color = "black", size = 3)
    )
  } else {
    stop("Plot type must be one of: 'bar', 'radar', or 'heatmap'")
  }
}
