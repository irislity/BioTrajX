#' Compute DOPE metrics for a single trajectory (linear)
#'
#' @description
#' High-level convenience wrapper that computes the four DOPE components—
#' Directionality (D), Order consistency (O), Program coherence (P), and
#' Endpoints validity (E)—for a single trajectory. Works with either a
#' **Seurat** object (pulling assay data via `GetAssayData`) or an expression
#' matrix-like input (dense `matrix`, `dgCMatrix`, or `data.frame` with genes in
#' rows and cells in columns). Optionally handles **branched** trajectories by
#' subsetting a specified branch before scoring.
#'
#' @param expr_or_seurat A Seurat object, numeric matrix, `dgCMatrix`, or
#'   `data.frame` of gene expression (genes x cells).
#' @param pseudotime Numeric vector of length equal to the number of cells after
#'   any subsetting; names should be cell IDs. If unnamed, it must be in the
#'   same order as the columns of `expr_or_seurat`.
#' @param naive_markers Character vector of gene symbols/IDs expected to be high
#'   in the naïve/early state.
#' @param terminal_markers Character vector of gene symbols/IDs expected to be
#'   high in the terminal/late state.
#' @param cluster_labels Either (i) a per-cell vector (same length/order as
#'   cells) giving cluster/branch labels, or (ii) **for Seurat inputs only**, a
#'   single character string giving the name of a metadata column to use.
#' @param naive_clusters Character vector of cluster names defining naïve
#'   clusters for `E` when `E_method %in% c("clusters","combined")`.
#' @param terminal_clusters Character vector of cluster names defining terminal
#'   clusters for `E` when `E_method %in% c("clusters","combined")`.
#' @param pathways Optional pathway resource passed to `metrics_P()`; typically a
#'   list or database handle understood by that function.
#' @param gene_sets Optional named list of gene sets (each a character vector)
#'   used by `metrics_P()` to compute program coherence.
#' @param E_method One of `"gmm"`, `"clusters"`, or `"combined"`. Controls how
#'   endpoint validity is estimated inside `metrics_E()`.
#' @param id_type Identifier type for genes used in pathway/gene-set lookups;
#'   one of `"SYMBOL"`, `"ENSEMBL"`, or `"ENTREZID"`. Default `"SYMBOL"`.
#' @param species Species string for pathway/gene-set lookup (default
#'   `"Homo sapiens"`).
#' @param assay Assay name to pull from a Seurat object (defaults to
#'   `Seurat::DefaultAssay(x)` when `NULL`).
#' @param slot Slot to pull from a Seurat assay; typically `"data"` (default),
#'   `"counts"`, or `"scale.data"`.
#' @param min_remaining Minimum number of genes remaining in a module/gene set
#'   after filtering inside `metrics_P()` (default `10`).
#' @param min_fraction Minimum fraction of genes from a module that must remain
#'   after filtering in `metrics_P()` (default `0.20`).
#' @param min_genes_per_module Minimum size of a gene set/module considered by
#'   `metrics_P()` (default `3`).
#' @param plot_E Logical; if `TRUE`, allow `metrics_E()` to produce diagnostic
#'   density plots when using the GMM mode (default `TRUE`).
#' @param verbose Logical; if `TRUE`, print diagnostic messages when component
#'   metrics fail and are caught (default `FALSE`).
#' @param drop_unused_levels Logical; if `TRUE`, drop unused factor levels in
#'   the subsetted cluster labels (default `TRUE`).
#' @param tol Numeric tolerance forwarded to `metrics_O()` when that function
#'   supports a `tol` argument (default `1e-8`).
#'
#' @details
#' The wrapper:
#' \itemize{
#'   \item Aligns `pseudotime` (and `cluster_labels` when provided) to the
#'         current cell order.
#'   \item Optionally subsets to a specified branch for `"branched"` trajectories.
#'   \item Computes D via `metrics_D()`, O via `metrics_O()` (passing `tol`
#'         when supported), P via `metrics_P()` (using `pathways`/`gene_sets`),
#'         and E via `metrics_E()` (controlled by `E_method`, `plot_E`, and
#'         optional cluster specification).
#'   \item Aggregates a scalar `DOPE_score` as the mean of available component
#'         scores (ignoring `NA`s).
#' }
#'
#' Component failures are caught; the corresponding entry is set to `NA` and the
#' error message is recorded under `errors`. Per-cell naïve/terminal scores used
#' by `E` are computed internally as mean expression over the provided marker
#' sets.
#'
#' @return An object of class `"dope_results"`, a list with elements:
#' \describe{
#'   \item{$D$}{List with `D_naive`, `D_term`.}
#'   \item{$O$}{List with scalar `O`.}
#'   \item{$P$}{List with scalar `P`.}
#'   \item{$E$}{List with `E_naive`, `E_term`, and composite `E_comp`.}
#'   \item{$DOPE_score$}{Scalar mean over available component scores.}
#'   \item{$errors$}{List of caught error messages for D/O/P/E (may be `NULL`).}
#' }
#'
#' @note
#' This wrapper relies on the availability and behavior of the helper functions
#' `metrics_D()`, `metrics_O()`, `metrics_P()`, and `metrics_E()`, as well as
#' optional plotting utilities (e.g. `plot_metrics_*`). Ensure these are loaded
#' in your package namespace. If working with Seurat inputs, the **Seurat**
#' package must be installed and the requested assay/slot must exist.
#'
#' @seealso
#' `metrics_D()`, `metrics_O()`, `metrics_P()`, `metrics_E()`,
#' `subset_by_clusters()`
#'
#' @examples
#' \dontrun{
#'
#' res <- compute_single_DOPE_linear(
#'   expr_or_seurat   = expr,
#'   pseudotime       = pt,
#'   naive_markers    = naive,
#'   terminal_markers = term,
#'   E_method         = "gmm",
#'   species          = "Homo sapiens"
#' )
#' str(res)
#'
#' }
#'
#' @export


# ---------- Small utilities ----------.
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
compute_single_DOPE_linear <- function(expr_or_seurat,
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
                                       branch_include = NULL,
                                       branch_exclude = NULL,
                                       branch_min_cells = 10,
                                       drop_unused_levels = TRUE,
                                       tol = 1e-8) {
  E_method   <- match.arg(E_method)
  id_type    <- match.arg(id_type)

  expr <- .get_expr(expr_or_seurat, assay = assay, slot = slot)
  cells_all <- if (inherits(expr_or_seurat, "Seurat")) Seurat::Cells(expr_or_seurat) else colnames(expr)

  # align provided vectors to cell order
  pseudotime <- .align_to_cells(pseudotime, cells_all, what = "pseudotime")
  if (!is.null(cluster_labels)) {
    if (!(inherits(expr_or_seurat, "Seurat") && is.character(cluster_labels) && length(cluster_labels) == 1)) {
      cluster_labels <- .align_to_cells(cluster_labels, cells_all, what = "cluster_labels")
    }
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


  plot_metrics_D(pseudotime,D_res)
  plot_metrics_O(expr,pseudotime,O_res,
                 naive_markers,terminal_markers)
  plot_metrics_P(expr,pseudotime, P_res)
  plot_metrics_E(E_res)

  return (out = out)

}




# ------create_comparison_summary---------
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

compute_multi_DOPE_linear <- function(
    expr_or_seurat,
    pseudotime_list,
    naive_markers,
    terminal_markers,
    cluster_labels      = NULL,
    naive_clusters      = NULL,
    terminal_clusters   = NULL,
    pathways            = NULL,
    gene_sets           = NULL,
    E_method            = c("gmm", "clusters", "combined"),
    id_type             = c("SYMBOL", "ENSEMBL", "ENTREZID"),
    species             = "Homo sapiens",
    assay               = NULL,
    slot                = "data",
    min_remaining       = 10,
    min_fraction        = 0.20,
    min_genes_per_module = 3,
    plot_E              = TRUE,
    verbose             = TRUE,
    parallel            = FALSE,
    n_cores             = NULL,
    branch_include      = NULL,
    branch_exclude      = NULL,
    branch_min_cells    = 10,
    drop_unused_levels  = TRUE,
    tol                 = 1e-8
) {
  E_method   <- match.arg(E_method)
  id_type    <- match.arg(id_type)

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

  # ---- Parallel setup -------------------------------------------------------
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

  # ---- Worker: single trajectory -------------------------------------------
  compute_single_dope <- function(trajectory_name, pseudotime_vector) {
    if (verbose) cat("\n--- Processing trajectory:", trajectory_name, "---\n")

    tryCatch(
      {
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
          branch_include   = branch_include,
          branch_exclude   = branch_exclude,
          branch_min_cells = branch_min_cells,
          drop_unused_levels = drop_unused_levels,
          tol              = tol
        )
        res$trajectory_name <- trajectory_name
        res
      },
      error = function(e) {
        warning(paste("Error processing trajectory", trajectory_name, ":", e$message))
        list(
          trajectory_name = trajectory_name,
          D = list(D_naive = NA, D_term = NA),
          O = list(O = NA),
          P = list(P = NA),
          E = list(E_naive = NA, E_term = NA, E_comp = NA),
          DOPE_score = NA,
          errors = list(
            D_error = e$message,
            O_error = e$message,
            P_error = e$message,
            E_error = e$message
          )
        )
      }
    )
  }

  # ---- Run over trajectories ------------------------------------------------
  if (use_parallel) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "compute_single_dope", "metrics_D", "metrics_O", "metrics_P", "metrics_E",
        ".get_expr", ".per_cell_score", ".align_to_cells", "subset_by_clusters",
        "%||%", ".has_formal"
      ),
      envir = environment()
    )
    # Load packages on workers if needed, e.g.:
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
    trajectory_results <- Map(
      compute_single_dope,
      names(pseudotime_list),
      pseudotime_list
    )
    names(trajectory_results) <- names(pseudotime_list)
  }

  # ---- Summarize & pick best ------------------------------------------------
  comparison_summary <- create_comparison_summary(trajectory_results)

  valid_scores <- comparison_summary$DOPE_score[!is.na(comparison_summary$DOPE_score)]
  best_trajectory <- if (length(valid_scores) > 0) {
    comparison_summary$trajectory[which.max(comparison_summary$DOPE_score)]
  } else {
    NA_character_
  }

  multi_results <- list(
    results           = trajectory_results,
    comparison_summary = comparison_summary,
    best_trajectory   = best_trajectory,
    n_trajectories    = n_trajectories,
    method_info = list(
      E_method         = E_method,
      id_type          = id_type,
      species          = species,
      parallel         = use_parallel,
      n_cores          = if (use_parallel) n_cores else NA_integer_,
      branch_include   = branch_include,
      branch_exclude   = branch_exclude,
      branch_min_cells = branch_min_cells
    )
  )
  class(multi_results) <- "multi_dope_results"

  if (verbose) {
    cat("\n==========================================================\n")
    cat("Multi-trajectory DOPE analysis complete!\n")
    if (!is.na(best_trajectory)) {
      best_score <- comparison_summary$DOPE_score[
        comparison_summary$trajectory == best_trajectory
      ]
      cat(sprintf("Best trajectory: %s (DOPE score: %.3f)\n",
                  best_trajectory, best_score))
    }
  }

  multi_results
}





#' Plot DOPE Metrics Comparison Across Trajectories
#'
#' Generate visualizations of DOPE metrics from a multi-trajectory comparison.
#' Supports bar plots, radar charts, and heatmaps for intuitive comparison
#' across multiple trajectory inference methods.
#'
#' @param multi_dope_results A list-like object containing a `comparison_summary`
#'   data frame. Must include a column named `"trajectory"` and columns for the
#'   specified metrics.
#' @param metrics Character vector of metric names to plot. Defaults to
#'   \code{c("D_naive","D_term","O","P","E_naive","E_term","DOPE_score")}.
#' @param type Character string indicating the plot type. One of:
#'   \code{"bar"} (default), \code{"radar"}, or \code{"heatmap"}.
#'
#' @details
#' - **Bar plots**: show metric scores per trajectory using grouped bars.
#'   Requires the \pkg{reshape2} package.
#' - **Radar charts**: display metrics in a radial layout for each trajectory.
#'   Requires the \pkg{fmsb} package.
#' - **Heatmaps**: present metric values in a matrix form with annotated scores.
#'   Requires the \pkg{reshape2} package.
#'
#' Values are automatically constrained between 0 and 1 for radar charts.
#' Missing values are displayed as 0 in radar plots and labeled as "NA" in heatmaps.
#'
#' @return
#' - For \code{type = "bar"} or \code{"heatmap"}: a \pkg{ggplot2} object.
#' - For \code{type = "radar"}: draws the plot and invisibly returns the plotting function.
#'
#' @examples
#' \dontrun{
#' # Example with comparison_summary data
#' plot.multi_dope_results(multi_dope_results, type = "bar")
#' plot.multi_dope_results(multi_dope_results, type = "radar")
#' plot.multi_dope_results(multi_dope_results, type = "heatmap")
#' }
#'
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
