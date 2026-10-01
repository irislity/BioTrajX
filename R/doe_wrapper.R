#' Compute DOE metrics for a single trajectory (linear)
#'
#' @name compute_single_doe_linear
#' @description
#' High-level convenience wrapper that computes the DOE components—
#' Directionality (D), Order consistency (O), and
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
#' @param early_markers Character vector of gene symbols/IDs expected to be high
#'   in the early state.
#' @param terminal_markers Character vector of gene symbols/IDs expected to be
#'   high in the terminal/late state.
#' @param cluster_labels Either (i) a per-cell vector (same length/order as
#'   cells) giving cluster/branch labels, or (ii) **for Seurat inputs only**, a
#'   single character string giving the name of a metadata column to use.
#' @param early_clusters Character vector of cluster names defining early
#'   clusters for `E` when `E_method %in% c("clusters","combined")`.
#' @param terminal_clusters Character vector of cluster names defining terminal
#'   clusters for `E` when `E_method %in% c("clusters","combined")`.
#' @param E_method One of `"gmm"`, `"clusters"`, or `"combined"`. Controls how
#'   endpoint validity is estimated inside `metrics_e()`.
#' @param assay Assay name to pull from a Seurat object (defaults to
#'   `Seurat::DefaultAssay(x)` when `NULL`).
#' @param slot Slot to pull from a Seurat assay; typically `"data"` (default),
#'   `"counts"`, or `"scale.data"`.
#' @param plot_E Logical; if `TRUE`, allow `metrics_e()` to produce diagnostic
#'   density plots when using the GMM mode (default `TRUE`). Only takes effect
#'   when `plot_metrics` is `TRUE`; `plot_metrics` gates all plotting, and
#'   `plot_E` fine-tunes this one plot within that.
#' @param plot_metrics Logical; if `TRUE`, draw the D/O/E diagnostic plots
#'   (`plot_metrics_d()`, `plot_metrics_o()`, `plot_metrics_e()`, and, when
#'   `plot_E` is also `TRUE`, the internal GMM density plot from
#'   `metrics_e()`) after computing the metrics. Defaults to `NULL`, which is
#'   treated as `FALSE` (no plotting at all).
#' @param verbose Logical; if `TRUE`, print diagnostic messages when component
#'   metrics fail and are caught (default `FALSE`).
#' @param drop_unused_levels Logical; if `TRUE`, drop unused factor levels in
#'   the subsetted cluster labels (default `TRUE`).
#' @param pseudotime_rescale Logical; if `TRUE` (default), min-max normalize
#'   `pseudotime` to \[0,1\] before computing D/O/E. Set to `FALSE` to use the
#'   raw pseudotime values as provided. D/O/E are each invariant to monotonic
#'   transforms of pseudotime, so this only affects the scale shown in
#'   diagnostic plots, not the computed scores.
#'
#' @details
#' The wrapper:
#' \itemize{
#'   \item Aligns `pseudotime` (and `cluster_labels` when provided) to the
#'         current cell order.
#'   \item Optionally subsets to a specified branch for `"branched"` trajectories.
#'   \item Computes D via `metrics_d()`, O via `metrics_o()`,
#'         and E via `metrics_e()` (controlled by `E_method`, `plot_E`, and
#'         optional cluster specification).
#'   \item Aggregates a scalar `DOE_score` as the mean of available component
#'         scores (ignoring `NA`s).
#' }
#'
#' Component failures are caught; the corresponding entry is set to `NA` and the
#' error message is recorded under `errors`. Per-cell early/terminal scores used
#' by `E` are computed internally as mean expression over the provided marker
#' sets.
#'
#' @return An object of class `"doe_results"`, a list with elements:
#' \describe{
#'   \item{$D$}{List with `D_early`, `D_term`.}
#'   \item{$O$}{List with scalar `O`.}
#'   \item{$E$}{List with `E_early`, `E_term`, and composite `E_comp`.}
#'   \item{$DOE_score$}{Scalar mean over available component scores.}
#'   \item{$errors$}{List of caught error messages for D/O/E (may be `NULL`).}
#' }
#'
#' @note
#' This wrapper relies on the availability and behavior of the helper functions
#' `metrics_d()`, `metrics_o()`, and `metrics_e()`, as well as
#' optional plotting utilities (e.g. `plot_metrics_*`). Ensure these are loaded
#' in your package namespace. If working with Seurat inputs, the **Seurat**
#' package must be installed and the requested assay/slot must exist.
#'
#' @seealso
#' `metrics_d()`, `metrics_o()`, `metrics_e()`,
#' `subset_by_clusters()`
#'
#' @examples
#' \dontrun{
#'
#' res <- compute_single_doe_linear(
#'   expr_or_seurat   = expr,
#'   pseudotime       = pt,
#'   early_markers    = early,
#'   terminal_markers = term,
#'   E_method         = "gmm",
#'   species          = "Homo sapiens"
#' )
#' str(res)
#'
#' }
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
    return(as.matrix(Seurat::GetAssayData(x, assay = assay, layer = slot)))
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
#' @rdname compute_single_doe_linear
#' @export
compute_single_doe_linear <- function(expr_or_seurat,
                                       pseudotime,
                                       early_markers,
                                       terminal_markers,
                                       cluster_labels = NULL,
                                       early_clusters = NULL,
                                       terminal_clusters = NULL,
                                       E_method = c("gmm", "clusters", "combined"),
                                       assay = NULL,
                                       slot = "data",
                                       plot_E = TRUE,
                                       plot_metrics = NULL,
                                       verbose = FALSE,
                                       branch_include = NULL,
                                       branch_exclude = NULL,
                                       branch_min_cells = 10,
                                       drop_unused_levels = TRUE,
                                       pseudotime_rescale = TRUE,
                                       tol = 1e-8) {
  E_method   <- match.arg(E_method)

  expr <- .get_expr(expr_or_seurat, assay = assay, slot = slot)
  cells_all <- if (inherits(expr_or_seurat, "Seurat")) Seurat::Cells(expr_or_seurat) else colnames(expr)

  # align provided vectors to cell order
  pseudotime <- .align_to_cells(pseudotime, cells_all, what = "pseudotime")
  if (!is.null(cluster_labels)) {
    if (inherits(expr_or_seurat, "Seurat") && is.character(cluster_labels) && length(cluster_labels) == 1) {
      # cluster_labels is a metadata column name — extract actual per-cell labels
      col_name <- cluster_labels
      if (!col_name %in% colnames(expr_or_seurat@meta.data))
        stop("cluster_labels column '", col_name, "' not found in Seurat metadata.")
      cluster_labels <- setNames(expr_or_seurat@meta.data[[col_name]], colnames(expr_or_seurat))[cells_all]
    } else {
      cluster_labels <- .align_to_cells(cluster_labels, cells_all, what = "cluster_labels")
    }
  }

  if (length(pseudotime) != ncol(expr)) {
    stop("Length of pseudotime must match number of cells in expression data after alignment/subsetting.")
  }
  if (!is.null(cluster_labels) && length(cluster_labels) != ncol(expr)) {
    stop("Length of cluster_labels must match number of cells in expression data after alignment/subsetting.")
  }

  # Some TI methods (e.g. TSCAN) only assign pseudotime to a subset of cells
  # (those on the main path through their trajectory graph), leaving the rest
  # NA. Drop those cells here — once, before any metric runs — rather than
  # letting NA propagate into metrics_o's range calculation or metrics_e's
  # hard NA check.
  valid <- !is.na(pseudotime)
  if (!all(valid)) {
    if (verbose)
      message(sprintf("  Dropping %d/%d cells with NA pseudotime before computing D/O/E",
                      sum(!valid), length(valid)))
    expr       <- expr[, valid, drop = FALSE]
    pseudotime <- pseudotime[valid]
    if (!is.null(cluster_labels)) cluster_labels <- cluster_labels[valid]
  }

  # Rescale pseudotime to [0,1] via min-max normalization so trajectories from
  # different TI methods (which can emit wildly different raw ranges) are
  # comparable, e.g. in the diagnostic plots below. D/O/E are each invariant
  # to monotonic transforms of pseudotime, so this never changes their scores.
  if (isTRUE(pseudotime_rescale)) {
    pseudotime <- .minmax_normalize(pseudotime, label = "pseudotime values")
  }

  # ---- D metric
  D_res <- tryCatch({
    metrics_d(expr, early_markers, terminal_markers, pseudotime)
  }, error = function(e) {
    if (verbose) message("D metric failed: ", e$message)
    list(D_early = NA_real_, D_term = NA_real_, error = e$message)
  })

  # Precompute marker scores (E metric)
  early_scores <- .per_cell_score(expr, early_markers)
  term_scores  <- .per_cell_score(expr, terminal_markers)

  # ---- O metric
  O_res <- tryCatch(
    metrics_o(
      expr = expr,
      early_markers = early_markers,
      terminal_markers = terminal_markers,
      pseudotime = pseudotime
    ),
    error = function(e) {
      if (verbose) message("O metric failed: ", e$message)
      list(
        O = NA_real_,
        O_f = numeric(0),
        genes_used = character(0),
        orientation = NA_character_,
        error = e$message
      )
    }
  )
  # ---- E metric
  # plot_metrics gates all plotting; plot_E only fine-tunes the internal GMM
  # density plot when plot_metrics is TRUE.
  E_res <- tryCatch({
    metrics_e(pseudotime,
              cluster_labels = cluster_labels,
              early_marker_scores = early_scores,
              terminal_marker_scores = term_scores,
              early_clusters = early_clusters,
              terminal_clusters = terminal_clusters,
              method = E_method,
              plot = isTRUE(plot_metrics) && isTRUE(plot_E))
  }, error = function(e) {
    if (verbose) message("E metric failed: ", e$message)
    list(E_early = NA_real_, E_term = NA_real_, E_comp = NA_real_, error = e$message)
  })

  D_comp <- mean(c(D_res$D_early %||% NA_real_, D_res$D_term %||% NA_real_), na.rm = TRUE)
  if (is.nan(D_comp)) D_comp <- NA_real_

  comp_vec <- c(D_comp, O_res$O %||% NA_real_, E_res$E_comp %||% NA_real_)
  DOE_score <- if (all(is.na(comp_vec))) NA_real_ else mean(comp_vec, na.rm = TRUE)

  out <- list(
    D = list(D_early = D_res$D_early %||% NA_real_, D_term = D_res$D_term %||% NA_real_,
             D_comp = D_comp),
    O = list(O = O_res$O %||% NA_real_, orientation = O_res$orientation %||% NA_character_),
    E = list(E_early = E_res$E_early %||% NA_real_, E_term = E_res$E_term %||% NA_real_,
             E_comp = E_res$E_comp %||% NA_real_),
    DOE_score = DOE_score,
    errors = list(D_error = D_res$error, O_error = O_res$error, E_error = E_res$error)
  )
  class(out) <- "doe_results"

  if (isTRUE(plot_metrics)) {
    plot_metrics_d(pseudotime,D_res)
    plot_metrics_o(expr,pseudotime,O_res,
                   early_markers,terminal_markers)
    plot_metrics_e(E_res)
  }

  return(out)

}




# ------create_comparison_summary---------
create_comparison_summary <- function(trajectory_results) {
  if (length(trajectory_results) == 0) stop("trajectory_results cannot be empty")

  summary_data <- data.frame(
    trajectory = character(),
    D_early = numeric(),
    D_term = numeric(),
    D_comp = numeric(),
    O = numeric(),
    O_orientation = character(),
    E_early = numeric(),
    E_term = numeric(),
    E_comp = numeric(),
    DOE_score = numeric(),
    has_error = logical(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(trajectory_results)) {
    result <- trajectory_results[[i]]
    trajectory_name <- names(trajectory_results)[i] %||% paste0("Trajectory_", i)
    has_error <- if (!is.null(result$errors)) any(vapply(result$errors, function(x) !is.null(x), logical(1))) else FALSE

    row_data <- data.frame(
      trajectory = trajectory_name,
      D_early = result$D$D_early %||% NA_real_,
      D_term  = result$D$D_term  %||% NA_real_,
      D_comp  = result$D$D_comp  %||% NA_real_,
      O       = result$O$O       %||% NA_real_,
      O_orientation = result$O$orientation %||% NA_character_,
      E_early = result$E$E_early %||% NA_real_,
      E_term  = result$E$E_term  %||% NA_real_,
      E_comp  = result$E$E_comp  %||% NA_real_,
      DOE_score = result$DOE_score %||% NA_real_,
      has_error = has_error,
      stringsAsFactors = FALSE
    )
    summary_data <- rbind(summary_data, row_data)
  }

  metrics <- c("D_early","D_term","D_comp","O","E_early","E_term","E_comp","DOE_score")
  for (metric in metrics) {
    rank_col <- paste0(metric, "_rank")
    summary_data[[rank_col]] <- rank(-summary_data[[metric]], na.last = "keep", ties.method = "min")
  }

  summary_data <- summary_data[order(summary_data$DOE_score, decreasing = TRUE, na.last = TRUE), ]
  rownames(summary_data) <- NULL
  summary_data
}



# ===============================
# Multi-trajectory wrapper
# ===============================
#' Compute DOE metrics for multiple linear trajectories
#'
#' @description
#' Runs the full DOE pipeline (Directionality **D**, Order consistency **O**,
#' Endpoints validity **E**, and the combined
#' `DOE_score`) across **multiple linear** pseudotime vectors sharing the
#' same expression object. Results are collected per trajectory and summarized
#' in a comparison table with the best trajectory highlighted.
#'
#' @details
#' This function is a convenience wrapper that repeatedly calls
#' [compute_single_doe_linear()] for each linear trajectory provided in
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
#' @param early_markers Character vector of gene symbols/IDs for early
#'   programs (used by D/O/E as relevant).
#' @param terminal_markers Character vector of gene symbols/IDs for terminal/late
#'   programs.
#' @param cluster_labels Optional character/factor vector of cluster labels
#'   (one per cell) used by E (when `E_method = "clusters"` or `"combined"`).
#' @param early_clusters Optional character vector of cluster names considered
#'   early (E Case 2 / cluster-based).
#' @param terminal_clusters Optional character vector of cluster names considered
#'   terminal (E Case 2 / cluster-based).
#' @param E_method Endpoints validity method. One of `"gmm"`, `"clusters"`,
#'   or `"combined"`. See [compute_single_doe_linear()] for details.
#' @param assay If `expr_or_seurat` is a Seurat object, the assay to use.
#'   If `NULL`, defaults to `Seurat::DefaultAssay()`.
#' @param slot If `expr_or_seurat` is a Seurat object, the assay slot to
#'   extract (e.g., `"data"`, `"counts"`). Default `"data"`.
#' @param plot_E Logical; if `TRUE`, produce diagnostic density plots when
#'   `E_method` uses GMM. Only takes effect when `plot_metrics` is `TRUE`.
#' @param plot_metrics Logical; if `TRUE`, draw the D/O/E diagnostic plots for
#'   each trajectory, including the `plot_E` GMM plot when also `TRUE` (see
#'   [compute_single_doe_linear()]). Defaults to `NULL`, which is treated as
#'   `FALSE` (no plotting at all).
#' @param verbose Logical; print progress.
#' @param parallel Logical; if `TRUE` and multiple trajectories are provided,
#'   use **parallel** workers.
#' @param n_cores Integer number of cores. If `NULL`, uses `detectCores()-1`.
#' @param drop_unused_levels Logical; drop unused factor levels (where relevant).
#' @param pseudotime_rescale Logical; if `TRUE` (default), min-max normalize
#'   each trajectory's pseudotime to \[0,1\] before computing D/O/E. See
#'   [compute_single_doe_linear()] for details.
#' @param tol Numerical tolerance for internal numerical checks.
#'
#' @return
#' An object of class `"multi_doe_results"` with components:
#' \itemize{
#'   \item \code{results}: named list of per-trajectory DOE results
#'         (each as returned by \code{compute_single_doe_linear()}).
#'   \item \code{comparison_summary}: \code{data.frame} summarizing
#'         D/O/E and \code{DOE_score} per trajectory.
#'   \item \code{best_trajectory}: character scalar with the name of the
#'         highest-scoring trajectory (or \code{NA} if none valid).
#'   \item \code{n_trajectories}: integer count.
#'   \item \code{method_info}: list of key settings used.
#' }
#'
#' @seealso [compute_single_doe_linear()], plotting helpers like
#'   \code{plot.multi_doe_results ()}.
#'
#' @examples
#' \dontrun{
#' # expr: genes x cells matrix; pt1, pt2: numeric pseudotime vectors (length = ncol(expr))
#' res <- compute_multi_doe_linear(
#'   expr_or_seurat = expr,
#'   pseudotime_list = list(linear_a = pt1, linear_b = pt2),
#'   early_markers = c("TCF7","LEF1"),
#'   terminal_markers = c("GZMB","PRF1"),
#'   E_method = "combined",
#'   parallel = TRUE
#' )
#' res$comparison_summary
#' res$best_trajectory
#' }
#'
#' @export

compute_multi_doe_linear <- function(
    expr_or_seurat,
    pseudotime_list,
    early_markers,
    terminal_markers,
    cluster_labels      = NULL,
    early_clusters      = NULL,
    terminal_clusters   = NULL,
    E_method            = c("gmm", "clusters", "combined"),
    assay               = NULL,
    slot                = "data",
    plot_E              = TRUE,
    plot_metrics        = NULL,
    verbose             = TRUE,
    parallel            = FALSE,
    n_cores             = NULL,
    branch_include      = NULL,
    branch_exclude      = NULL,
    branch_min_cells    = 10,
    drop_unused_levels  = TRUE,
    pseudotime_rescale  = TRUE,
    tol                 = 1e-8
) {
  E_method   <- match.arg(E_method)

  if (!is.list(pseudotime_list)) {
    stop("pseudotime_list must be a list of pseudotime vectors")
  }
  if (is.null(names(pseudotime_list))) {
    names(pseudotime_list) <- paste0("trajectory_", seq_along(pseudotime_list))
  }

  n_trajectories <- length(pseudotime_list)

  if (verbose) {
    cat("Computing DOE metrics for", n_trajectories, "trajectories...\n")
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
  compute_single_doe <- function(trajectory_name, pseudotime_vector) {
    if (verbose) cat("\n--- Processing trajectory:", trajectory_name, "---\n")

    tryCatch(
      {
        res <- compute_single_doe_linear(
          expr_or_seurat   = expr_or_seurat,
          pseudotime       = pseudotime_vector,
          early_markers    = early_markers,
          terminal_markers = terminal_markers,
          cluster_labels   = cluster_labels,
          early_clusters   = early_clusters,
          terminal_clusters = terminal_clusters,
          E_method         = E_method,
          assay            = assay,
          slot             = slot,
          plot_E           = plot_E,
          plot_metrics     = plot_metrics,
          verbose          = verbose,
          branch_include   = branch_include,
          branch_exclude   = branch_exclude,
          branch_min_cells = branch_min_cells,
          drop_unused_levels = drop_unused_levels,
          pseudotime_rescale = pseudotime_rescale,
          tol               = tol
        )
        res$trajectory_name <- trajectory_name
        res
      },
      error = function(e) {
        warning(paste("Error processing trajectory", trajectory_name, ":", e$message))
        list(
          trajectory_name = trajectory_name,
          D = list(D_early = NA, D_term = NA),
          O = list(O = NA, orientation = NA_character_),
          E = list(E_early = NA, E_term = NA, E_comp = NA),
          DOE_score = NA,
          errors = list(
            D_error = e$message,
            O_error = e$message,
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
        "compute_single_doe", "compute_single_doe_linear",
        "metrics_d", "metrics_o", "metrics_e",
        ".get_expr", ".per_cell_score", ".align_to_cells", "subset_by_clusters",
        "%||%", ".has_formal", ".minmax_normalize", "tol"
      ),
      envir = environment()
    )
    parallel::clusterExport(
      cl,
      varlist = c("metrics_d", "metrics_o", "metrics_e",
                  ".get_expr", ".per_cell_score", ".align_to_cells",
                  "subset_by_clusters", "%||%", ".has_formal", ".minmax_normalize"),
      envir = asNamespace("BioTrajX")
    )

    idxs <- seq_along(pseudotime_list)
    trajectory_results <- parallel::parLapply(
      cl, idxs,
      function(i) {
        nm <- names(pseudotime_list)[i]
        compute_single_doe(nm, pseudotime_list[[i]])
      }
    )
    names(trajectory_results) <- names(pseudotime_list)

  } else {
    trajectory_results <- Map(
      compute_single_doe,
      names(pseudotime_list),
      pseudotime_list
    )
    names(trajectory_results) <- names(pseudotime_list)
  }

  # ---- Summarize & pick best ------------------------------------------------
  comparison_summary <- create_comparison_summary(trajectory_results)

  valid_scores <- comparison_summary$DOE_score[!is.na(comparison_summary$DOE_score)]
  best_trajectory <- if (length(valid_scores) > 0) {
    comparison_summary$trajectory[which.max(comparison_summary$DOE_score)]
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
      parallel         = use_parallel,
      n_cores          = if (use_parallel) n_cores else NA_integer_,
      branch_include   = branch_include,
      branch_exclude   = branch_exclude,
      branch_min_cells = branch_min_cells
    )
  )
  class(multi_results) <- "multi_doe_results"

  if (verbose) {
    cat("\n==========================================================\n")
    cat("Multi-trajectory DOE analysis complete!\n")
    if (!is.na(best_trajectory)) {
      best_score <- comparison_summary$DOE_score[
        comparison_summary$trajectory == best_trajectory
      ]
      cat(sprintf("Best trajectory: %s (DOE score: %.3f)\n",
                  best_trajectory, best_score))
    }
  }

  multi_results
}





#' Plot DOE Metrics Comparison Across Trajectories
#'
#' Generate visualizations of DOE metrics from a multi-trajectory comparison.
#' Supports bar plots, radar charts, and heatmaps for intuitive comparison
#' across multiple trajectory inference methods.
#'
#' @param multi_doe_results A list-like object containing a `comparison_summary`
#'   data frame. Must include a column named `"trajectory"` and columns for the
#'   specified metrics.
#' @param metrics Character vector of metric names to plot. Defaults to
#'   \code{c("D_early","D_term","O","E_early","E_term","DOE_score")}.
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
#' plot.multi_doe_results (multi_doe_results, type = "bar")
#' plot.multi_doe_results (multi_doe_results, type = "radar")
#' plot.multi_doe_results (multi_doe_results, type = "heatmap")
#' }
#' @method plot multi_doe_results
#' @export
plot.multi_doe_results <- function(multi_doe_results,
                                    metrics = c("D_early","D_term","O","E_early","E_term","DOE_score"),
                                    type = "bar") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package required for plotting")
  plot_data <- multi_doe_results$comparison_summary[, c("trajectory", metrics), drop = FALSE]

  if (type == "bar") {
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("reshape2 package required for bar plots")
    plot_data_long <- reshape2::melt(plot_data, id.vars = "trajectory", variable.name = "metric", value.name = "score")
    return(
      ggplot2::ggplot(plot_data_long, ggplot2::aes(x = trajectory, y = score, fill = metric)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "DOE Metrics Comparison Across Trajectories", x = "Trajectory Method", y = "Score") +
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
                       cglcol = "grey80", cglty = 1, axislabcol = "grey40",
                       caxislabels = c("0","0.25","0.5","0.75","1"), calcex = 0.72,
                       cglwd = 0.6, vlcex = 0.9, title = "DOE Metrics")
      par(mar = c(0,0,0,0)); plot.new()
      legend("center", legend = plot_data$trajectory, col = line_colors, lty = 1, lwd = 2, cex = 0.8, bty = "n")
      layout(1)
    }
    plot_radar(); return(invisible(plot_radar))
  } else if (type == "heatmap") {
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("reshape2 package required for heatmap plots")

    # Full column order: sub-components + composites grouped by letter
    heatmap_cols <- c("D_early","D_term","O","E_early","E_term","DOE_score")
    avail_cols   <- intersect(heatmap_cols,
                              colnames(multi_doe_results$comparison_summary))
    hm_data <- multi_doe_results$comparison_summary[, c("trajectory", avail_cols), drop = FALSE]

    composite_metrics <- c("DOE_score")

    heatmap_data <- reshape2::melt(hm_data, id.vars = "trajectory",
                                   variable.name = "metric", value.name = "score")
    orient_map <- setNames(
      multi_doe_results$comparison_summary$O_orientation,
      multi_doe_results$comparison_summary$trajectory
    )
    heatmap_data$cell_label <- ifelse(
      is.na(heatmap_data$score), "NA",
      ifelse(
        heatmap_data$metric == "O" & !is.na(orient_map[as.character(heatmap_data$trajectory)]),
        sprintf("%.2f\n(%s)", heatmap_data$score, orient_map[as.character(heatmap_data$trajectory)]),
        sprintf("%.2f", heatmap_data$score)
      )
    )
    # Rank trajectories by DOE_score: highest at top of plot (last factor level)
    doe_rank <- hm_data$trajectory[order(hm_data$DOE_score, decreasing = FALSE, na.last = FALSE)]
    heatmap_data$metric      <- factor(heatmap_data$metric, levels = avail_cols)
    heatmap_data$trajectory  <- factor(heatmap_data$trajectory, levels = doe_rank)
    heatmap_data$is_composite <- heatmap_data$metric %in% composite_metrics

    ggplot2::ggplot(heatmap_data,
                    ggplot2::aes(x = metric, y = trajectory, fill = score)) +
      ggplot2::geom_tile(ggplot2::aes(color = is_composite, linewidth = is_composite)) +
      ggplot2::scale_color_manual(values = c(`TRUE` = "#555555", `FALSE` = "white"),
                                  guide = "none") +
      ggplot2::scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4),
                                      guide = "none") +
      ggplot2::scale_fill_gradient2(low = "#d95f4f", mid = "#f3e96b", high = "#5fbf7a",
                                    midpoint = 0.5, limits = c(0, 1), name = "Score",
                                    na.value = "grey90") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
                       angle = 45, hjust = 1,
                       face  = ifelse(avail_cols %in% composite_metrics, "bold", "plain")),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = "DOE Metrics Heatmap", x = "Metric", y = "Trajectory Method") +
      ggplot2::geom_text(ggplot2::aes(label = cell_label),
                         color = "black", size = 3, lineheight = 0.9)
  } else {
    stop("Plot type must be one of: 'bar', 'radar', or 'heatmap'")
  }
}

#' Plot DOE metrics for a single linear trajectory
#'
#' Visualize the Directionality (D), Order (O), and Endpoint (E)
#' component metrics — along with the combined DOE score —
#' from a single trajectory analysis.
#'
#' @param doe_results A result object returned by
#'   [compute_single_doe_linear()], typically containing elements
#'   `$D`, `$O`, `$E`, and `$DOE_score`.
#' @param type Character string indicating plot type:
#'   one of `"bar"`, `"radar"`, or `"heatmap"`.
#'
#' @return
#' - For `"bar"` or `"heatmap"`: a \pkg{ggplot2} object.
#' - For `"radar"`: draws the plot and invisibly returns `NULL`.
#'
#' @details
#' This function summarizes all DOE components into an easy-to-read visual.
#' Each metric is normalized between 0 and 1 where applicable.
#' Missing (`NA`) values are handled gracefully.
#'
#' @examples
#' \dontrun{
#' res <- compute_single_doe_linear(expr, pt, early, term)
#' plot.doe_results(res, type = "bar")
#' plot.doe_results(res, type = "radar")
#' plot.doe_results(res, type = "heatmap")
#' }
#'
#' @method plot doe_results
#' @export
plot.doe_results <- function(doe_results, type = c("bar","radar","heatmap")) {
  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required for plotting")

  metrics <- data.frame(
    metric = c("D_early","D_term","O","E_early","E_term","DOE_score"),
    value = c(
      doe_results$D$D_early %||% NA_real_,
      doe_results$D$D_term  %||% NA_real_,
      doe_results$O$O       %||% NA_real_,
      doe_results$E$E_early %||% NA_real_,
      doe_results$E$E_term  %||% NA_real_,
      doe_results$DOE_score %||% NA_real_
    )
  )

  if (type == "bar") {
    p <- ggplot2::ggplot(metrics, ggplot2::aes(x = metric, y = value, fill = metric)) +
      ggplot2::geom_col() +
      ggplot2::theme_minimal() +
      ggplot2::coord_cartesian(ylim = c(0,1)) +
      ggplot2::labs(title = "Single Trajectory DOE Metrics",
                    x = "Metric", y = "Score") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    return(p)
  }

  if (type == "heatmap") {
    if (!requireNamespace("reshape2", quietly = TRUE))
      stop("reshape2 required for heatmap")
    metrics$value <- pmax(0, pmin(1, metrics$value))
    p <- ggplot2::ggplot(metrics, ggplot2::aes(x = 1, y = metric, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradient2(low = "#d95f4f", mid = "#f3e96b", high = "#5fbf7a",
                                    midpoint = 0.5, name = "Score", na.value = "grey90") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = "Single Trajectory DOE Heatmap", x = NULL, y = "Metric") +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(value), "NA", sprintf("%.2f", value))),
                         color = "black", size = 3)
    return(p)
  }

  if (type == "radar") {
    if (!requireNamespace("fmsb", quietly = TRUE))
      stop("fmsb required for radar plots: install.packages('fmsb')")
    if (!requireNamespace("scales", quietly = TRUE))
      stop("scales required for radar plots")

    radar_data <- as.data.frame(t(metrics$value))
    colnames(radar_data) <- metrics$metric
    radar_data[is.na(radar_data)] <- 0
    radar_data <- rbind(rep(1, ncol(radar_data)), rep(0, ncol(radar_data)), radar_data)
    rownames(radar_data) <- c("Max","Min","Trajectory")

    fmsb::radarchart(radar_data, axistype = 1,
                     pcol = "blue", pfcol = scales::alpha("skyblue", 0.4),
                     plwd = 2, plty = 1,
                     cglcol = "grey", cglty = 1, axislabcol = "grey",
                     cglwd = 0.5, vlcex = 0.8, caxislabels = rep("", 5),
                     title = "DOE Metrics Radar (Single Trajectory)")
    invisible(NULL)
  }
}


