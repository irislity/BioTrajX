# =========================
# Branched DOE
# =========================

# ---------- Utilities ----------
.stop_if_missing <- function(lst, required_names, label) {
  miss <- setdiff(required_names, names(lst))
  if (length(miss)) stop("Missing names in ", label, ": ", paste(miss, collapse = ", "))
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

# Detect whether a function supports an argument by name
.has_formal <- function(fun, arg) {
  is.function(fun) && !is.null(names(formals(fun))) && arg %in% names(formals(fun))
}

# ---------- Row flattener for branched/linear ----------
.flatten_traj_result <- function(result, traj_name) {
  has_branches <- !is.null(result$branches) && length(result$branches) > 0

  pull_row <- function(res, branch = NA_character_) {
    data.frame(
      trajectory   = traj_name,
      branch_data       = branch,
      D_naive      = res$D$D_naive %||% NA_real_,
      D_term       = res$D$D_term  %||% NA_real_,
      O            = res$O$O       %||% NA_real_,
      E_naive      = res$E$E_naive %||% NA_real_,
      E_term       = res$E$E_term  %||% NA_real_,
      E_comp       = res$E$E_comp  %||% NA_real_,
      DOE_score   = res$DOE_score %||% NA_real_,
      has_error    = if (!is.null(res$errors)) any(vapply(res$errors, function(x) !is.null(x), logical(1))) else FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (!has_branches) {
    pull_row(result, NA_character_)
  } else {
    br_names <- names(result$branches)
    do.call(rbind, lapply(br_names, function(b) pull_row(result$branches[[b]], branch = b)))
  }
}

.add_metric_ranks <- function(df, metric_cols) {
  for (metric in metric_cols) {
    rank_col <- paste0(metric, "_rank")
    df[[rank_col]] <- rank(-df[[metric]], na.last = "keep", ties.method = "min")
  }
  df
}

# ---------- Core: single trajectory across multiple branches ----------
#' Compute DOE metrics for a branched trajectory
#'
#' Evaluate Directionality, Order, and Endpoint (DOE) metrics for a single
#' pseudotime trajectory that contains multiple branches. Marker sets are
#' provided per branch and the function optionally subsets cells before
#' computing the component metrics.
#'
#' @param expr_or_seurat Expression matrix (genes C cells) or Seurat object.
#' @param pseudotime Numeric vector of pseudotime values (named by cell).
#' @param naive_markers_list Named list of naC/ve-state markers per branch.
#' @param terminal_markers_list Named list of terminal-state markers per branch.
#' @param cluster_labels Optional per-cell cluster/branch labels or metadata
#'   column name (for Seurat inputs).
#' @param branch_filters Optional named list of branch-specific subsetting
#'   rules (each element can contain `include`, `exclude`, `min_cells`,
#'   `drop_unused_levels`).
#' @param E_method Endpoint scoring approach; one of `"gmm"`, `"clusters"`,
#'   or `"combined"`.
#' @param assay Assay to use when `expr_or_seurat` is a Seurat object.
#' @param slot Assay slot to use when `expr_or_seurat` is a Seurat object.
#' @param plot_E Logical; pass through to `metrics_E()` plotting helpers.
#' @param verbose Logical; print progress messages.
#' @param tol Numerical tolerance forwarded to `metrics_O()` when supported.
#'
#' @return A list containing branch-level DOE results, aggregated scores, and
#'   metadata. Each branch entry includes component metrics and a scalar
#'   `DOE_score`.
#'
#' @export
compute_single_DOE_branched <- function(expr_or_seurat,
                                         pseudotime,
                                         naive_markers_list,
                                         terminal_markers_list,
                                         cluster_labels = NULL,
                                         branch_filters = NULL,
                                         E_method = c("gmm","clusters","combined"),
                                         assay = NULL,
                                         slot = "data",
                                         plot_E = FALSE,
                                         verbose = TRUE,
                                         tol = 1e-8) {
  E_method <- match.arg(E_method)

  # ---- Branch sanity ----
  branch_names <- union(names(naive_markers_list), names(terminal_markers_list))
  if (is.null(branch_names) || any(branch_names == "")) {
    stop("Branch lists must be named.")
  }
  .stop_if_missing(naive_markers_list,    branch_names, "naive_markers_list")
  .stop_if_missing(terminal_markers_list, branch_names, "terminal_markers_list")
  # ---- Expression extraction ----
  is_seurat <- inherits(expr_or_seurat, "Seurat")
  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Seurat object provided but Seurat package is not available.")
    if (is.null(assay)) assay <- Seurat::DefaultAssay(expr_or_seurat)
    expr <- as.matrix(Seurat::GetAssayData(expr_or_seurat, assay = assay, slot = slot))
    cells_all <- Seurat::Cells(expr_or_seurat)
  } else {
    expr <- as.matrix(expr_or_seurat)
    cells_all <- colnames(expr)
  }
  if (is.null(cells_all) || length(cells_all) == 0)
    stop("No column (cell) names found in expression matrix.")

  # ---- Align helpers ----
  .align_to_cells <- function(x, target_cells, what="vector") {
    if (!is.null(names(x))) {
      miss <- setdiff(target_cells, names(x))
      if (length(miss)) {
        stop("Named ", what, " missing cells: ", paste(head(miss, 5), collapse = ", "),
             if (length(miss) > 5) " ...")
      }
      x[target_cells]
    } else {
      if (length(x) != length(target_cells))
        stop("Length of ", what, " (", length(x), ") must match #cells (", length(target_cells), ").")
      x
    }
  }

  pt_all <- .align_to_cells(pseudotime, cells_all, "pseudotime")

  if (!is.null(cluster_labels)) {
    if (is_seurat && is.character(cluster_labels) && length(cluster_labels) == 1) {
      if (!cluster_labels %in% colnames(expr_or_seurat[[]]))
        stop("Cluster column '", cluster_labels, "' not found in Seurat meta.data.")
      cl_all <- expr_or_seurat[[cluster_labels]][,1,drop=TRUE]
      names(cl_all) <- cells_all
    } else {
      cl_all <- .align_to_cells(cluster_labels, cells_all, "cluster_labels")
    }
  } else {
    cl_all <- NULL
  }

  # ---- Branch subsetting ----
  subset_for_branch <- function(include=NULL, exclude=NULL, min_cells=10, drop_unused_levels=TRUE) {
    mask <- rep(TRUE, length(cells_all))
    if (!is.null(cl_all)) {
      if (!is.null(include)) mask <- cl_all %in% include
      if (!is.null(exclude)) mask <- mask & !(cl_all %in% exclude)
    }
    idx <- which(mask)
    if (length(idx) < min_cells)
      stop("Branch subsetting left ", length(idx), " cells (< ", min_cells, ").")

    list(
      expr  = expr[, idx, drop=FALSE],
      pt    = pt_all[idx],
      cl    = if (is.null(cl_all)) NULL else as.character(if (drop_unused_levels) factor(cl_all[idx]) else cl_all[idx]),
      cells = cells_all[idx]
    )
  }

  # ---- Per-branch compute ----
  run_branch <- function(b) {
    if (isTRUE(verbose)) message(">>> Branch: ", b)
    filt    <- branch_filters[[b]] %||% list()
    include <- filt$include %||% NULL
    exclude <- filt$exclude %||% NULL
    minc    <- filt$min_cells %||% 10
    drop_lv <- filt$drop_unused_levels %||% TRUE

    sub <- subset_for_branch(include, exclude, minc, drop_lv)

    naive_markers    <- naive_markers_list[[b]]
    terminal_markers <- terminal_markers_list[[b]]
    # ---- D ----
    D_res <- tryCatch({
      metrics_D(sub$expr, naive_markers, terminal_markers, sub$pt)
    }, error = function(e) list(D_naive=NA_real_, D_term=NA_real_, error=e$message))

    # ---- O ----
    O_res <- tryCatch({
      if (.has_formal(metrics_O)) {
        metrics_O(sub$expr, naive_markers, terminal_markers, sub$pt)
      } else {
        metrics_O(sub$expr, naive_markers, terminal_markers, sub$pt)
      }
    }, error = function(e) list(O=NA_real_, error=e$message))

    # ---- E ----
    # robust colMeans: handle empty overlaps to avoid numeric(0)
    n_idx <- intersect(naive_markers, rownames(sub$expr))
    t_idx <- intersect(terminal_markers, rownames(sub$expr))
    naive_scores <- if (length(n_idx)) colMeans(sub$expr[n_idx, , drop=FALSE]) else rep(NA_real_, ncol(sub$expr))
    term_scores  <- if (length(t_idx)) colMeans(sub$expr[t_idx, , drop=FALSE]) else rep(NA_real_, ncol(sub$expr))

    E_res <- tryCatch({
      metrics_E(
        sub$pt,
        cluster_labels = sub$cl,
        naive_marker_scores = naive_scores,
        terminal_marker_scores = term_scores,
        naive_clusters = NULL, terminal_clusters = NULL,
        method = E_method, plot = isTRUE(plot_E)
      )
    }, error = function(e) list(E_naive=NA_real_, E_term=NA_real_, E_comp=NA_real_, error=e$message))



    # ---- Aggregate per-branch score ----
    comp_vec <- c(D_res$D_naive %||% NA_real_,
                  D_res$D_term  %||% NA_real_,
                  O_res$O       %||% NA_real_,
                  E_res$E_comp  %||% NA_real_)
    score <- if (all(is.na(comp_vec))) NA_real_ else mean(comp_vec, na.rm = TRUE)

    # Branch payload
    list(
      branch = b,
      n_cells = ncol(sub$expr),
      D = D_res, O = O_res, E = E_res,
      DOE_score = score,
      cells = sub$cells
    )
  }

  branch_results <- lapply(branch_names, run_branch)
  names(branch_results) <- branch_names

  # ---- Overall aggregate (weighted by branch cell count) ----
  w <- vapply(branch_results, function(x) x$n_cells, numeric(1))
  s <- vapply(branch_results, function(x) x$DOE_score, numeric(1))
  agg <- if (all(is.na(s))) NA_real_ else stats::weighted.mean(s, w, na.rm = TRUE)

  # ---- Return structure expected by the multi-trajectory wrapper ----
  list(
    branches = branch_results,            # each [[b]] has $DOE_score
    aggregate_DOE = agg,                 # overall
    weights = w,                          # branch cell counts
    n_cells_total = sum(w),
    params = list(E_method=E_method)
  )

}


# ---------- Multi-trajectory (branched) ----------
#' Compute DOE across multiple branched trajectories
#'
#' Run `compute_single_DOE_branched()` over a set of pseudotime trajectories
#' that share an expression object. The function collates branch-level results
#' and generates summary tables for downstream plotting and ranking.
#'
#' @inheritParams compute_single_DOE_branched
#' @param expr_or_seurat Expression matrix or Seurat object shared across
#'   trajectories.
#' @param pseudotime_list Named list of pseudotime vectors (one per trajectory).
#' @param parallel Logical; evaluate trajectories in parallel when possible.
#' @param n_cores Number of workers to use for parallel evaluation.
#' @param tol Numerical tolerance forwarded to `compute_single_DOE_branched()`.
#'
#' @return A list with class `"multi_doe_branched"` containing per-trajectory
#'   results, overall/branch comparisons, and metadata about the run.
#'
#' @export
compute_multi_DOE_branched <- function(expr_or_seurat,
                                        pseudotime_list,
                                        naive_markers_list,
                                        terminal_markers_list,
                                        cluster_labels = NULL,
                                        branch_filters = NULL,
                                        E_method = c("gmm","clusters","combined"),
                                        assay = NULL,
                                        slot  = "data",
                                        plot_E = FALSE,
                                        verbose = TRUE,
                                        parallel = FALSE,
                                        n_cores = NULL,
                                        tol = 1e-8) {

  E_method <- match.arg(E_method)

  if (!is.list(pseudotime_list)) stop("pseudotime_list must be a named list.")
  if (is.null(names(pseudotime_list)))
    names(pseudotime_list) <- paste0("trajectory_", seq_along(pseudotime_list))

  run_one <- function(name, pt) {
    if (isTRUE(verbose)) message("\n=== Trajectory: ", name, " ===")
    out <- compute_single_DOE_branched(
      expr_or_seurat = expr_or_seurat,
      pseudotime = pt,
      naive_markers_list = naive_markers_list,
      terminal_markers_list = terminal_markers_list,
      cluster_labels = cluster_labels,
      branch_filters = branch_filters,
      E_method = E_method,
      assay = assay,
      slot = slot,
      plot_E = plot_E,
      verbose = verbose,
      tol = tol
    )
    out$trajectory <- name
    out
  }

  # --- Run (simple parallel optional) ---
  if (isTRUE(parallel) && length(pseudotime_list) > 1 &&
      requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res_list <- parallel::parLapply(
      cl, seq_along(pseudotime_list),
      function(i) {
        name <- names(pseudotime_list)[i]
        pt   <- pseudotime_list[[i]]
        run_one(name, pt)  # NOTE: compute_single_DOE_branched must be visible on workers
      }
    )
    names(res_list) <- names(pseudotime_list)
  } else {
    res_list <- Map(run_one, names(pseudotime_list), pseudotime_list)
    names(res_list) <- names(pseudotime_list)
  }

  # --- Overall comparison ---
  comparison_overall <- data.frame(
    trajectory = names(res_list),
    aggregate_DOE = vapply(res_list, function(x) {
      if (!is.null(x$aggregate_DOE)) as.numeric(x$aggregate_DOE) else NA_real_
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  comparison_overall <- comparison_overall[order(comparison_overall$aggregate_DOE,
                                                 decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  best_overall <- if (all(is.na(comparison_overall$aggregate_DOE))) NA_character_
  else comparison_overall$trajectory[1]

  # --- Per-branch ranking via branches[[br]]$DOE_score ---
  all_branches <- sort(unique(unlist(lapply(res_list, function(x) names(x$branches)))))

  branch_comparisons <- setNames(vector("list", length(all_branches)), all_branches)
  best_per_branch <- vector("list", length(all_branches)); names(best_per_branch) <- all_branches
  branch_data <- setNames(vector("list", length(all_branches)), all_branches)

  for (br in all_branches) {
    scores <- vapply(res_list, function(x) {
      if (!is.null(x$branches) && !is.null(x$branches[[br]]) &&
          !is.null(x$branches[[br]]$DOE_score)) {
        as.numeric(x$branches[[br]]$DOE_score)
      } else {
        NA_real_
      }
    }, numeric(1))

    df <- data.frame(
      trajectory = names(res_list),
      DOE_score = scores,
      stringsAsFactors = FALSE
    )
    df <- df[order(df$DOE_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    branch_comparisons[[br]] <- df

    bt <- if (all(is.na(df$DOE_score))) NA_character_ else df$trajectory[1]
    bs <- if (is.na(bt)) NA_real_ else df$DOE_score[1]
    best_per_branch[[br]] <- data.frame(
      branch = br, best_trajectory = bt, best_score = bs,
      stringsAsFactors = FALSE
    )

    branch_data[[br]] <- lapply(res_list, function(x) {
      if (!is.null(x$branches)) x$branches[[br]] else NULL
    })
    names(branch_data[[br]]) <- names(res_list)
  }

  best_per_branch <- do.call(rbind, best_per_branch)

  comparison_by_branch <- do.call(
    rbind,
    lapply(all_branches, function(br) {
      df <- branch_comparisons[[br]]
      if (!nrow(df)) return(NULL)
      data.frame(branch = br, df, row.names = NULL, stringsAsFactors = FALSE)
    })
  )

  res <- list(
    results = res_list,
    comparison_overall = comparison_overall,
    best_overall = best_overall,
    comparison_by_branch = comparison_by_branch,  # columns: branch, trajectory, DOE_score
    best_per_branch = best_per_branch,            # columns: branch, best_trajectory, best_score
    branch_comparisons = branch_comparisons,      # per-branch sorted tables (DOE_score)
    branch_data = branch_data,                    # branch -> trajectory -> raw branch object
    args = list(E_method=E_method,
                parallel=parallel, n_cores=n_cores)
  )
  class(res) <- "multi_doe_branched"

  if (isTRUE(verbose)) {
    cat("\n==========================================================\n")
    cat("Multi-trajectory branched DOE analysis complete!\n")
    if (!is.na(res$best_overall)) {
      bo <- res$comparison_overall$aggregate_DOE[
        res$comparison_overall$trajectory == res$best_overall
      ]
      cat(sprintf("Best overall trajectory: %s (DOE: %.3f)\n", res$best_overall, bo))
    } else {
      cat("Best overall trajectory: NA (no valid scores)\n")
    }
    if (!is.null(res$best_per_branch) && nrow(res$best_per_branch)) {
      cat("\nBest trajectory per branch:\n")
      for (i in seq_len(nrow(res$best_per_branch))) {
        br <- res$best_per_branch$branch[i]
        bt <- res$best_per_branch$best_trajectory[i]
        sc <- res$best_per_branch$best_score[i]
        if (is.na(bt)) {
          cat(sprintf("  - %s: NA (no valid scores)\n", br))
        } else {
          cat(sprintf("  - %s: %s (DOE: %.3f)\n", br, bt, sc))
        }
      }
    }
  }

  return(res)
}

# ---------- Plot ----------
#' Plot DOE metrics for a single branched trajectory
#'
#' @description
#' Visualize results from [compute_single_DOE_branched()] across branches
#' within one trajectory. Supports bar, heatmap, and radar visualizations.
#'
#' @param single_doe_branched A result list from [compute_single_DOE_branched()].
#' @param metrics Character vector of metric names to include.
#'   Defaults to c("D_naive","D_term","O","E_naive","E_term","DOE_score").
#' @param type Plot type: `"bar"`, `"heatmap"`, or `"radar"`.
#'
#' @return For `"bar"` and `"heatmap"`, a `ggplot` object.
#'   For `"radar"`, draws a base R radar plot and returns `NULL` invisibly.
#' @method plot single_doe_branched
#' @export
plot.single_doe_branched <- function(single_doe_branched,
                                     metrics = c("D_naive","D_term","O",
                                                 "E_naive","E_term","DOE_score"),
                                     type = c("bar","heatmap","radar")) {

  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required. Please install it.")

  # Optional dependencies
  if (type %in% c("bar","heatmap") && !requireNamespace("reshape2", quietly = TRUE))
    stop("reshape2 required (install.packages('reshape2'))")
  if (type == "radar" && !requireNamespace("fmsb", quietly = TRUE))
    stop("fmsb required (install.packages('fmsb'))")
  if (type == "radar" && !requireNamespace("scales", quietly = TRUE))
    stop("scales required (install.packages('scales'))")

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  num1 <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x)) x[1] else NA_real_
  }

  branches <- names(single_doe_branched$branches)
  if (is.null(branches) || !length(branches))
    stop("No branches found in object.")

  rows <- lapply(branches, function(br) {
    br_obj <- single_doe_branched$branches[[br]]
    vals <- setNames(numeric(length(metrics)), metrics)
    vals["D_naive"] <- num1(br_obj$D$D_naive)
    vals["D_term"]  <- num1(br_obj$D$D_term)
    vals["O"]       <- num1(br_obj$O$O)
    vals["E_naive"] <- num1(br_obj$E$E_naive)
    vals["E_term"]  <- num1(br_obj$E$E_term)
    vals["DOE_score"] <- num1(br_obj$DOE_score)
    data.frame(branch = br, as.list(vals), check.names = FALSE, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  metrics <- intersect(metrics, names(df))

  if (type == "bar") {
    long <- reshape2::melt(df, id.vars = "branch", variable.name = "metric", value.name = "score")
    p <- ggplot2::ggplot(long, ggplot2::aes(x = branch, y = score, fill = metric)) +
      ggplot2::geom_col(position = "dodge", na.rm = TRUE) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "DOE Metrics per Branch", x = "Branch", y = "Score") +
      ggplot2::coord_cartesian(ylim = c(0, 1))
    return(p)
  }

  if (type == "heatmap") {
    long <- reshape2::melt(df, id.vars = "branch", variable.name = "metric", value.name = "score")
    p <- ggplot2::ggplot(long, ggplot2::aes(x = metric, y = branch, fill = score)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                                    midpoint = 0.5, name = "Score", na.value = "grey90") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = "DOE Metrics Heatmap (Branches)", x = "Metric", y = "Branch") +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(score), "NA", sprintf("%.2f", score))),
                         color = "black", size = 3)
    return(p)
  }

  if (type == "radar") {
    mat <- df[, metrics, drop = FALSE]
    mat[is.na(mat)] <- 0
    for (m in metrics) mat[[m]] <- pmax(0, pmin(1, mat[[m]]))
    rad <- rbind(rep(1, length(metrics)), rep(0, length(metrics)), mat)
    rownames(rad) <- c("Max","Min", df$branch)
    n <- nrow(df)
    fmsb::radarchart(rad, axistype = 1,
                     pcol = seq_len(n),
                     pfcol = scales::alpha(seq_len(n), 0.25),
                     plwd = 2, plty = 1,
                     cglcol = "grey", cglty = 1, axislabcol = "grey",
                     caxislabels = rep("",5), cglwd = 0.5, vlcex = 0,
                     title = "DOE Metrics Radar (Single Branched Trajectory)")
    graphics::legend("topright", legend = df$branch,
                     col = seq_len(n), lty = 1, lwd = 2, bty = "n", cex = 0.8)
    return(invisible(NULL))
  }

  stop("Unknown type.")
}



#' Plot DOE metrics for multi-trajectory branched analysis
#'
#' This function visualizes the results from
#' [compute_multi_DOE_branched()] across trajectories and branches.
#' It supports bar plots, heatmaps, and radar charts, and can display
#' either overall aggregate DOE scores or per-branch metrics.
#'
#' @param multi_doe_branched A result object returned by
#'   [compute_multi_DOE_branched()], containing overall comparisons,
#'   per-branch comparisons, and branch-level data.
#' @param scope Character scalar, either `"overall"` to plot aggregate
#'   DOE scores across trajectories, or `"branch"` to plot per-branch
#'   metrics. Default is `"branch"`.
#' @param metrics Character vector of metric names to plot. Common options
#'   include `"D_naive"`, `"D_term"`, `"O"`, `"E_naive"`,
#'   `"E_term"`, and `"DOE_score"`. Defaults to all of these.
#' @param type Character scalar, plot type: `"bar"`, `"heatmap"`,
#'   or `"radar"`.
#' @param branch_mode For branch scope only. Controls how branches are
#'   displayed:
#'   * `"facet"` ??? facet each branch into a separate panel,
#'   * `"stack"` ??? combine trajectory and branch labels
#'     (e.g., `"traj1 [B]"`),
#'   * `"separate"` ??? for radar plots only, draw one radar per branch.
#'   Default is `"facet"`.
#' @param branches Optional character vector of branch names to include
#'   (for `scope = "branch"`). Defaults to all branches present.
#'
#' @return For `"bar"` and `"heatmap"`, returns a `ggplot` object.
#' For `"radar"`, draws the plot using base graphics and returns `NULL`
#' invisibly.
#'
#' @details
#' * Overall plots use `comparison_overall` from the branched DOE result,
#'   showing each trajectory???s aggregate DOE score.
#' * Branch plots iterate through `branch_data` and extract requested
#'   metrics for each (trajectory, branch) pair. Metrics can be stored as
#'   numeric values or nested lists/doubles inside each branch object; the
#'   function extracts the first valid numeric value.
#'
#' @examples NULL
#' @method plot multi_doe_branched
#' @export
plot.multi_doe_branched <- function(multi_doe_branched,
                                    scope = c("branch","overall"),
                                    metrics = c("D_naive","D_term","O","E_naive","E_term","DOE_score"),
                                    type = c("bar","heatmap","radar"),
                                    branch_mode = c("facet","stack","separate"),
                                    branches = NULL) {
  scope       <- match.arg(scope)
  type        <- match.arg(type)
  branch_mode <- match.arg(branch_mode)
  
  # -------- Dependencies --------
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (type %in% c("bar","heatmap") && !requireNamespace("reshape2", quietly = TRUE))
    stop("reshape2 required for long data (install.packages('reshape2'))")
  if (type == "radar" && !requireNamespace("fmsb", quietly = TRUE))
    stop("fmsb required for radar (install.packages('fmsb'))")
  if (type == "radar" && !requireNamespace("scales", quietly = TRUE))
    stop("scales required for radar (install.packages('scales'))")
  
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
  num1 <- function(x) {
    if (is.null(x)) return(NA_real_)
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x)) x[1] else NA_real_
  }
  
  # extractor supports common nested fields inside each branch object
  extract_metric <- function(obj, m) {
    if (is.null(obj)) return(NA_real_)
    switch(m,
           "D_naive"   = num1(obj$D$D_naive %||% obj$D_naive %||% obj$D$naive),
           "D_term"    = num1(obj$D$D_term  %||% obj$D_term  %||% obj$D$terminal),
           "O"         = num1(obj$O$O       %||% obj$O),
           "E_naive"   = num1(obj$E$E_naive %||% obj$E_naive %||% obj$E$naive),
           "E_term"    = num1(obj$E$E_term  %||% obj$E_term  %||% obj$E$terminal),
           "DOE_score"= num1(obj$DOE_score %||% obj$DOE %||% obj$score),
           # fallback: try direct name
           num1(obj[[m]])
    )
  }
  
  # -------- Build tidy data --------
  if (scope == "overall") {
    # Use the aggregate DOE by trajectory (mirrors your non-branched plot logic)
    plot_data <- multi_doe_branched$comparison_overall
    if (is.null(plot_data) || !nrow(plot_data)) stop("No overall comparison data found.")
    # Permit plotting a single metric: rename to DOE_score to reuse code paths
    names(plot_data)[names(plot_data) == "aggregate_DOE"] <- "DOE_score"
    metrics <- intersect(metrics, names(plot_data))
    if (!length(metrics)) metrics <- "DOE_score"
    
    if (type == "bar") {
      long <- reshape2::melt(plot_data[, c("trajectory", metrics), drop = FALSE],
                             id.vars = "trajectory", variable.name = "metric", value.name = "score")
      p <- ggplot2::ggplot(long, ggplot2::aes(x = trajectory, y = score, fill = metric)) +
        ggplot2::geom_col(position = "dodge", na.rm = TRUE) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Aggregate DOE Across Trajectories",
                      x = "Trajectory", y = "Score") +
        ggplot2::coord_cartesian(ylim = c(0, 1))
      return(p)
    }
    
    if (type == "heatmap") {
      long <- reshape2::melt(plot_data[, c("trajectory", metrics), drop = FALSE],
                             id.vars = "trajectory", variable.name = "metric", value.name = "score")
      p <- ggplot2::ggplot(long, ggplot2::aes(x = metric, y = trajectory, fill = score)) +
        ggplot2::geom_tile(color = "white", size = 0.5) +
        ggplot2::scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                                      midpoint = 0.5, name = "Score", na.value = "grey90") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title = "Aggregate DOE Heatmap", x = "Metric", y = "Trajectory") +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(score), "NA", sprintf("%.2f", score))),
                           color = "black", size = 3)
      return(p)
    }
    
    if (type == "radar") {
      rd <- plot_data[, metrics, drop = FALSE]
      rd[is.na(rd)] <- 0
      for (m in metrics) rd[[m]] <- pmax(0, pmin(1, rd[[m]]))
      rd <- rbind(rep(1, length(metrics)), rep(0, length(metrics)), rd)
      rownames(rd) <- c("Max","Min", plot_data$trajectory)
      n <- nrow(plot_data)
      fmsb::radarchart(rd, axistype = 1,
                       pcol = seq_len(n),
                       pfcol = scales::alpha(seq_len(n), 0.25),
                       plwd = 2, plty = 1, caxislabels = rep("", 5),
                       cglcol = rgb(0, 0, 0, alpha = 0), cglty = 1,vlcex = 0,
                       cglwd = 0.5, vlcex = 0.8, title = "Aggregate DOE Radar")
      graphics::legend("topright", legend = plot_data$trajectory,
                       col = seq_len(n), lty = 1, lwd = 2, bty = "n", cex = 0.8)
      return(invisible(NULL))
    }
  } else {
    # scope == "branch": assemble [branch, trajectory, metrics...] from branch_data
    bd <- multi_doe_branched$branch_data
    if (is.null(bd) || !length(bd)) stop("No branch_data found.")
    # choose branches to include
    if (is.null(branches)) branches <- names(bd)
    branches <- intersect(branches, names(bd))
    if (!length(branches)) stop("No matching branches to plot.")
    
    rows <- list()
    for (br in branches) {
      tr_map <- bd[[br]]                 # trajectory -> branch_object
      if (is.null(tr_map) || !length(tr_map)) next
      for (traj in names(tr_map)) {
        obj <- tr_map[[traj]]
        vals <- setNames(numeric(length(metrics)), metrics)
        for (m in metrics) vals[m] <- extract_metric(obj, m)
        rows[[length(rows) + 1]] <- data.frame(
          branch = br,
          trajectory = traj,
          as.list(vals),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
    if (!length(rows)) stop("No branch rows assembled.")
    plot_df <- do.call(rbind, rows)
    # keep only present metrics
    metrics <- intersect(metrics, setdiff(names(plot_df), c("branch","trajectory")))
    if (!length(metrics)) stop("Requested metrics not found in branch results.")
    
    # Label for stack mode
    plot_df$traj_label <- paste0(plot_df$trajectory, " [", plot_df$branch, "]")
    
    if (type == "bar") {
      long <- reshape2::melt(plot_df[, c("branch","trajectory","traj_label", metrics), drop = FALSE],
                             id.vars = c("branch","trajectory","traj_label"),
                             variable.name = "metric", value.name = "score")
      # x label depends on branch_mode
      x_var <- if (branch_mode == "stack") "traj_label" else "trajectory"
      p <- ggplot2::ggplot(long, ggplot2::aes_string(x = x_var, y = "score", fill = "metric")) +
        ggplot2::geom_col(position = "dodge", na.rm = TRUE) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "DOE Metrics by Branch", x = "Trajectory", y = "Score") +
        ggplot2::coord_cartesian(ylim = c(0, 1))
      if (branch_mode == "facet") {
        p <- p + ggplot2::facet_wrap(~ branch, scales = "free_x")
      }
      return(p)
    }
    
    if (type == "heatmap") {
      long <- reshape2::melt(plot_df[, c("branch","trajectory","traj_label", metrics), drop = FALSE],
                             id.vars = c("branch","trajectory","traj_label"),
                             variable.name = "metric", value.name = "score")
      y_var <- if (branch_mode == "stack") "traj_label" else "trajectory"
      p <- ggplot2::ggplot(long, ggplot2::aes_string(x = "metric", y = y_var, fill = "score")) +
        ggplot2::geom_tile(color = "white", size = 0.5) +
        ggplot2::scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                                      midpoint = 0.5, name = "Score", na.value = "grey90") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title = "DOE Metrics Heatmap (Branches)", x = "Metric", y = "Trajectory") +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(score), "NA", sprintf("%.2f", score))),
                           color = "black", size = 3)
      if (branch_mode == "facet") {
        p <- p + ggplot2::facet_wrap(~ branch, scales = "free_y")
      }
      return(p)
    }
    
    if (type == "radar") {
      # separate: one radar per branch; others: a single radar containing all rows
      to_radar <- function(df_rows, title = "DOE Metrics Radar") {
        mat <- as.data.frame(df_rows[, metrics, drop = FALSE])
        mat[is.na(mat)] <- 0
        for (m in metrics) mat[[m]] <- pmax(0, pmin(1, mat[[m]]))
        rad <- rbind(rep(1, length(metrics)), rep(0, length(metrics)), mat)
        rownames(rad) <- c("Max","Min", rownames(df_rows))
        n <- nrow(df_rows)
        fmsb::radarchart(rad, axistype = 1,
                         pcol = seq_len(n),
                         pfcol = scales::alpha(seq_len(n), 0.25),
                         plwd = 2, plty = 1,
                         cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = rep("",5),
                         cglwd = 0.5, vlcex = 0, title = title)
        graphics::legend("topright", legend = rownames(df_rows),
                         col = seq_len(n), lty = 1, lwd = 2, bty = "n", cex = 0.8)
      }
      
      if (branch_mode == "separate") {
        opar <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(opar))
        brs <- unique(plot_df$branch)
        n <- length(brs); nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
        graphics::par(mfrow = c(nr, nc), mar = c(1,1,2,1))
        for (b in brs) {
          dfb <- plot_df[plot_df$branch == b, c("trajectory", metrics), drop = FALSE]
          rownames(dfb) <- dfb$trajectory
          dfb$trajectory <- NULL
          to_radar(dfb, title = paste("Branch:", b))
        }
        return(invisible(NULL))
      } else {
        # one radar: rows are either trajectory or trajectory [branch]
        df <- plot_df[, c(if (branch_mode == "stack") "traj_label" else "trajectory", metrics), drop = FALSE]
        rownames(df) <- df[[if (branch_mode == "stack") "traj_label" else "trajectory"]]
        df[[if (branch_mode == "stack") "traj_label" else "trajectory"]] <- NULL
        to_radar(df, title = "DOE Metrics Radar (Branches)")
        return(invisible(NULL))
      }
    }
  }
  
  stop("Unknown configuration.")
}