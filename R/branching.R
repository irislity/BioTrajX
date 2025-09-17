# =========================
# Branched DOPE (rewritten)
# =========================

# ---------- Utilities ----------
.stop_if_missing <- function(lst, required_names, label) {
  miss <- setdiff(required_names, names(lst))
  if (length(miss)) stop("Missing names in ", label, ": ", paste(miss, collapse = ", "))
}

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
      branch       = branch,
      D_naive      = res$D$D_naive %||% NA_real_,
      D_term       = res$D$D_term  %||% NA_real_,
      O            = res$O$O       %||% NA_real_,
      P            = res$P$P       %||% NA_real_,
      E_naive      = res$E$E_naive %||% NA_real_,
      E_term       = res$E$E_term  %||% NA_real_,
      E_comp       = res$E$E_comp  %||% NA_real_,
      DOPE_score   = res$DOPE_score %||% NA_real_,
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
#' Compute DOPE for one pseudotime with multiple branches
compute_DOPE_single_branched <- function(expr_or_seurat,
                                         pseudotime,
                                         naive_markers_list,
                                         terminal_markers_list,
                                         pathways_list = NULL,
                                         gene_sets_list = NULL,
                                         cluster_labels = NULL,
                                         branch_filters = NULL,
                                         E_method = c("gmm","clusters","combined"),
                                         id_type = c("SYMBOL","ENSEMBL","ENTREZID"),
                                         species = "Homo sapiens",
                                         assay = NULL,
                                         slot = "data",
                                         min_remaining = 10,
                                         min_fraction  = 0.20,
                                         min_genes_per_module = 3,
                                         tol = 1e-8,
                                         plot_E = FALSE,
                                         verbose = TRUE) {
  E_method <- match.arg(E_method)
  id_type  <- match.arg(id_type)

  # Branch sanity
  branch_names <- union(names(naive_markers_list), names(terminal_markers_list))
  if (is.null(branch_names) || any(branch_names == "")) stop("Branch lists must be named.")
  .stop_if_missing(naive_markers_list,    branch_names, "naive_markers_list")
  .stop_if_missing(terminal_markers_list, branch_names, "terminal_markers_list")
  if (!is.null(pathways_list))  .stop_if_missing(pathways_list,  branch_names, "pathways_list")
  if (!is.null(gene_sets_list)) .stop_if_missing(gene_sets_list, branch_names, "gene_sets_list")

  if (!is.null(pathways_list) && !is.null(gene_sets_list)) {
    warning("Both pathways_list and gene_sets_list provided; gene_sets_list will take precedence for P metric.")
  }

  # Expression extraction
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
  if (is.null(cells_all) || length(cells_all) == 0) stop("No column (cell) names found in expression matrix.")

  # Align helper
  .align_to_cells <- function(x, target_cells, what="vector") {
    if (!is.null(names(x))) {
      miss <- setdiff(target_cells, names(x))
      if (length(miss)) stop("Named ", what, " missing cells: ", paste(head(miss, 5), collapse = ", "), if (length(miss) > 5) " ...")
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

  # --- Per-branch compute ---
  run_branch <- function(b) {
    if (isTRUE(verbose)) message(">>> Branch: ", b)
    filt   <- branch_filters[[b]] %||% list()
    include <- filt$include %||% NULL
    exclude <- filt$exclude %||% NULL
    minc    <- filt$min_cells %||% 10
    drop_lv <- filt$drop_unused_levels %||% TRUE

    sub <- subset_for_branch(include, exclude, minc, drop_lv)

    naive_markers    <- naive_markers_list[[b]]
    terminal_markers <- terminal_markers_list[[b]]
    gene_sets        <- if (!is.null(gene_sets_list))  gene_sets_list[[b]]  else NULL
    pathways         <- if (!is.null(pathways_list))   pathways_list[[b]]   else NULL

    # D
    D_res <- tryCatch({
      metrics_D(sub$expr, naive_markers, terminal_markers, sub$pt)
    }, error = function(e) list(D_naive=NA_real_, D_term=NA_real_, error=e$message))

    # O (tol passed only if supported)
    O_res <- tryCatch({
      if (.has_formal(metrics_O, "tol")) {
        metrics_O(sub$expr, naive_markers, terminal_markers, sub$pt, tol = tol)
      } else {
        metrics_O(sub$expr, naive_markers, terminal_markers, sub$pt)
      }
    }, error = function(e) list(O=NA_real_, error=e$message))

    # P (gene_sets wins if both given)
    P_res <- tryCatch({
      metrics_P(
        sub$expr, sub$pt,
        pathways = if (is.null(gene_sets)) pathways else NULL,
        gene_sets = gene_sets,
        naive_markers = naive_markers,
        terminal_markers = terminal_markers,
        id_type = id_type,
        species = species,
        assay = NULL, slot = "data",
        min_remaining = min_remaining,
        min_fraction = min_fraction,
        min_genes_per_module = min_genes_per_module,
        verbose = verbose
      )
    }, error = function(e) list(P=NA_real_, error=e$message))

    # E
    naive_scores <- colMeans(sub$expr[intersect(naive_markers, rownames(sub$expr)), , drop=FALSE])
    term_scores  <- colMeans(sub$expr[intersect(terminal_markers, rownames(sub$expr)), , drop=FALSE])
    E_res <- tryCatch({
      metrics_E(
        sub$pt,
        cluster_labels = sub$cl,
        naive_marker_scores = naive_scores,
        terminal_marker_scores = term_scores,
        naive_clusters = NULL, terminal_clusters = NULL,
        method = E_method, plot = FALSE
      )
    }, error = function(e) list(E_naive=NA_real_, E_term=NA_real_, E_comp=NA_real_, error=e$message))

    comp_vec <- c(D_res$D_naive %||% NA_real_,
                  D_res$D_term  %||% NA_real_,
                  O_res$O       %||% NA_real_,
                  P_res$P       %||% NA_real_,
                  E_res$E_comp  %||% NA_real_)
    score <- if (all(is.na(comp_vec))) NA_real_ else mean(comp_vec, na.rm = TRUE)

    structure(list(
      class = "dope_results",
      branch = b,
      n_cells = ncol(sub$expr),
      D = D_res, O = O_res, P = P_res, E = E_res,
      DOPE_score = score
    ), class = "dope_results")
  }

  branch_results <- lapply(branch_names, run_branch)
  names(branch_results) <- branch_names

  # Aggregate (weighted by cells)
  w <- vapply(branch_results, function(x) x$n_cells, numeric(1))
  s <- vapply(branch_results, function(x) x$DOPE_score, numeric(1))
  agg <- if (all(is.na(s))) NA_real_ else stats::weighted.mean(s, w, na.rm = TRUE)

  list(
    branches = branch_results,
    aggregate_DOPE = agg,
    weights = w,
    params = list(E_method=E_method, id_type=id_type, species=species,
                  min_remaining=min_remaining, min_fraction=min_fraction,
                  min_genes_per_module=min_genes_per_module, tol=tol)
  )
}

# ---------- Multi-trajectory (branched) ----------
#' Compute DOPE across multiple pseudotime trajectories with branches
compute_multi_DOPE_branched <- function(expr_or_seurat,
                                        pseudotime_list,
                                        naive_markers_list,
                                        terminal_markers_list,
                                        pathways_list = NULL,
                                        gene_sets_list = NULL,
                                        cluster_labels = NULL,
                                        branch_filters = NULL,
                                        E_method = c("gmm","clusters","combined"),
                                        id_type = c("SYMBOL","ENSEMBL","ENTREZID"),
                                        species = "Homo sapiens",
                                        assay = NULL,
                                        slot  = "data",
                                        min_remaining = 10,
                                        min_fraction  = 0.20,
                                        min_genes_per_module = 3,
                                        tol = 1e-8,
                                        plot_E = FALSE,
                                        verbose = TRUE,
                                        parallel = FALSE,
                                        n_cores = NULL) {
  E_method <- match.arg(E_method)
  id_type  <- match.arg(id_type)

  if (!is.list(pseudotime_list)) stop("pseudotime_list must be a named list.")
  if (is.null(names(pseudotime_list)))
    names(pseudotime_list) <- paste0("trajectory_", seq_along(pseudotime_list))

  run_one <- function(name, pt) {
    if (isTRUE(verbose)) message("\n=== Trajectory: ", name, " ===")
    out <- compute_DOPE_single_branched(
      expr_or_seurat = expr_or_seurat,
      pseudotime = pt,
      naive_markers_list = naive_markers_list,
      terminal_markers_list = terminal_markers_list,
      pathways_list = pathways_list,
      gene_sets_list = gene_sets_list,
      cluster_labels = cluster_labels,
      branch_filters = branch_filters,
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
    out$trajectory <- name
    out
  }

  if (parallel && length(pseudotime_list) > 1 && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export up-to-date objects and functions
    parallel::clusterExport(
      cl,
      varlist = c(
        # current environment helpers
        ".stop_if_missing", "%||%", ".has_formal", ".flatten_traj_result",
        # core functions and external metric symbols
        "compute_DOPE_single_branched", "metrics_D", "metrics_O", "metrics_P", "metrics_E"
      ),
      envir = environment()
    )
    # Load packages on workers if needed (uncomment as appropriate):
    # parallel::clusterEvalQ(cl, { library(Seurat); library(Matrix) })

    res_list <- parallel::parLapply(cl, seq_along(pseudotime_list), function(i) {
      name <- names(pseudotime_list)[i]; pt <- pseudotime_list[[i]]
      run_one(name, pt)
    })
    names(res_list) <- names(pseudotime_list)
  } else {
    res_list <- Map(run_one, names(pseudotime_list), pseudotime_list)
    names(res_list) <- names(pseudotime_list)
  }

  comp <- data.frame(
    trajectory = names(res_list),
    aggregate_DOPE = vapply(res_list, function(x) x$aggregate_DOPE, numeric(1)),
    stringsAsFactors = FALSE
  )
  comp <- comp[order(comp$aggregate_DOPE, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  best <- if (all(is.na(comp$aggregate_DOPE))) NA_character_ else comp$trajectory[1]

  list(
    results = res_list,
    comparison = comp,
    best_trajectory = best,
    args = list(E_method=E_method, id_type=id_type, species=species, parallel=parallel, n_cores=n_cores)
  )
}

# ---------- Summary / Print / Plot ----------
create_comparison_summary_branched <- function(trajectory_results) {
  if (length(trajectory_results) == 0) stop("trajectory_results cannot be empty")

  rows <- vector("list", length(trajectory_results$results))
  for (i in seq_along(trajectory_results$results)) {
    traj_name <- names(trajectory_results$results)[i] %||% paste0("Trajectory_", i)
    res_i     <- trajectory_results$results[[i]]
    rows[[i]] <- .flatten_traj_result(res_i, traj_name)
  }
  summary_data <- do.call(rbind, rows)

  metric_cols <- c("D_naive","D_term","O","P","E_naive","E_term","E_comp","DOPE_score")
  metric_cols <- intersect(metric_cols, names(summary_data))
  summary_data <- .add_metric_ranks(summary_data, metric_cols)

  ord <- order(summary_data$DOPE_score, decreasing = TRUE, na.last = TRUE,
               summary_data$trajectory, summary_data$branch)
  summary_data <- summary_data[ord, , drop = FALSE]
  rownames(summary_data) <- NULL
  summary_data
}

print.multi_dope_results <- function(x, ...) {
  cat("Multi-Trajectory DOPE Analysis\n")
  cat("==============================\n\n")
  cat("Analysis info:\n")
  cat(sprintf("  Trajectories analyzed: %d\n", x$n_trajectories %||% length(x$results)))
  if (!is.null(x$method_info$E_method)) cat(sprintf("  E method: %s\n", x$method_info$E_method))
  if (isTRUE(x$method_info$parallel)) cat(sprintf("  Parallel processing: %d cores\n", x$method_info$n_cores))
  cat(sprintf("  Trajectory mode: %s\n\n", x$method_info$trajectory %||% "branched"))

  cat("Top Trajectories/Branches:\n")
  cat("--------------------------\n")
  top_n <- min(5, nrow(x$comparison_summary))
  top_rows <- x$comparison_summary[seq_len(top_n), ]

  for (i in seq_len(nrow(top_rows))) {
    row <- top_rows[i, ]
    rank_suffix <- c("st","nd","rd","th")[pmin(i,4)]
    name <- if (!is.na(row$branch)) sprintf("%s ⟂ %s", row$trajectory, row$branch) else row$trajectory
    dope_score <- if (is.na(row$DOPE_score)) "NA" else sprintf("%.3f", row$DOPE_score)
    cat(sprintf("%d%s: %s (DOPE: %s)\n", i, rank_suffix, name, dope_score))
  }

  if (!is.na(x$best_trajectory)) cat(sprintf("\nBest: %s\n", x$best_trajectory))
  cat("\nUse summary() for detailed comparison table\n")
}

summary.multi_dope_results <- function(object, ...) {
  cat("\nDetailed Comparison Table:\n")
  cat("=========================\n")
  display_cols <- c("trajectory","branch","DOPE_score","DOPE_score_rank","D_term","O","P","E_comp")
  display_cols <- display_cols[display_cols %in% names(object$comparison_summary)]
  display_data <- object$comparison_summary[, display_cols, drop = FALSE]
  numeric_cols <- vapply(display_data, is.numeric, logical(1))
  display_data[numeric_cols] <- lapply(display_data[numeric_cols], function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x)))
  print(display_data, row.names = FALSE)

  if (!is.na(object$best_trajectory)) {
    cat(sprintf("\nDetailed results for best (%s):\n", object$best_trajectory))
    cat("================================\n")
    bt <- object$best_trajectory
    best_result <- NULL
    if (!is.null(object$results[[bt]])) {
      best_result <- object$results[[bt]]
    } else {
      parts <- strsplit(bt, "::", fixed = TRUE)[[1]]
      if (length(parts) == 2 && !is.null(object$results[[parts[1]]])) {
        best_result <- object$results[[parts[1]]]$branches[[parts[2]]]
      }
    }
    if (!is.null(best_result) && inherits(best_result, "dope_results")) print(best_result)
  }
  invisible(object)
}

plot.multi_dope_results <- function(multi_dope_results,
                                    metrics = c("D_naive","D_term","O","P","E_naive","E_term","DOPE_score"),
                                    type = "bar",
                                    branch_mode = c("facet","stack","separate")) {
  branch_mode <- match.arg(branch_mode)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package required for plotting")
  cs <- multi_dope_results$comparison_summary

  metrics <- intersect(metrics, names(cs))
  if (length(metrics) == 0L) stop("No valid metrics to plot.")

  has_branch <- "branch" %in% names(cs) && any(!is.na(cs$branch))
  label_var <- "trajectory"
  if (has_branch && branch_mode == "stack") {
    cs$traj_label <- ifelse(is.na(cs$branch), cs$trajectory, paste0(cs$trajectory, " [", cs$branch, "]"))
    label_var <- "traj_label"
  }

  if (type == "bar") {
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("reshape2 package required for bar plots")
    plot_df <- cs[, c(label_var, "trajectory", "branch", metrics), drop = FALSE]
    plot_long <- reshape2::melt(plot_df, id.vars = c(label_var, "trajectory", "branch"),
                                variable.name = "metric", value.name = "score")

    p <- ggplot2::ggplot(plot_long, ggplot2::aes_string(x = label_var, y = "score", fill = "metric")) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "DOPE Metrics Comparison", x = "Trajectory", y = "Score") +
      ggplot2::ylim(0, 1)

    if (has_branch && branch_mode == "facet") p <- p + ggplot2::facet_wrap(~ branch, scales = "free_x")
    return(p)

  } else if (type == "heatmap") {
    if (!requireNamespace("reshape2", quietly = TRUE)) stop("reshape2 package required for heatmap plots")
    plot_df <- cs[, c(label_var, "trajectory", "branch", metrics), drop = FALSE]
    plot_long <- reshape2::melt(plot_df, id.vars = c(label_var, "trajectory", "branch"),
                                variable.name = "metric", value.name = "score")

    p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = metric, y = !!as.name(label_var), fill = score)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                                    midpoint = 0.5, name = "Score", na.value = "grey90") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = "DOPE Metrics Heatmap", x = "Metric", y = "Trajectory") +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(score), "NA", sprintf("%.2f", score))),
                         color = "black", size = 3)

    if (has_branch && branch_mode == "facet") p <- p + ggplot2::facet_wrap(~ branch, scales = "free_y")
    return(p)

  } else if (type == "radar") {
    if (!requireNamespace("fmsb", quietly = TRUE))
      stop("fmsb package required for radar plots. Install with: install.packages('fmsb')")
    prep_radar <- function(df) {
      mat <- as.data.frame(df[, metrics, drop = FALSE])
      mat[is.na(mat)] <- 0
      for (m in metrics) mat[[m]] <- pmax(0, pmin(1, mat[[m]]))
      rbind(rep(1, length(metrics)), rep(0, length(metrics)), mat)
    }

    if (has_branch && branch_mode == "separate") {
      opar <- par(no.readonly = TRUE); on.exit(par(opar))
      brs <- unique(cs$branch[!is.na(cs$branch)])
      n <- length(brs); nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
      par(mfrow = c(nr, nc), mar = c(1,1,2,1))
      for (b in brs) {
        dfb <- cs[cs$branch == b, c(label_var, metrics), drop = FALSE]
        rownames(dfb) <- dfb[[label_var]]
        dfb[[label_var]] <- NULL
        rad <- prep_radar(dfb)
        rownames(rad) <- c("Max","Min", rownames(dfb))
        fmsb::radarchart(rad, axistype = 1, pcol = 1:nrow(dfb), pfcol = scales::alpha(1:nrow(dfb), .2),
                         plwd = 2, plty = 1,
                         cglcol = "grey", cglty = 1, axislabcol = "grey",
                         cglwd = 0.5, vlcex = 0.8, title = paste("Branch:", b))
        legend("topright", legend = rownames(dfb), col = 1:nrow(dfb), lty = 1, lwd = 2, bty = "n", cex = 0.8)
      }
      return(invisible(NULL))
    } else {
      df <- cs[, c(label_var, metrics), drop = FALSE]
      rownames(df) <- df[[label_var]]
      df[[label_var]] <- NULL
      rad <- prep_radar(df)
      rownames(rad) <- c("Max","Min", rownames(df))
      fmsb::radarchart(rad, axistype = 1, pcol = 1:nrow(df), pfcol = scales::alpha(1:nrow(df), .2),
                       plwd = 2, plty = 1,
                       cglcol = "grey", cglty = 1, axislabcol = "grey",
                       cglwd = 0.5, vlcex = 0.8, title = "DOPE Metrics Radar")
      legend("topright", legend = rownames(df), col = 1:nrow(df), lty = 1, lwd = 2, bty = "n", cex = 0.8)
      return(invisible(NULL))
    }
  } else {
    stop("Plot type must be one of: 'bar', 'radar', or 'heatmap'")
  }
}

# ---------- Helper to assemble a nice S3-like container ----------
as.multi_dope_results <- function(multi_out) {
  comp_sum <- create_comparison_summary_branched(multi_out)
  structure(list(
    results = multi_out$results,
    comparison = multi_out$comparison,
    comparison_summary = comp_sum,
    best_trajectory = multi_out$best_trajectory %||% NA_character_,
    n_trajectories = length(multi_out$results),
    method_info = list(
      E_method = multi_out$args$E_method %||% NA_character_,
      parallel = multi_out$args$parallel %||% FALSE,
      n_cores  = multi_out$args$n_cores %||% NA_integer_,
      trajectory = "branched"
    )
  ), class = "multi_dope_results")
}
