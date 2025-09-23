plot.multi_dope_results_branched <- function(multi_dope_branched,
                                             scope = c("branch","overall"),
                                             metrics = c("D_naive","D_term","O","P","E_naive","E_term","DOPE_score"),
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
           "P"         = num1(obj$P$P       %||% obj$P),
           "E_naive"   = num1(obj$E$E_naive %||% obj$E_naive %||% obj$E$naive),
           "E_term"    = num1(obj$E$E_term  %||% obj$E_term  %||% obj$E$terminal),
           "DOPE_score"= num1(obj$DOPE_score %||% obj$DOPE %||% obj$score),
           # fallback: try direct name
           num1(obj[[m]])
    )
  }

  # -------- Build tidy data --------
  if (scope == "overall") {
    # Use the aggregate DOPE by trajectory (mirrors your non-branched plot logic)
    plot_data <- multi_dope_branched$comparison_overall
    if (is.null(plot_data) || !nrow(plot_data)) stop("No overall comparison data found.")
    # Permit plotting a single metric: rename to DOPE_score to reuse code paths
    names(plot_data)[names(plot_data) == "aggregate_DOPE"] <- "DOPE_score"
    metrics <- intersect(metrics, names(plot_data))
    if (!length(metrics)) metrics <- "DOPE_score"

    if (type == "bar") {
      long <- reshape2::melt(plot_data[, c("trajectory", metrics), drop = FALSE],
                             id.vars = "trajectory", variable.name = "metric", value.name = "score")
      p <- ggplot2::ggplot(long, ggplot2::aes(x = trajectory, y = score, fill = metric)) +
        ggplot2::geom_col(position = "dodge", na.rm = TRUE) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Aggregate DOPE Across Trajectories",
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
        ggplot2::labs(title = "Aggregate DOPE Heatmap", x = "Metric", y = "Trajectory") +
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
                       plwd = 2, plty = 1,
                       cglcol = "grey", cglty = 1, axislabcol = "grey",
                       cglwd = 0.5, vlcex = 0.8, title = "Aggregate DOPE Radar")
      graphics::legend("topright", legend = plot_data$trajectory,
                       col = seq_len(n), lty = 1, lwd = 2, bty = "n", cex = 0.8)
      return(invisible(NULL))
    }
  } else {
    # scope == "branch": assemble [branch, trajectory, metrics...] from branch_data
    bd <- multi_dope_branched$branch_data
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
        ggplot2::labs(title = "DOPE Metrics by Branch", x = "Trajectory", y = "Score") +
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
        ggplot2::labs(title = "DOPE Metrics Heatmap (Branches)", x = "Metric", y = "Trajectory") +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(score), "NA", sprintf("%.2f", score))),
                           color = "black", size = 3)
      if (branch_mode == "facet") {
        p <- p + ggplot2::facet_wrap(~ branch, scales = "free_y")
      }
      return(p)
    }

    if (type == "radar") {
      # separate: one radar per branch; others: a single radar containing all rows
      to_radar <- function(df_rows, title = "DOPE Metrics Radar") {
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
                         cglcol = "grey", cglty = 1, axislabcol = "grey",
                         cglwd = 0.5, vlcex = 0.8, title = title)
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
        to_radar(df, title = "DOPE Metrics Radar (Branches)")
        return(invisible(NULL))
      }
    }
  }

  stop("Unknown configuration.")
}
