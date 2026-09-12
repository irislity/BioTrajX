# =============================================================================
# test-doe_wrapper.R  —  Unit tests for compute_single_doe_linear() and
#                         compute_multi_doe_linear()
# =============================================================================
library(testthat)

# ── Helper: cluster labels aligned to fixture ─────────────────────────────────
make_cluster_labels <- function(f) {
  n <- ncol(f$expr)
  cl <- c(rep("early_cl", 15), rep("mid_cl", n - 30), rep("term_cl", 15))
  setNames(cl, colnames(f$expr))
}

# =============================================================================
# compute_single_doe_linear
# =============================================================================

test_that("compute_single_doe_linear: returns doe_results object", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  expect_s3_class(res, "doe_results")
})

test_that("compute_single_doe_linear: result has D, O, E, DOE_score, errors", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  expect_true(all(c("D", "O", "E", "DOE_score", "errors") %in% names(res)))
})

test_that("compute_single_doe_linear: D sub-list has D_early and D_term", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  expect_named(res$D, c("D_early", "D_term", "D_comp"), ignore.order = TRUE)
})

test_that("compute_single_doe_linear: O sub-list has scalar O", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  expect_length(res$O$O, 1)
})

test_that("compute_single_doe_linear: E sub-list has E_early, E_term, E_comp", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  expect_true(all(c("E_early", "E_term", "E_comp") %in% names(res$E)))
})

test_that("compute_single_doe_linear: DOE_score is scalar in [0, 1]", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  expect_length(res$DOE_score, 1)
  if (!is.na(res$DOE_score)) {
    expect_gte(res$DOE_score, 0)
    expect_lte(res$DOE_score, 1)
  }
})

test_that("compute_single_doe_linear: DOE_score is mean of available components", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  comp_vec <- c(res$D$D_comp, res$O$O, res$E$E_comp)
  expected <- mean(comp_vec, na.rm = TRUE)
  expect_equal(res$DOE_score, expected, tolerance = 1e-10)
})

test_that("compute_single_doe_linear: pseudotime length mismatch throws error", {
  f <- make_linear_fixture()
  expect_error(
    compute_single_doe_linear(
      f$expr, f$pseudotime[-1], f$early_markers, f$term_markers,
      E_method = "clusters", plot_E = FALSE
    )
  )
})

test_that("compute_single_doe_linear: unnamed pseudotime of correct length works", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt_unnamed <- unname(f$pseudotime)
  expect_no_error(
    suppressWarnings(
      compute_single_doe_linear(
        f$expr, pt_unnamed, f$early_markers, f$term_markers,
        cluster_labels    = unname(cl),
        early_clusters    = "early_cl",
        terminal_clusters = "term_cl",
        E_method = "clusters", plot_E = FALSE
      )
    )
  )
})

test_that("compute_single_doe_linear: data.frame input is accepted", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  df <- as.data.frame(f$expr)
  expect_no_error(
    suppressWarnings(
      compute_single_doe_linear(
        df, f$pseudotime, f$early_markers, f$term_markers,
        cluster_labels    = cl,
        early_clusters    = "early_cl",
        terminal_clusters = "term_cl",
        E_method = "clusters", plot_E = FALSE
      )
    )
  )
})

# =============================================================================
# compute_multi_doe_linear
# =============================================================================

test_that("compute_multi_doe_linear: returns multi_doe_results object", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_s3_class(res, "multi_doe_results")
})

test_that("compute_multi_doe_linear: has results, comparison_summary, best_trajectory", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_true(all(c("results", "comparison_summary",
                    "best_trajectory", "n_trajectories") %in% names(res)))
})

test_that("compute_multi_doe_linear: n_trajectories matches input list length", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_equal(res$n_trajectories, 2)
})

test_that("compute_multi_doe_linear: comparison_summary has one row per trajectory", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_equal(nrow(res$comparison_summary), 2)
})

test_that("compute_multi_doe_linear: best_trajectory is one of the input names", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_true(res$best_trajectory %in% c("traj_A", "traj_B"))
})

test_that("compute_multi_doe_linear: well-ordered trajectory scores higher than reversed", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt_rev <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(good = f$pseudotime, bad = pt_rev),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_equal(res$best_trajectory, "good")
})

test_that("compute_multi_doe_linear: auto-names unnamed pseudotime_list", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(f$pseudotime),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_true(grepl("trajectory_1", res$best_trajectory))
})

test_that("compute_multi_doe_linear: non-list pseudotime_list throws error", {
  f <- make_linear_fixture()
  expect_error(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list  = f$pseudotime,   # vector, not list
      early_markers    = f$early_markers,
      terminal_markers = f$term_markers,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    ),
    "list"
  )
})

# ── plot.doe_results ──────────────────────────────────────────────────────────

test_that("plot.doe_results bar: returns ggplot object", {
  skip_if_not_installed("ggplot2")
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  p <- plot(res, type = "bar")
  expect_s3_class(p, "ggplot")
})

# ── Error propagation ─────────────────────────────────────────────────────────

test_that("compute_single_doe_linear: absent markers produce NA components, not error", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime,
      early_markers     = c("nonexistent_A", "nonexistent_B"),
      terminal_markers  = c("nonexistent_C", "nonexistent_D"),
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  # D and O should be NA (no marker overlap), E should still compute from clusters
  expect_true(is.na(res$D$D_early))
  expect_true(is.na(res$D$D_term))
  expect_false(is.na(res$E$E_early))
})

test_that("compute_single_doe_linear: DOE_score is mean of non-NA components only", {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  res <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime,
      early_markers     = c("nonexistent_A", "nonexistent_B"),
      terminal_markers  = c("nonexistent_C", "nonexistent_D"),
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  comp_vec <- c(res$D$D_comp, res$O$O, res$E$E_comp)
  expected <- mean(comp_vec, na.rm = TRUE)
  if (!is.na(res$DOE_score)) {
    expect_equal(res$DOE_score, expected, tolerance = 1e-10)
  }
})

# ── Parallel execution ────────────────────────────────────────────────────────

test_that("compute_multi_doe_linear: parallel=TRUE returns same structure as sequential", {
  skip_on_cran()
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))

  res_seq <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE,
      parallel = FALSE
    )
  )
  res_par <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE,
      parallel = TRUE, n_cores = 2
    )
  )
  expect_s3_class(res_par, "multi_doe_results")
  expect_equal(res_par$n_trajectories, res_seq$n_trajectories)
  expect_equal(res_par$best_trajectory, res_seq$best_trajectory)
})

test_that("plot.multi_doe_results bar: returns ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  p <- plot(res, type = "bar")
  expect_s3_class(p, "ggplot")
})

# ── plot.doe_results additional types ─────────────────────────────────────────

.make_single_res <- function() {
  f  <- make_linear_fixture()
  cl <- make_cluster_labels(f)
  suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
}

test_that("plot.doe_results heatmap: returns ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  p <- plot(.make_single_res(), type = "heatmap")
  expect_s3_class(p, "ggplot")
})

test_that("plot.doe_results radar: runs without error (returns invisible)", {
  skip_if_not_installed("fmsb")
  skip_if_not_installed("scales")
  expect_no_error(plot(.make_single_res(), type = "radar"))
})

test_that("plot.doe_results invalid type: throws error", {
  expect_error(plot(.make_single_res(), type = "volcano"), "arg")
})

# ── plot.multi_doe_results additional types ───────────────────────────────────

.make_multi_res <- function() {
  f   <- make_linear_fixture()
  cl  <- make_cluster_labels(f)
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))
  suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
}

test_that("plot.multi_doe_results heatmap: returns ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  p <- plot(.make_multi_res(), type = "heatmap")
  expect_s3_class(p, "ggplot")
})

test_that("plot.multi_doe_results radar: runs without error", {
  skip_if_not_installed("fmsb")
  expect_no_error(plot(.make_multi_res(), type = "radar"))
})

# ── print / summary for DOE result objects ────────────────────────────────────

test_that("print.doe_results: produces output without error", {
  res <- .make_single_res()
  # DOE result objects use the default list print — just verify no error
  expect_no_error(capture.output(print(res)))
})

test_that("print.multi_doe_results: produces output without error", {
  res <- .make_multi_res()
  expect_no_error(capture.output(print(res)))
})
