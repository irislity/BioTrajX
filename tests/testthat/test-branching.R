# =============================================================================
# test-branching.R  —  Unit tests for compute_single_doe_branched() and
#                       compute_multi_doe_branched()
# =============================================================================
library(testthat)

# =============================================================================
# compute_single_doe_branched
# =============================================================================

test_that("compute_single_doe_branched: returns a list with branches and aggregate_DOE", {
  f   <- make_branched_fixture()
  res <- suppressWarnings(
    compute_single_doe_branched(
      f$expr, f$pseudotime,
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_type(res, "list")
  expect_true(all(c("branches", "aggregate_DOE", "weights",
                    "n_cells_total") %in% names(res)))
})

test_that("compute_single_doe_branched: branches list has one entry per branch", {
  f   <- make_branched_fixture()
  res <- suppressWarnings(
    compute_single_doe_branched(
      f$expr, f$pseudotime,
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_equal(sort(names(res$branches)), sort(c("AB", "AC")))
})

test_that("compute_single_doe_branched: each branch has D, O, E, DOE_score", {
  f   <- make_branched_fixture()
  res <- suppressWarnings(
    compute_single_doe_branched(
      f$expr, f$pseudotime,
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  for (b in names(res$branches)) {
    br <- res$branches[[b]]
    expect_true(all(c("D", "O", "E", "DOE_score") %in% names(br)),
                info = paste("Branch", b, "missing fields"))
  }
})

test_that("compute_single_doe_branched: per-branch DOE_score is in [0, 1]", {
  f   <- make_branched_fixture()
  res <- suppressWarnings(
    compute_single_doe_branched(
      f$expr, f$pseudotime,
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  for (b in names(res$branches)) {
    s <- res$branches[[b]]$DOE_score
    if (!is.na(s)) {
      expect_gte(s, 0, label = paste("Branch", b, "DOE_score >= 0"))
      expect_lte(s, 1, label = paste("Branch", b, "DOE_score <= 1"))
    }
  }
})

test_that("compute_single_doe_branched: aggregate_DOE is weighted mean of branch scores", {
  f   <- make_branched_fixture()
  res <- suppressWarnings(
    compute_single_doe_branched(
      f$expr, f$pseudotime,
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  w <- res$weights
  s <- vapply(res$branches, function(x) x$DOE_score, numeric(1))
  expected <- stats::weighted.mean(s, w, na.rm = TRUE)
  expect_equal(res$aggregate_DOE, expected, tolerance = 1e-10)
})

test_that("compute_single_doe_branched: n_cells_total equals sum of branch cell counts", {
  f   <- make_branched_fixture()
  res <- suppressWarnings(
    compute_single_doe_branched(
      f$expr, f$pseudotime,
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_equal(res$n_cells_total, sum(res$weights))
})

test_that("compute_single_doe_branched: unnamed branch lists throw error", {
  f <- make_branched_fixture()
  expect_error(
    compute_single_doe_branched(
      f$expr, f$pseudotime,
      early_markers_list    = unname(f$early_markers_list),
      terminal_markers_list = unname(f$term_markers_list),
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    ),
    "named"
  )
})

test_that("compute_single_doe_branched: branch filter leaving < min_cells throws error", {
  f <- make_branched_fixture()
  bad_filters <- list(
    AB = list(include = c("A", "B"), min_cells = 1000),  # impossible threshold
    AC = list(include = c("A", "C"))
  )
  expect_error(
    compute_single_doe_branched(
      f$expr, f$pseudotime,
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = bad_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
})

# =============================================================================
# compute_multi_doe_branched
# =============================================================================

test_that("compute_multi_doe_branched: returns multi_doe_branched object", {
  f    <- make_branched_fixture()
  pt2  <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res  <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_s3_class(res, "multi_doe_branched")
})

test_that("compute_multi_doe_branched: has all required top-level fields", {
  f    <- make_branched_fixture()
  pt2  <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res  <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expected_fields <- c("results", "comparison_overall", "best_overall",
                       "comparison_by_branch", "best_per_branch",
                       "branch_comparisons", "branch_data")
  expect_true(all(expected_fields %in% names(res)))
})

test_that("compute_multi_doe_branched: comparison_overall has one row per trajectory", {
  f    <- make_branched_fixture()
  pt2  <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res  <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_equal(nrow(res$comparison_overall), 2)
})

test_that("compute_multi_doe_branched: best_overall is one of the trajectory names", {
  f    <- make_branched_fixture()
  pt2  <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res  <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_true(res$best_overall %in% c("traj_A", "traj_B"))
})

test_that("compute_multi_doe_branched: best_per_branch has one row per branch", {
  f    <- make_branched_fixture()
  pt2  <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res  <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_equal(nrow(res$best_per_branch), 2)  # branches A and B
})

test_that("compute_multi_doe_branched: auto-names unnamed pseudotime_list", {
  f   <- make_branched_fixture()
  res <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(f$pseudotime),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_true(grepl("trajectory_1", res$best_overall))
})

# ── plot.multi_doe_branched ───────────────────────────────────────────────────

test_that("plot.multi_doe_branched overall bar: returns ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  f    <- make_branched_fixture()
  pt2  <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res  <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  p <- plot(res, scope = "overall", type = "bar")
  expect_s3_class(p, "ggplot")
})

test_that("plot.multi_doe_branched branch heatmap: returns ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  f    <- make_branched_fixture()
  pt2  <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res  <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  p <- plot(res, scope = "branch", type = "heatmap", branch_mode = "stack")
  expect_s3_class(p, "ggplot")
})

# ── plot.multi_doe_branched additional types / scopes ────────────────────────

.make_branched_multi <- function() {
  f   <- make_branched_fixture()
  pt2 <- setNames(rev(f$pseudotime), names(f$pseudotime))
  suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(traj_A = f$pseudotime, traj_B = pt2),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
}

test_that("plot.multi_doe_branched overall heatmap: returns ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  p <- plot(.make_branched_multi(), scope = "overall", type = "heatmap")
  expect_s3_class(p, "ggplot")
})

test_that("plot.multi_doe_branched branch bar facet: returns ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  p <- plot(.make_branched_multi(), scope = "branch", type = "bar",
            branch_mode = "facet")
  expect_s3_class(p, "ggplot")
})

test_that("plot.multi_doe_branched branch radar: runs without error", {
  skip_if_not_installed("fmsb")
  skip_if_not_installed("scales")
  expect_no_error(
    plot(.make_branched_multi(), scope = "branch", type = "radar",
         branch_mode = "separate")
  )
})

# ── Branched biological correctness ───────────────────────────────────────────

test_that("compute_multi_doe_branched: correct trajectory scores higher than reversed", {
  f      <- make_branched_fixture()
  pt_rev <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(good = f$pseudotime, bad = pt_rev),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  # The correctly-ordered trajectory should have higher aggregate_DOE
  good_score <- res$comparison_overall$aggregate_DOE[
    res$comparison_overall$trajectory == "good"]
  bad_score  <- res$comparison_overall$aggregate_DOE[
    res$comparison_overall$trajectory == "bad"]
  if (!is.na(good_score) && !is.na(bad_score)) {
    expect_gte(good_score, bad_score)
  }
})

test_that("compute_multi_doe_branched: best_overall is the correct trajectory", {
  f      <- make_branched_fixture()
  pt_rev <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res <- suppressWarnings(
    compute_multi_doe_branched(
      f$expr,
      pseudotime_list       = list(good = f$pseudotime, bad = pt_rev),
      early_markers_list    = f$early_markers_list,
      terminal_markers_list = f$term_markers_list,
      cluster_labels        = f$cluster_labels,
      branch_filters        = f$branch_filters,
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )
  expect_equal(res$best_overall, "good")
})

# ── Mismatched branch names ───────────────────────────────────────────────────

test_that("compute_single_doe_branched: mismatched branch names between early and terminal lists errors", {
  f <- make_branched_fixture()
  # Give terminal_markers_list a different branch name
  bad_term <- f$term_markers_list
  names(bad_term) <- c("A", "C")   # "C" instead of "B"
  expect_error(
    suppressWarnings(
      compute_single_doe_branched(
        f$expr, f$pseudotime,
        early_markers_list    = f$early_markers_list,   # A, B
        terminal_markers_list = bad_term,                # A, C
        cluster_labels        = f$cluster_labels,
        branch_filters        = f$branch_filters,
        E_method = "clusters", plot_E = FALSE, verbose = FALSE
      )
    )
  )
})
