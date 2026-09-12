# =============================================================================
# test-metric_directionality.R  —  Unit tests for metrics_d()
# =============================================================================
library(testthat)

# ── Return structure ──────────────────────────────────────────────────────────

test_that("metrics_d: returns a list with required names", {  
  f   <- make_linear_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_type(res, "list")
  expect_named(res, c("D_term", "D_early", "s_i_early", "s_i_term"),
               ignore.order = TRUE)
})

test_that("metrics_d: D_early and D_term are scalars in [0, 1]", {  
  f   <- make_linear_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_length(res$D_early, 1)
  expect_length(res$D_term,  1)
  expect_gte(res$D_early, 0); expect_lte(res$D_early, 1)
  expect_gte(res$D_term,  0); expect_lte(res$D_term,  1)
})

test_that("metrics_d: per-cell score vectors have length == ncol(expr)", {  
  f   <- make_linear_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_length(res$s_i_early, ncol(f$expr))
  expect_length(res$s_i_term,  ncol(f$expr))
})

# ── Biological direction ──────────────────────────────────────────────────────

test_that("metrics_d: D_early is HIGH for well-ordered early markers (decreasing)", { # nolint
  # early genes decrease along pseudotime → rho_early < 0 → D_early close to 1
  # Threshold is 0.65 (vs 0.75 for Gaussian) to accommodate ZINB overdispersion.
  f   <- make_linear_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gt(res$D_term,  0.65)
  expect_gt(res$D_early, 0.65)
})

test_that("metrics_d: D_term is HIGH for well-ordered terminal markers (increasing)", { # nolint
  f   <- make_linear_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gt(res$D_term, 0.65)
})

test_that("metrics_d: reversed pseudotime flips D_early and D_term", {
  f    <- make_linear_fixture()
  res1 <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  res2 <- metrics_d(f$expr, f$early_markers, f$term_markers, rev(f$pseudotime))
  # Reversing pseudotime should flip the Spearman correlations and create D_early and D_term with lower score
  expect_lt(res2$D_early, res1$D_early)
  expect_lt(res2$D_term,  res1$D_term)
})

test_that("metrics_d: flat expression gives D near 0 (ReLU clips random noise)", {
  # With ReLU formulation D = max(0, rho), flat noise yields rho ≈ 0 so D ≈ 0
  f   <- make_flat_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_lt(res$D_early, 0.2)
  expect_lt(res$D_term,  0.2)
})

# ── Edge cases ────────────────────────────────────────────────────────────────

test_that("metrics_d: warns and returns NA scores when no markers overlap", { 
  f <- make_linear_fixture()
  expect_warning(
    res <- metrics_d(f$expr, c("FAKE1", "FAKE2"), f$term_markers, f$pseudotime),
    "No markers found"
  )
  expect_true(all(is.na(res$s_i_early)))
})

test_that("metrics_d: partial marker overlap uses only present genes", { 
  f <- make_linear_fixture()
  mixed_early <- c(f$early_markers[1:5], "NOTREAL1", "NOTREAL2")
  # Should not error, just use the 5 real genes
  expect_no_error(
    res <- metrics_d(f$expr, mixed_early, f$term_markers, f$pseudotime)
  )
  expect_false(is.na(res$D_early))
})

test_that("metrics_d: single-cell input does not crash", {  
  f <- make_linear_fixture()
  expect_no_error(
    metrics_d(f$expr[, 1, drop = FALSE],
              f$early_markers, f$term_markers,
              f$pseudotime[1])
  )
})

test_that("metrics_d_early: per-cell scores are means of marker expression", {  
  f     <- make_linear_fixture()
  genes <- intersect(f$early_markers, rownames(f$expr))
  expected_early <- colSums(f$expr[genes, , drop = FALSE]) / length(genes)
  res   <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_equal(res$s_i_early, expected_early, tolerance = 1e-10)
})


test_that("metrics_d_term: per-cell scores are means of marker expression", {
  f     <- make_linear_fixture()
  genes <- intersect(f$term_markers, rownames(f$expr))
  expected_term <- colSums(f$expr[genes, , drop = FALSE]) / length(genes)
  res   <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_equal(res$s_i_term, expected_term, tolerance = 1e-10)
})

# ── plot_metrics_d ────────────────────────────────────────────────────────────

test_that("plot_metrics_d: runs without error on valid input", {
  f   <- make_linear_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_no_error(
    suppressWarnings(plot_metrics_d(f$pseudotime, res))
  )
})

test_that("plot_metrics_d: returns invisible NULL when all scores are NA", {
  f   <- make_linear_fixture()
  res <- suppressWarnings(
    metrics_d(f$expr, c("FAKE1"), c("FAKE2"), f$pseudotime)
  )
  expect_warning(out <- plot_metrics_d(f$pseudotime, res), "NA")
  expect_null(out)
})

# ── Ceiling tests (Gaussian fixture) ─────────────────────────────────────────
# Gaussian noise preserves the full score range, so D should be near 1.
# These tests verify mathematical correctness, not just direction.

test_that("metrics_d ceiling: D_early near 1 on Gaussian linear fixture", {
  f   <- make_gaussian_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gte(res$D_early, 0.95)
})

test_that("metrics_d ceiling: D_term near 1 on Gaussian linear fixture", {
  f   <- make_gaussian_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gte(res$D_term, 0.95)
})

test_that("metrics_d floor: D_early near 0 on Gaussian reversed fixture", {
  f   <- make_gaussian_reversed_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_lte(res$D_early, 0.05)
})

test_that("metrics_d floor: D_term near 0 on Gaussian reversed fixture", {
  f   <- make_gaussian_reversed_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_lte(res$D_term, 0.05)
})

test_that("metrics_d ceiling: D_early + D_term near 2 on Gaussian linear fixture", {
  f   <- make_gaussian_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gte(res$D_early + res$D_term, 1.9)
})

test_that("metrics_d ceiling: D_early + D_term near 0 on Gaussian reversed fixture", {
  f   <- make_gaussian_reversed_fixture()
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_lte(res$D_early + res$D_term, 0.1)
})
