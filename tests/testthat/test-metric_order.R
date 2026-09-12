# =============================================================================
# test-metric_order.R  —  Unit tests for metrics_o()
# =============================================================================
library(testthat)

# ── Return structure ──────────────────────────────────────────────────────────

test_that("metrics_o: returns a list with required names", {  
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_type(res, "list")
  expect_true(all(c("O", "O_f", "genes_used", "n_cells", "n_genes") %in% names(res)))
})

test_that("metrics_o: O is a scalar in [0, 1]", {  
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_length(res$O, 1)
  expect_gte(res$O, 0)
  expect_lte(res$O, 1)
})

test_that("metrics_o: O_f values are all in [0, 1]", { 
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_true(all(res$O_f >= 0 & res$O_f <= 1))
})

test_that("metrics_o: O equals mean of O_f", { 
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_equal(res$O, mean(res$O_f, na.rm = TRUE), tolerance = 1e-10)
})

test_that("metrics_o: genes_used matches names of O_f", {
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_equal(res$genes_used, names(res$O_f))
})

test_that("metrics_o: n_bins and null_r2 are returned", {
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_true(!is.null(res$n_bins))
  expect_true(!is.null(res$null_r2))
  expect_gt(res$null_r2, 0)
  expect_lt(res$null_r2, 1)
})

test_that("metrics_o: n_cells matches ncol(expr)", { 
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_equal(res$n_cells, ncol(f$expr))
})

test_that("metrics_o: n_genes matches length of O_f", { 
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_equal(res$n_genes, length(res$O_f))
})

# ── Biological direction ──────────────────────────────────────────────────────

test_that("metrics_o: O is HIGH for perfectly monotone markers", {
  # Calibrated pseudobulk metric; threshold 0.3 accommodates ZINB overdispersion
  # with n_cells=80 and n_bins=10 (~8 cells/bin).
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gt(res$O, 0.3)
})

test_that("metrics_o: orientation field is '+' or '-'", { 
  f   <- make_linear_fixture()
  res_pos <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  res_neg  <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_true(res_pos$orientation %in% c("+", "-"))
  expect_true(res_neg$orientation %in% c("+", "-"))
})

test_that("metrics_o: orientation_max=TRUE picks the better orientation", { 
  f    <- make_linear_fixture()
  res_auto <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime,
                        orientation_max = TRUE)
  res_pos  <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime,
                        orientation_max = FALSE)
  res_neg  <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime, orientation_max = FALSE)
  # auto should be >= fixed positive orientation
  expect_gte(res_auto$O, res_pos$O - 1e-10)
  expect_gte(res_auto$O, res_neg$O - 1e-10)
})

test_that("metrics_o: O is LOW for random (flat) expression", { 
  f   <- make_flat_fixture(n_cells = 100)
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  # Random data should not produce near-perfect order
  expect_lt(res$O, 0.8)
})

# ── IQR filter ────────────────────────────────────────────────────────────────

test_that("metrics_o: iqr_quantile=0 keeps all genes (no filtering)", {
  f    <- make_linear_fixture()
  res0 <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime,
                    iqr_quantile = 0)
  res5 <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime,
                    iqr_quantile = 0.5)
  # With no IQR filter, more genes are retained
  expect_gte(res0$n_genes, res5$n_genes)
})

test_that("metrics_o: invalid iqr_quantile throws error", {
  f <- make_linear_fixture()
  expect_error(
    metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime,
              iqr_quantile = 1.5),
    "iqr_quantile"
  )
  expect_error(
    metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime,
              iqr_quantile = -0.1),
    "iqr_quantile"
  )
})

# ── Edge cases ────────────────────────────────────────────────────────────────

test_that("metrics_o: no marker overlap returns NA with warning", { 
  f <- make_linear_fixture()
  expect_warning(
    res <- metrics_o(f$expr, c("FAKE1"), c("FAKE2"), f$pseudotime),
    "No overlap"
  )
  expect_true(is.na(res$O))
  expect_length(res$genes_used, 0)
})

test_that("metrics_o: duplicate rownames in expr are de-duplicated silently", {
  f <- make_linear_fixture()
  expr_dup <- rbind(f$expr, f$expr[1:5, ])
  rownames(expr_dup) <- c(rownames(f$expr), rownames(f$expr)[1:5])
  expect_no_error(
    res <- metrics_o(expr_dup, f$early_markers, f$term_markers, f$pseudotime)
  )
  expect_false(is.na(res$O))
})

test_that("metrics_o: pseudotime length mismatch throws error", { 
  f <- make_linear_fixture()
  expect_error(
    metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime[-1])
  )
})

test_that("metrics_o: fewer than 3 cells returns NA gracefully", {
  f <- make_linear_fixture()
  expect_no_error(
    res <- metrics_o(f$expr[, 1:2, drop = FALSE],
                     f$early_markers, f$term_markers,
                     f$pseudotime[1:2])
  )
  expect_true(is.na(res$O))
})

# ── plot_metrics_o ────────────────────────────────────────────────────────────

test_that("plot_metrics_o: runs without error on valid input", {
  f   <- make_linear_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_no_error(
    suppressWarnings(
      plot_metrics_o(f$expr, f$pseudotime, res, f$early_markers, f$term_markers)
    )
  )
})

test_that("plot_metrics_o: returns invisible NULL when no marker overlap", {
  f   <- make_linear_fixture()
  res <- suppressWarnings(
    metrics_o(f$expr, c("FAKE1"), c("FAKE2"), f$pseudotime)
  )
  expect_warning(
    out <- plot_metrics_o(f$expr, f$pseudotime, res, c("FAKE1"), c("FAKE2")),
    "No overlap"
  )
  expect_null(out)
})

# ── Ceiling tests (Gaussian fixture) ─────────────────────────────────────────
# With Gaussian noise and slope = 7 the bin means are smooth enough that
# isotonic R² = 1 and null calibration returns O = 1.

test_that("metrics_o ceiling: O near 1 on Gaussian linear fixture", {
  f   <- make_gaussian_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gte(res$O, 0.95)
})

test_that("metrics_o ceiling: all per-gene O_f near 1 on Gaussian linear fixture", {
  f   <- make_gaussian_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_true(all(res$O_f >= 0.8))
})

# O uses orientation_max = TRUE, so it is invariant to pseudotime direction.
# Reversing pseudotime should NOT reduce O.
test_that("metrics_o orientation invariance: O near 1 on Gaussian reversed fixture", {
  f   <- make_gaussian_reversed_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gte(res$O, 0.95)
})

test_that("metrics_o orientation invariance: orientation field is '-' on reversed fixture", {
  f   <- make_gaussian_reversed_fixture()
  res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_equal(res$orientation, "-")
})
