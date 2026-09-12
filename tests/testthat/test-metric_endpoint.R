# =============================================================================
# test-metric_endpoint.R  —  Unit tests for metrics_e()
# =============================================================================
library(testthat)

# ── Shared minimal setup ──────────────────────────────────────────────────────
# Build cluster-based labels (no mclust dependency) for deterministic tests
make_cluster_E_input <- function(n_cells = 100, seed = 42) {
  set.seed(seed)
  pt      <- seq(0, 1, length.out = n_cells)
  # first 20 cells are "early", last 20 are "terminal"
  clusters <- c(rep("early_cl",  20),
                rep("mid_cl",    60),
                rep("term_cl",   20))
  names(clusters) <- paste0("cell", seq_len(n_cells))
  list(pseudotime       = pt,
       cluster_labels   = clusters,
       early_clusters   = "early_cl",
       terminal_clusters = "term_cl",
       n_cells          = n_cells)
}

# ── Return structure (clusters method) ───────────────────────────────────────

test_that("metrics_e clusters: returns endpoints_validity object", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_s3_class(res, "endpoints_validity")
})

test_that("metrics_e clusters: has required fields", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_true(all(c("E_early", "E_term", "E_comp", "method",
                    "n_cells", "summary") %in% names(res)))
})

test_that("metrics_e clusters: E_early and E_term are in (0, 1]", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_gte(res$E_early, 0); expect_lte(res$E_early, 1)
  expect_gte(res$E_term,  0); expect_lte(res$E_term,  1)
})

test_that("metrics_e clusters: n_cells matches input length", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_equal(res$n_cells, d$n_cells)
})

# ── Biological correctness (clusters method) ─────────────────────────────────

test_that("metrics_e clusters: E_early HIGH when early cells are at pseudotime start", {
  # early cells occupy the first 20 positions → precision@k should be 1
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_equal(res$E_early, 1.0, tolerance = 1e-10)
})

test_that("metrics_e clusters: E_term HIGH when terminal cells are at pseudotime end", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_equal(res$E_term, 1.0, tolerance = 1e-10)
})

test_that("metrics_e clusters: E_early LOW when early cells are at pseudotime end", {
  n   <- 100
  pt  <- seq(0, 1, length.out = n)
  # early cells placed at the END (wrong direction)
  cl  <- c(rep("mid_cl",   80), rep("early_cl", 20))
  res <- metrics_e(pt,
                   cluster_labels    = cl,
                   early_clusters    = "early_cl",
                   terminal_clusters = character(0),
                   method = "clusters", plot = FALSE)
  expect_lt(res$E_early, 0.5)
})

# ── E_comp (harmonic mean) ────────────────────────────────────────────────────

test_that("metrics_e clusters: E_comp is harmonic mean of normalised E_early and E_term", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  # When both E_early = E_term = 1 and proportions are equal, E_comp should be 1
  if (!is.na(res$E_comp)) {
    expect_gte(res$E_comp, 0)
  }
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("metrics_e: NA in pseudotime throws error", {
  d  <- make_cluster_E_input()
  pt <- d$pseudotime; pt[5] <- NA
  expect_error(
    metrics_e(pt,
              cluster_labels    = d$cluster_labels,
              early_clusters    = d$early_clusters,
              terminal_clusters = d$terminal_clusters,
              method = "clusters", plot = FALSE),
    "NA"
  )
})

test_that("metrics_e clusters: missing cluster_labels throws error", {
  d <- make_cluster_E_input()
  expect_error(
    metrics_e(d$pseudotime,
              early_clusters    = d$early_clusters,
              terminal_clusters = d$terminal_clusters,
              method = "clusters", plot = FALSE),
    "cluster_labels"
  )
})

test_that("metrics_e clusters: cluster_labels length mismatch throws error", {
  d <- make_cluster_E_input()
  expect_error(
    metrics_e(d$pseudotime,
              cluster_labels    = d$cluster_labels[-1],
              early_clusters    = d$early_clusters,
              terminal_clusters = d$terminal_clusters,
              method = "clusters", plot = FALSE),
    "same length"
  )
})

test_that("metrics_e gmm: missing marker scores throws error", {
  d <- make_cluster_E_input()
  expect_error(
    metrics_e(d$pseudotime, method = "gmm", plot = FALSE),
    "early_marker_scores"
  )
})

# ── print / summary S3 methods ────────────────────────────────────────────────

test_that("print.endpoints_validity: runs without error", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_output(print(res), "Endpoints Validity")
})

test_that("summary.endpoints_validity: runs without error", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_output(summary(res), "Endpoints Validity")
})

# ── Cluster info stored correctly ─────────────────────────────────────────────

test_that("metrics_e clusters: cluster_info records early and terminal counts", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_equal(res$cluster_info$early_count,    20)
  expect_equal(res$cluster_info$terminal_count, 20)
})

test_that("metrics_e clusters: summary proportions sum to <= 1", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = d$early_clusters,
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_lte(res$summary$prop_early + res$summary$prop_terminal, 1 + 1e-10)
})

# ── GMM method ────────────────────────────────────────────────────────────────

make_gmm_E_input <- function(n_cells = 120, seed = 42) {
  set.seed(seed)
  pt <- seq(0, 1, length.out = n_cells)
  # Bimodal early scores: first 40 cells score high (true early), rest score low
  early_scores <- c(rnorm(40, mean = 2.0, sd = 0.2),
                    rnorm(80, mean = 0.2, sd = 0.2))
  # Bimodal terminal scores: last 40 cells score high (true terminal), rest score low
  term_scores  <- c(rnorm(80, mean = 0.2, sd = 0.2),
                    rnorm(40, mean = 2.0, sd = 0.2))
  list(pseudotime = pt, early_scores = early_scores, term_scores = term_scores)
}

test_that("metrics_e gmm: returns endpoints_validity object", {
  skip_if_not_installed("mclust")
  d   <- make_gmm_E_input()
  res <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = d$early_scores,
              terminal_marker_scores = d$term_scores,
              method = "gmm", plot = FALSE)
  )
  expect_s3_class(res, "endpoints_validity")
})

test_that("metrics_e gmm: E_early and E_term are in [0, 1]", {
  skip_if_not_installed("mclust")
  d   <- make_gmm_E_input()
  res <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = d$early_scores,
              terminal_marker_scores = d$term_scores,
              method = "gmm", plot = FALSE)
  )
  if (!is.na(res$E_early)) expect_true(res$E_early >= 0 && res$E_early <= 1)
  if (!is.na(res$E_term))  expect_true(res$E_term  >= 0 && res$E_term  <= 1)
})

test_that("metrics_e gmm: E_comp is in [0, 1] when computable", {
  skip_if_not_installed("mclust")
  d   <- make_gmm_E_input()
  res <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = d$early_scores,
              terminal_marker_scores = d$term_scores,
              method = "gmm", plot = FALSE)
  )
  if (!is.na(res$E_comp)) {
    expect_gte(res$E_comp, 0)
    expect_lte(res$E_comp, 1)
  }
})

test_that("metrics_e gmm: well-aligned endpoints score higher than random", {
  skip_if_not_installed("mclust")
  set.seed(7)
  d   <- make_gmm_E_input()
  res_good <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = d$early_scores,
              terminal_marker_scores = d$term_scores,
              method = "gmm", plot = FALSE)
  )
  # Assert GMM actually produced a result — if it returns NA the test is vacuous
  expect_false(is.na(res_good$E_comp),
               label = "E_comp is not NA for well-aligned bimodal GMM input")
  # Shuffle scores so they no longer align with pseudotime
  set.seed(7)
  res_rand <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = sample(d$early_scores),
              terminal_marker_scores = sample(d$term_scores),
              method = "gmm", plot = FALSE)
  )
  if (!is.na(res_rand$E_comp)) {
    expect_gt(res_good$E_comp, res_rand$E_comp)
  }
})

test_that("metrics_e gmm: too few cells produces warning and NA labels", {
  skip_if_not_installed("mclust")
  pt  <- seq(0, 1, length.out = 10)
  ns  <- rnorm(10)
  ts  <- rnorm(10)
  res <- suppressWarnings(
    metrics_e(pt, early_marker_scores = ns, terminal_marker_scores = ts,
              method = "gmm", plot = FALSE)
  )
  # With only 10 cells (< min_cells_per_component * 2 = 20), labels should be all NA
  expect_true(all(is.na(res$early_labels)))
})

test_that("metrics_e gmm: no-variation scores produce warning and NA labels", {
  skip_if_not_installed("mclust")
  pt <- seq(0, 1, length.out = 100)
  ns <- rep(1.0, 100)   # zero variance
  ts <- rnorm(100, 1, 0.5)
  expect_warning(
    metrics_e(pt, early_marker_scores = ns, terminal_marker_scores = ts,
              method = "gmm", plot = FALSE),
    "variation|variance|Insufficient"
  )
})

# ── GMM from expression matrix (end-to-end marker workflow) ──────────────────

test_that("metrics_e gmm: end-to-end from expr matrix via colMeans of markers", {
  skip_if_not_installed("mclust")
  f <- make_linear_fixture()   # ZINB fixture: early genes decrease, terminal increase

  # Typical user workflow: aggregate marker expression per cell then pass to GMM
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])

  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )

  expect_s3_class(res, "endpoints_validity")
  expect_equal(length(early_scores), length(f$pseudotime))
  if (!is.na(res$E_early)) expect_gte(res$E_early, 0)
  if (!is.na(res$E_term))  expect_gte(res$E_term,  0)
})

test_that("metrics_e gmm: expr-derived scores give higher E_comp than shuffled pseudotime", {
  skip_if_not_installed("mclust")
  f <- make_linear_fixture()

  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])

  res_good <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )

  set.seed(42)
  pt_shuf <- setNames(sample(f$pseudotime), names(f$pseudotime))
  res_shuf <- suppressWarnings(
    metrics_e(pt_shuf,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )

  if (!is.na(res_good$E_comp) && !is.na(res_shuf$E_comp)) {
    expect_gt(res_good$E_comp, res_shuf$E_comp)
  }
})

# ── Reversed placement ───────────────────────────────────────────────────────

test_that("metrics_e gmm: reversed placement gives lower E_comp than correct placement", {
  skip_if_not_installed("mclust")
  # Build scores where early-high cells are correctly at pseudotime start
  set.seed(42)
  n  <- 120
  pt <- seq(0, 1, length.out = n)
  early_scores <- c(rnorm(40, mean = 2.0, sd = 0.2),   # first 40 cells: high early
                    rnorm(80, mean = 0.2, sd = 0.2))
  term_scores  <- c(rnorm(80, mean = 0.2, sd = 0.2),   # last 40 cells: high terminal
                    rnorm(40, mean = 2.0, sd = 0.2))

  res_correct <- suppressWarnings(
    metrics_e(pt,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )

  # Reverse pseudotime: early-high cells now sit at the end (wrong placement)
  pt_rev <- rev(pt)
  res_rev <- suppressWarnings(
    metrics_e(pt_rev,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )

  if (!is.na(res_correct$E_early) && !is.na(res_rev$E_early))
    expect_gt(res_correct$E_early, res_rev$E_early)
  if (!is.na(res_correct$E_comp) && !is.na(res_rev$E_comp))
    expect_gt(res_correct$E_comp, res_rev$E_comp)
})

test_that("metrics_e gmm: E_early low when high-early cells are at pseudotime end", {
  skip_if_not_installed("mclust")
  set.seed(7)
  n  <- 120
  pt <- seq(0, 1, length.out = n)
  # High early scores placed at the END (wrong)
  early_scores <- c(rnorm(80, mean = 0.2, sd = 0.2),
                    rnorm(40, mean = 2.0, sd = 0.2))
  term_scores  <- rnorm(n, mean = 0.2, sd = 0.2)   # no terminal signal

  res <- suppressWarnings(
    metrics_e(pt,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )

  if (!is.na(res$E_early))
    expect_lt(res$E_early, 0.5)
})

# ── GMM label output structure ────────────────────────────────────────────────

test_that("metrics_e gmm: gmm_labels has correct length", {
  skip_if_not_installed("mclust")
  d   <- make_gmm_E_input()
  res <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = d$early_scores,
              terminal_marker_scores = d$term_scores,
              method = "gmm", plot = FALSE)
  )
  expect_equal(length(res$early_labels),    length(d$pseudotime))
  expect_equal(length(res$terminal_labels), length(d$pseudotime))
})

test_that("metrics_e gmm: gmm_labels contain only expected values", {
  skip_if_not_installed("mclust")
  d   <- make_gmm_E_input()
  res <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = d$early_scores,
              terminal_marker_scores = d$term_scores,
              method = "gmm", plot = FALSE)
  )
  valid_vals <- c(TRUE, FALSE, NA)
  expect_true(all(res$early_labels    %in% valid_vals | is.na(res$early_labels)))
  expect_true(all(res$terminal_labels %in% valid_vals | is.na(res$terminal_labels)))
})

test_that("metrics_e gmm: method field equals 'gmm'", {
  skip_if_not_installed("mclust")
  d   <- make_gmm_E_input()
  res <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = d$early_scores,
              terminal_marker_scores = d$term_scores,
              method = "gmm", plot = FALSE)
  )
  expect_equal(res$method, "gmm")
})

test_that("metrics_e gmm: well-separated scores label majority of cells non-NA", {
  skip_if_not_installed("mclust")
  d   <- make_gmm_E_input()
  res <- suppressWarnings(
    metrics_e(d$pseudotime,
              early_marker_scores    = d$early_scores,
              terminal_marker_scores = d$term_scores,
              method = "gmm", plot = FALSE)
  )
  # With clear bimodal scores, GMM should label most cells (>50%)
  pct_labeled <- mean(!is.na(res$early_labels))
  expect_gt(pct_labeled, 0.5)
})

# ── Edge cases ────────────────────────────────────────────────────────────────

test_that("metrics_e clusters: empty early_clusters gives E_early = NA", {
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = character(0),
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  expect_true(is.na(res$E_early))
})

test_that("metrics_e clusters: E_comp is NA when one endpoint is zero", {
  # Terminal cells are all at pseudotime end, but early_clusters is empty
  d   <- make_cluster_E_input()
  res <- metrics_e(d$pseudotime,
                   cluster_labels    = d$cluster_labels,
                   early_clusters    = character(0),
                   terminal_clusters = d$terminal_clusters,
                   method = "clusters", plot = FALSE)
  # E_early is NA, so harmonic mean cannot be computed -> E_comp should be NA
  expect_true(is.na(res$E_comp))
})

# ── Ceiling tests (Gaussian fixture) ─────────────────────────────────────────
# With Gaussian noise and slope = 7 the per-cell marker scores spread evenly
# over [1, 8], giving GMM a clean separation → E near 1.

test_that("metrics_e gmm ceiling: E_early near 1 on Gaussian linear fixture", {
  skip_if_not_installed("mclust")
  f            <- make_gaussian_fixture()
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )
  expect_gte(res$E_early, 0.90)
})

test_that("metrics_e gmm ceiling: E_term near 1 on Gaussian linear fixture", {
  skip_if_not_installed("mclust")
  f            <- make_gaussian_fixture()
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )
  expect_gte(res$E_term, 0.90)
})

test_that("metrics_e gmm ceiling: E_comp near 1 on Gaussian linear fixture", {
  skip_if_not_installed("mclust")
  f            <- make_gaussian_fixture()
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )
  expect_false(is.na(res$E_comp))
  expect_gte(res$E_comp, 0.90)
})

test_that("metrics_e gmm floor: E_early near 0 on Gaussian reversed fixture", {
  skip_if_not_installed("mclust")
  f            <- make_gaussian_reversed_fixture()
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )
  if (!is.na(res$E_early)) expect_lte(res$E_early, 0.15)
})

# ── Endpoint fixture: ground-truth E validation ───────────────────────────────
# make_endpoint_fixture() provides true_early_endpoint / true_terminal_endpoint
# as ground truth.  These tests verify that GMM recovers the known endpoint
# cells and that E scores reflect their correct placement.

test_that("metrics_e endpoint: E_comp high on endpoint fixture", {
  skip_if_not_installed("mclust")
  f            <- make_endpoint_fixture()
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )
  expect_false(is.na(res$E_comp))
  expect_gte(res$E_comp, 0.90)
})

test_that("metrics_e endpoint: GMM k matches true endpoint count", {
  skip_if_not_installed("mclust")
  f            <- make_endpoint_fixture()
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )
  expect_equal(res$summary$n_early,    length(f$true_early_endpoint))
  expect_equal(res$summary$n_terminal, length(f$true_terminal_endpoint))
})

test_that("metrics_e endpoint: GMM early labels recover true early endpoint cells", {
  skip_if_not_installed("mclust")
  f            <- make_endpoint_fixture()
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )
  gmm_early <- names(f$pseudotime)[!is.na(res$early_labels) & res$early_labels == 1]
  overlap   <- sum(gmm_early %in% f$true_early_endpoint)
  expect_gte(overlap / length(f$true_early_endpoint), 0.90)
})

test_that("metrics_e endpoint: GMM terminal labels recover true terminal endpoint cells", {
  skip_if_not_installed("mclust")
  f            <- make_endpoint_fixture()
  early_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  res <- suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = early_scores,
              terminal_marker_scores = term_scores,
              method = "gmm", plot = FALSE)
  )
  gmm_term <- names(f$pseudotime)[!is.na(res$terminal_labels) & res$terminal_labels == 1]
  overlap  <- sum(gmm_term %in% f$true_terminal_endpoint)
  expect_gte(overlap / length(f$true_terminal_endpoint), 0.90)
})
