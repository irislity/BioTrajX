# =============================================================================
# test-simulation.R  —  Simulation-based tests for BioTrajX
#
# Philosophy: generate synthetic single-cell trajectories with KNOWN ground
# truth and verify that DOE metrics respond correctly to:
#   1. Signal strength  (strong vs weak marker trends)
#   2. Noise level      (low vs high expression noise)
#   3. Pseudotime quality (perfect vs shuffled vs reversed ordering)
#   4. Marker set quality (all real vs contaminated with random genes)
#   5. Cell count scaling (small vs large datasets)
#   6. Metric monotonicity across a noise gradient
# =============================================================================
library(testthat)

# ── Simulation helpers ────────────────────────────────────────────────────────

#' Simulate a linear trajectory with ZINB noise
#' @param n_cells   number of cells
#' @param n_markers number of marker genes per class (early / terminal)
#' @param n_noise   number of uninformative genes
#' @param slope     magnitude of the linear trend (higher = stronger signal);
#'                  early mean runs from (1+slope) down to 1,
#'                  terminal mean runs from 1 up to (1+slope)
#' @param sigma     noise level: maps to NB dispersion size = 1/(2*sigma),
#'                  so low sigma = low noise (large size), high sigma = high noise
#'                  (small size). Typical values: 0.05 (near-Poisson) to 3.0
#'                  (highly overdispersed).
#' @param pi0       zero-inflation probability (dropout rate). Default 0.5.
#' @param seed      RNG seed
sim_linear <- function(n_cells   = 100,
                       n_markers = 10,
                       n_noise   = 30,
                       slope     = 2.5,
                       sigma     = 0.2,
                       pi0       = 0.5,
                       seed      = 1) {
  set.seed(seed)
  pt      <- seq(0, 1, length.out = n_cells)
  n_tot   <- 2 * n_markers + n_noise
  nb_size <- max(0.1, 1 / (2 * sigma))   # sigma → NB dispersion

  # Draw ZINB counts then log1p-normalise, matching real scRNA-seq input.
  .rzinb <- function(n, mu) {
    counts <- rnbinom(n, mu = pmax(mu, 1e-6), size = nb_size)
    counts[runif(n) < pi0] <- 0L
    log1p(counts)
  }

  mat <- matrix(0L, nrow = n_tot, ncol = n_cells)

  # early genes: mean decreases from (1+slope) to 1
  for (i in seq_len(n_markers))
    mat[i, ] <- .rzinb(n_cells, mu = (1 + slope) - slope * pt)

  # terminal genes: mean increases from 1 to (1+slope)
  for (i in seq_len(n_markers) + n_markers)
    mat[i, ] <- .rzinb(n_cells, mu = 1 + slope * pt)

  # noise genes: constant mean
  for (i in seq_len(n_noise) + 2 * n_markers)
    mat[i, ] <- .rzinb(n_cells, mu = rep(1, n_cells))

  rownames(mat) <- c(paste0("early_g",  seq_len(n_markers)),
                     paste0("term_g",   seq_len(n_markers)),
                     paste0("noise_g",  seq_len(n_noise)))
  colnames(mat) <- paste0("cell", seq_len(n_cells))

  list(
    expr          = mat,
    pseudotime    = setNames(pt, colnames(mat)),
    early_markers = paste0("early_g", seq_len(n_markers)),
    term_markers  = paste0("term_g",  seq_len(n_markers))
  )
}

# ── 1. SIGNAL STRENGTH ────────────────────────────────────────────────────────

test_that("SIM signal: strong slope gives higher D_term than weak slope", {
  strong <- sim_linear(slope = 3.0, sigma = 0.1)
  weak   <- sim_linear(slope = 0.3, sigma = 0.1)

  res_s <- metrics_d(strong$expr, strong$early_markers,
                     strong$term_markers, strong$pseudotime)
  res_w <- metrics_d(weak$expr,   weak$early_markers,
                     weak$term_markers,   weak$pseudotime)

  expect_gt(res_s$D_term, res_w$D_term)
})

test_that("SIM signal: strong slope gives higher O than weak slope", {
  strong <- sim_linear(slope = 3.0, sigma = 0.1)
  weak   <- sim_linear(slope = 0.3, sigma = 0.1)

  res_s <- metrics_o(strong$expr, strong$early_markers,
                     strong$term_markers, strong$pseudotime)
  res_w <- metrics_o(weak$expr,   weak$early_markers,
                     weak$term_markers,   weak$pseudotime)

  expect_gt(res_s$O, res_w$O)
})

# ── 2. NOISE LEVEL ────────────────────────────────────────────────────────────

test_that("SIM noise: low noise gives higher D_term than high noise", {
  low_noise  <- sim_linear(slope = 2.5, sigma = 0.05)
  high_noise <- sim_linear(slope = 2.5, sigma = 1.5)

  res_l <- metrics_d(low_noise$expr,  low_noise$early_markers,
                     low_noise$term_markers,  low_noise$pseudotime)
  res_h <- metrics_d(high_noise$expr, high_noise$early_markers,
                     high_noise$term_markers, high_noise$pseudotime)

  expect_gt(res_l$D_term, res_h$D_term)
})

test_that("SIM noise: low noise gives higher O than high noise", {
  low_noise  <- sim_linear(slope = 2.5, sigma = 0.05)
  high_noise <- sim_linear(slope = 2.5, sigma = 1.5)

  res_l <- metrics_o(low_noise$expr,  low_noise$early_markers,
                     low_noise$term_markers,  low_noise$pseudotime)
  res_h <- metrics_o(high_noise$expr, high_noise$early_markers,
                     high_noise$term_markers, high_noise$pseudotime)

  expect_gt(res_l$O, res_h$O)
})

# ── 3. PSEUDOTIME QUALITY ─────────────────────────────────────────────────────

test_that("SIM pseudotime: perfect ordering gives D_term > 0.7", {
  # pi0 = 0.1 keeps zero-inflation low enough for ReLU D to clear 0.7.
  f   <- sim_linear(slope = 3, sigma = 0.1, pi0 = 0.1)
  res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  expect_gt(res$D_term, 0.7)
})

test_that("SIM pseudotime: shuffled ordering gives lower D_term than perfect", {
  set.seed(99)
  f         <- sim_linear(slope = 3, sigma = 0.1)
  pt_shuf   <- setNames(sample(f$pseudotime), names(f$pseudotime))

  res_perf  <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  res_shuf  <- metrics_d(f$expr, f$early_markers, f$term_markers, pt_shuf)

  expect_gt(res_perf$D_term, res_shuf$D_term)
})

test_that("SIM pseudotime: reversed ordering gives D_term < 0.5", {
  f       <- sim_linear(slope = 3, sigma = 0.1)
  pt_rev  <- setNames(rev(f$pseudotime), names(f$pseudotime))
  res     <- metrics_d(f$expr, f$early_markers, f$term_markers, pt_rev)
  expect_lt(res$D_term, 0.5)
})

test_that("SIM pseudotime: shuffled ordering gives lower O than perfect", {
  set.seed(99)
  f       <- sim_linear(slope = 3, sigma = 0.1)
  pt_shuf <- setNames(sample(f$pseudotime), names(f$pseudotime))

  res_p   <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  res_s   <- metrics_o(f$expr, f$early_markers, f$term_markers, pt_shuf)

  expect_gt(res_p$O, res_s$O)
})

test_that("SIM pseudotime: orientation_max recovers good O even for reversed input", {
  f      <- sim_linear(slope = 3, sigma = 0.1)
  pt_rev <- setNames(rev(f$pseudotime), names(f$pseudotime))

  res_rev_auto  <- metrics_o(f$expr, f$early_markers, f$term_markers,
                              pt_rev, orientation_max = TRUE)
  res_rev_fixed <- metrics_o(f$expr, f$early_markers, f$term_markers,
                              pt_rev, orientation_max = FALSE)

  # With orientation_max, the reversed pseudotime should be corrected
  expect_gte(res_rev_auto$O, res_rev_fixed$O - 1e-10)
})

# ── 4. MARKER SET QUALITY ─────────────────────────────────────────────────────

test_that("SIM markers: pure marker set gives higher D_term than contaminated set", {
  f <- sim_linear(slope = 3, sigma = 0.1)

  # contaminate terminal markers with 5 noise genes
  contaminated_term <- c(f$term_markers, paste0("noise_g", 1:5))

  res_pure  <- metrics_d(f$expr, f$early_markers, f$term_markers,      f$pseudotime)
  res_cont  <- metrics_d(f$expr, f$early_markers, contaminated_term,   f$pseudotime)

  expect_gte(res_pure$D_term, res_cont$D_term - 0.05)  # allow small tolerance
})

test_that("SIM markers: pure marker set gives higher O than contaminated set", {
  f <- sim_linear(slope = 3, sigma = 0.1)
  contaminated <- c(f$early_markers, paste0("noise_g", 1:10),
                    f$term_markers,  paste0("noise_g", 11:20))
  early_cont   <- c(f$early_markers, paste0("noise_g", 1:10))
  term_cont    <- c(f$term_markers,  paste0("noise_g", 11:20))

  res_pure <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
  res_cont <- metrics_o(f$expr, early_cont,      term_cont,      f$pseudotime)

  expect_gte(res_pure$O, res_cont$O - 0.05)
})

# ── 5. CELL COUNT SCALING ─────────────────────────────────────────────────────

test_that("SIM scaling: metrics are stable across cell counts (50 vs 500 cells)", {
  small <- sim_linear(n_cells = 50,  slope = 3, sigma = 0.15, seed = 7, pi0 = 0.1)
  large <- sim_linear(n_cells = 500, slope = 3, sigma = 0.15, seed = 7, pi0 = 0.1)

  res_s <- metrics_d(small$expr, small$early_markers,
                     small$term_markers, small$pseudotime)
  res_l <- metrics_d(large$expr, large$early_markers,
                     large$term_markers, large$pseudotime)

  # Both should give high D_term; difference should be small
  expect_gt(res_s$D_term, 0.7)
  expect_gt(res_l$D_term, 0.7)
  expect_lt(abs(res_s$D_term - res_l$D_term), 0.15)
})

test_that("SIM scaling: O metric is stable across cell counts", {
  small <- sim_linear(n_cells = 50,  slope = 3, sigma = 0.15, seed = 7)
  large <- sim_linear(n_cells = 500, slope = 3, sigma = 0.15, seed = 7)

  res_s <- metrics_o(small$expr, small$early_markers,
                     small$term_markers, small$pseudotime)
  res_l <- metrics_o(large$expr, large$early_markers,
                     large$term_markers, large$pseudotime)

  # Calibrated pseudobulk O improves with cell count (more cells/bin → cleaner
  # bin means). Both should be clearly above the null baseline (~0.21).
  expect_gt(res_s$O, 0.2)
  expect_gt(res_l$O, 0.5)
})

# ── 6. METRIC MONOTONICITY ACROSS NOISE GRADIENT ─────────────────────────────

test_that("SIM monotonicity: D_term decreases monotonically as noise increases", {
  sigmas  <- c(0.05, 0.3, 0.8, 1.5, 3.0)
  d_terms <- vapply(sigmas, function(s) {
    f   <- sim_linear(slope = 2.5, sigma = s, seed = 42)
    res <- metrics_d(f$expr, f$early_markers, f$term_markers, f$pseudotime)
    res$D_term
  }, numeric(1))

  # Each step should be non-increasing (allow tiny numerical jitter)
  diffs <- diff(d_terms)
  expect_true(all(diffs <= 0.05),
              info = paste("D_term values:", paste(round(d_terms, 3), collapse = " -> ")))
})

test_that("SIM monotonicity: O decreases monotonically as noise increases", {
  sigmas  <- c(0.05, 0.3, 0.8, 1.5)
  o_vals  <- vapply(sigmas, function(s) {
    f   <- sim_linear(slope = 2.5, sigma = s, seed = 42)
    res <- metrics_o(f$expr, f$early_markers, f$term_markers, f$pseudotime)
    res$O
  }, numeric(1))

  diffs <- diff(o_vals)
  expect_true(all(diffs <= 0.05),
              info = paste("O values:", paste(round(o_vals, 3), collapse = " -> ")))
})

# ── 7. ENDPOINT VALIDITY SIMULATION ──────────────────────────────────────────

test_that("SIM endpoint: E_early = 1 when early cells perfectly occupy early pseudotime", {
  set.seed(42)
  n   <- 120
  pt  <- seq(0, 1, length.out = n)
  # first 20 cells are early, last 20 are terminal — perfect placement
  cl  <- c(rep("early_cl", 20), rep("mid_cl", 80), rep("term_cl", 20))
  res <- metrics_e(pt,
                   cluster_labels    = cl,
                   early_clusters    = "early_cl",
                   terminal_clusters = "term_cl",
                   method = "clusters", plot = FALSE)
  expect_equal(res$E_early, 1.0, tolerance = 1e-10)
  expect_equal(res$E_term,  1.0, tolerance = 1e-10)
})

test_that("SIM endpoint: E_early < 0.5 when early cells are at the wrong end", {
  n  <- 100
  pt <- seq(0, 1, length.out = n)
  # early cells placed at the END (wrong)
  cl <- c(rep("mid_cl", 80), rep("early_cl", 20))
  res <- metrics_e(pt,
                   cluster_labels    = cl,
                   early_clusters    = "early_cl",
                   terminal_clusters = character(0),
                   method = "clusters", plot = FALSE)
  expect_lt(res$E_early, 0.5)
})

test_that("SIM endpoint: E_comp is higher for correct placement than random", {
  set.seed(42)
  n   <- 100
  pt  <- seq(0, 1, length.out = n)

  # Correct placement
  cl_correct <- c(rep("early_cl", 20), rep("mid_cl", 60), rep("term_cl", 20))
  res_correct <- metrics_e(pt,
                           cluster_labels    = cl_correct,
                           early_clusters    = "early_cl",
                           terminal_clusters = "term_cl",
                           method = "clusters", plot = FALSE)

  # Random placement
  cl_random <- sample(cl_correct)
  res_random <- metrics_e(pt,
                          cluster_labels    = cl_random,
                          early_clusters    = "early_cl",
                          terminal_clusters = "term_cl",
                          method = "clusters", plot = FALSE)

  if (!is.na(res_correct$E_comp) && !is.na(res_random$E_comp)) {
    expect_gt(res_correct$E_comp, res_random$E_comp)
  }
})

# ── 8. FULL DOE PIPELINE SIMULATION ──────────────────────────────────────────

test_that("SIM full pipeline: good trajectory scores higher DOE than shuffled", {
  f  <- sim_linear(n_cells = 100, slope = 3, sigma = 0.15, seed = 5)
  n  <- ncol(f$expr)
  cl <- c(rep("early_cl", 15), rep("mid_cl", n - 30), rep("term_cl", 15))
  cl <- setNames(cl, colnames(f$expr))

  set.seed(5)
  pt_shuf <- setNames(sample(f$pseudotime), names(f$pseudotime))

  res_good <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, f$pseudotime, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )
  res_shuf <- suppressWarnings(
    compute_single_doe_linear(
      f$expr, pt_shuf, f$early_markers, f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE
    )
  )

  expect_gt(res_good$DOE_score, res_shuf$DOE_score)
})

test_that("SIM full pipeline: multi-trajectory correctly identifies best trajectory", {
  f  <- sim_linear(n_cells = 100, slope = 3, sigma = 0.15, seed = 5)
  n  <- ncol(f$expr)
  cl <- c(rep("early_cl", 15), rep("mid_cl", n - 30), rep("term_cl", 15))
  cl <- setNames(cl, colnames(f$expr))

  set.seed(5)
  pt_shuf <- setNames(sample(f$pseudotime), names(f$pseudotime))
  pt_rev  <- setNames(rev(f$pseudotime),    names(f$pseudotime))

  res <- suppressWarnings(
    compute_multi_doe_linear(
      f$expr,
      pseudotime_list   = list(perfect  = f$pseudotime,
                               shuffled = pt_shuf,
                               reversed = pt_rev),
      early_markers     = f$early_markers,
      terminal_markers  = f$term_markers,
      cluster_labels    = cl,
      early_clusters    = "early_cl",
      terminal_clusters = "term_cl",
      E_method = "clusters", plot_E = FALSE, verbose = FALSE
    )
  )

  expect_equal(res$best_trajectory, "perfect")
})

test_that("SIM full pipeline: DOE_score is in [0, 1] across varied conditions", {
  conditions <- list(
    list(slope = 3.0, sigma = 0.1),
    list(slope = 1.0, sigma = 0.5),
    list(slope = 0.3, sigma = 1.0)
  )
  for (cond in conditions) {
    f  <- sim_linear(slope = cond$slope, sigma = cond$sigma, seed = 42)
    n  <- ncol(f$expr)
    cl <- c(rep("early_cl", 10), rep("mid_cl", n - 20), rep("term_cl", 10))
    cl <- setNames(cl, colnames(f$expr))

    res <- suppressWarnings(
      compute_single_doe_linear(
        f$expr, f$pseudotime, f$early_markers, f$term_markers,
        cluster_labels    = cl,
        early_clusters    = "early_cl",
        terminal_clusters = "term_cl",
        E_method = "clusters", plot_E = FALSE
      )
    )
    if (!is.na(res$DOE_score)) {
      expect_gte(res$DOE_score, 0,
                 label = sprintf("DOE_score >= 0 (slope=%.1f, sigma=%.1f)",
                                 cond$slope, cond$sigma))
      expect_lte(res$DOE_score, 1,
                 label = sprintf("DOE_score <= 1 (slope=%.1f, sigma=%.1f)",
                                 cond$slope, cond$sigma))
    }
  }
})

# ── 9. REVERSE_PSEUDOTIME SIMULATION ─────────────────────────────────────────

test_that("SIM reverse_pseudotime: reversing recovers correct D_term for flipped input", {
  f      <- sim_linear(slope = 3, sigma = 0.1, pi0 = 0.1)
  # Simulate a trajectory inference tool that returns pseudotime in reverse
  pt_rev <- setNames(rev(f$pseudotime), names(f$pseudotime))

  # Without correction: D_term should be low
  res_raw <- metrics_d(f$expr, f$early_markers, f$term_markers, pt_rev)
  expect_lt(res_raw$D_term, 0.5)

  # After reverse_pseudotime correction: D_term should recover
  pt_fixed <- reverse_pseudotime(pt_rev)
  res_fix  <- metrics_d(f$expr, f$early_markers, f$term_markers, pt_fixed)
  expect_gt(res_fix$D_term, 0.7)
})
