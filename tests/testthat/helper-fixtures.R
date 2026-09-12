# =============================================================================
# helper-fixtures.R
# Shared synthetic data for BioTrajX tests.
# Loaded automatically by testthat before every test file.
#
# Expression matrices are simulated with Zero-Inflated Negative Binomial (ZINB)
# noise to mimic the sparsity and overdispersion of real scRNA-seq count data.
# =============================================================================

# ── ZINB sampler ─────────────────────────────────────────────────────────────
# Draws ZINB counts and returns log1p-normalised values, matching the
# log-normalised input that real scRNA-seq users supply to BioTrajX.
# size : NB dispersion (lower = more overdispersed, typical scRNA-seq: 0.3-1)
# pi0  : zero-inflation probability (dropout rate, typical scRNA-seq: ~0.5)
rzinb <- function(n, mu, size = 0.5, pi0 = 0.5) {
  counts <- rnbinom(n, mu = pmax(mu, 1e-6), size = size)
  counts[runif(n) < pi0] <- 0L
  log1p(counts)   # log-normalise as in real scRNA-seq workflows
}

# ── ZINB parameter defaults ───────────────────────────────────────────────────
# Centralise the ZINB knobs so every fixture accepts them consistently.
# slope : magnitude of the linear mean trend (early: (1+slope) → 1, terminal: 1 → (1+slope))
# size  : NB dispersion parameter (lower = more overdispersed)
# pi0   : zero-inflation / dropout probability

# ── Reproducible RNG ─────────────────────────────────────────────────────────
set.seed(42)

# ── Small deterministic expression matrix (genes x cells) ────────────────────
#   50 genes, 80 cells (default)
#   early genes (rows 1-10)  : mean DECREASES along pseudotime (1+slope → 1)
#   terminal genes (rows 11-20): mean INCREASES along pseudotime (1 → 1+slope)
#   noise genes (rows 21-end): constant mean (2), no trend
#
#   slope : linear gradient magnitude (default 7 → mean range 1–8)
#   size  : NB dispersion (default 0.5, typical scRNA-seq)
#   pi0   : dropout probability (default 0.5)
make_linear_fixture <- function(n_cells = 80, n_genes = 50, seed = 42,
                                slope = 7, size = 1.5, pi0 = 0.05) {
  set.seed(seed)
  pt <- seq(0, 1, length.out = n_cells)

  # Background: constant-mean ZINB noise
  mat <- matrix(rzinb(n_genes * n_cells, mu = 2, size = size, pi0 = pi0),
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_cells))

  # early genes: mean decreases (1+slope) → 1
  for (i in 1:10)
    mat[i, ] <- rzinb(n_cells, mu = (1 + slope) - slope * pt, size = size, pi0 = pi0)

  # terminal genes: mean increases 1 → (1+slope)
  for (i in 11:20)
    mat[i, ] <- rzinb(n_cells, mu = 1 + slope * pt, size = size, pi0 = pi0)

  list(
    expr          = mat,
    pseudotime    = setNames(pt, colnames(mat)),
    early_markers = paste0("gene", 1:10),
    term_markers  = paste0("gene", 11:20),
    n_cells       = n_cells,
    n_genes       = n_genes,
    slope         = slope,
    size          = size,
    pi0           = pi0
  )
}

# ── Reversed / anti-correlated fixture (D should be low) ─────────────────────
make_reversed_fixture <- function(n_cells = 80, seed = 42) {
  f <- make_linear_fixture(n_cells = n_cells, seed = seed)
  f$pseudotime <- setNames(rev(f$pseudotime), names(f$pseudotime))
  f
}

# ── Shuffled pseudotime fixture (negative control) ────────────────────────────
#   Same ZINB linear expression data, but pseudotime is randomly permuted
#   across cells, destroying any real ordering signal.
#   Expected: D, O, E all collapse toward the noise floor (~0).
make_shuffled_fixture <- function(n_cells = 80, seed = 42) {
  f <- make_linear_fixture(n_cells = n_cells, seed = seed)
  set.seed(seed)
  f$pseudotime <- setNames(sample(f$pseudotime), names(f$pseudotime))
  f
}

# ── Flat / uninformative fixture (all metrics near 0.5) ──────────────────────
#   Constant mean → no pseudotime trend
make_flat_fixture <- function(n_cells = 80, n_genes = 50, seed = 42) {
  set.seed(seed)
  pt  <- seq(0, 1, length.out = n_cells)
  mat <- matrix(rzinb(n_genes * n_cells, mu = 2, pi0 = 0.05),
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_cells))
  list(
    expr          = mat,
    pseudotime    = setNames(pt, colnames(mat)),
    early_markers = paste0("gene", 1:10),
    term_markers  = paste0("gene", 11:20)
  )
}

# ── Branched fixture ──────────────────────────────────────────────────────────
#   A → B and A → C topology:
#     A cells (early progenitors): pt 0 → 0.5, genes 1-10 early (high→low)
#     B cells (terminal fate B):   pt 0.5 → 1,  genes 11-20 terminal (low→high)
#     C cells (terminal fate C):   pt 0.5 → 1,  genes 21-30 terminal (low→high)
#     genes 31-40: background noise
#   DOE is computed separately for the AB and AC trajectories.
make_branched_fixture <- function(n_cells_per_branch = 50, seed = 42) {
  set.seed(seed)
  n_total <- n_cells_per_branch * 3
  n_genes <- 40

  idx_a <- seq_len(n_cells_per_branch)
  idx_b <- seq(n_cells_per_branch + 1L,      2L * n_cells_per_branch)
  idx_c <- seq(2L * n_cells_per_branch + 1L, 3L * n_cells_per_branch)

  pt_a <- seq(0,   0.5, length.out = n_cells_per_branch)
  pt_b <- seq(0.5, 1.0, length.out = n_cells_per_branch)
  pt_c <- seq(0.5, 1.0, length.out = n_cells_per_branch)
  pt   <- c(pt_a, pt_b, pt_c)

  # ZINB parameters matched to make_linear_fixture (size=1.5, pi0=0.05) so that
  # Branched ZINB is a fair branched analogue of ZINB Linear, not an accidental
  # extreme-dropout scenario.
  size <- 1.5; pi0 <- 0.05

  # Background noise
  mat <- matrix(rzinb(n_genes * n_total, mu = 2, size = size, pi0 = pi0),
                nrow = n_genes, ncol = n_total)
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_total))

  # Shared early markers (genes 1-10): decrease along global pt for ALL cells
  for (i in 1:10)
    mat[i, ] <- rzinb(n_total, mu = 8 - 6 * pt, size = size, pi0 = pi0)

  # B-specific terminal markers (genes 11-20): increase only in B cells
  # (A cells remain at background — D_term degradation reflects branching, not noise)
  pt_b_rsc <- (pt_b - 0.5) / 0.5   # rescale to 0→1 within branch B
  for (i in 11:20)
    mat[i, idx_b] <- rzinb(n_cells_per_branch, mu = 2 + 6 * pt_b_rsc, size = size, pi0 = pi0)

  # C-specific terminal markers (genes 21-30): increase only in C cells
  pt_c_rsc <- (pt_c - 0.5) / 0.5   # rescale to 0→1 within branch C
  for (i in 21:30)
    mat[i, idx_c] <- rzinb(n_cells_per_branch, mu = 2 + 6 * pt_c_rsc, size = size, pi0 = pi0)

  cluster_labels <- c(rep("A", n_cells_per_branch),
                       rep("B", n_cells_per_branch),
                       rep("C", n_cells_per_branch))
  names(cluster_labels) <- colnames(mat)

  list(
    expr               = mat,
    pseudotime         = setNames(pt, colnames(mat)),
    cluster_labels     = cluster_labels,
    early_markers_list = list(AB = paste0("gene", 1:10),  AC = paste0("gene", 1:10)),
    term_markers_list  = list(AB = paste0("gene", 11:20), AC = paste0("gene", 21:30)),
    branch_filters     = list(
      AB = list(include = c("A", "B")),
      AC = list(include = c("A", "C"))
    )
  )
}

# =============================================================================
# CEILING / CORRECTNESS FIXTURES
# These use Gaussian noise instead of ZINB and are designed so that each metric
# has a known, assertable expected value (D ≈ 1, O ≈ 1, E ≈ 1).
# Use them to verify that metric implementations are mathematically correct.
# =============================================================================

# ── Gaussian linear fixture ───────────────────────────────────────────────────
#   Replaces ZINB with additive Gaussian noise so the full score range is
#   preserved after log-normalisation.  With slope = 7 and sigma = 0.3 the
#   signal-to-noise ratio is high enough that D ≈ 1, O ≈ 1, and the per-cell
#   marker scores spread evenly over [1, 8], giving GMM a clean split → E ≈ 1.
#
#   slope : mean range of marker genes along pseudotime (default 7 → 1 to 8)
#   sigma : Gaussian noise std dev (default 0.3)
make_gaussian_fixture <- function(n_cells = 100, n_genes = 50, seed = 42,
                                  slope = 7, sigma = 0.3) {
  set.seed(seed)
  pt <- seq(0, 1, length.out = n_cells)

  # Background noise genes: constant mean, no trend
  mat <- matrix(rnorm(n_genes * n_cells, mean = 2, sd = sigma),
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_cells))

  # early genes: mean decreases (1 + slope) → 1
  for (i in 1:10)
    mat[i, ] <- (1 + slope) - slope * pt + rnorm(n_cells, 0, sigma)

  # terminal genes: mean increases 1 → (1 + slope)
  for (i in 11:20)
    mat[i, ] <- 1 + slope * pt + rnorm(n_cells, 0, sigma)

  mat[mat < 0] <- 0

  list(
    expr          = mat,
    pseudotime    = setNames(pt, colnames(mat)),
    early_markers = paste0("gene", 1:10),
    term_markers  = paste0("gene", 11:20),
    n_cells       = n_cells,
    n_genes       = n_genes,
    slope         = slope,
    sigma         = sigma
  )
}

# ── Gaussian reversed fixture ─────────────────────────────────────────────────
#   Same Gaussian data, pseudotime flipped.
#   Expected: D_early ≈ 0, D_term ≈ 0, O ≈ 0 (floor reference).
make_gaussian_reversed_fixture <- function(n_cells = 100, seed = 42) {
  f <- make_gaussian_fixture(n_cells = n_cells, seed = seed)
  f$pseudotime <- setNames(rev(f$pseudotime), names(f$pseudotime))
  f
}

# =============================================================================
# ENDPOINT FIXTURE  (E-metric validation)
# Tests whether E correctly detects that known endpoint cells sit at the
# right ends of pseudotime.  Uses a step-function model: the first early_frac
# of cells are in a pure early state (high early, low terminal expression) and
# the last term_frac are in a pure terminal state, with middle cells at
# background level.  Provides ground-truth cell sets for direct assertion.
# =============================================================================

# ── Endpoint fixture ──────────────────────────────────────────────────────────
#   early_frac : fraction of cells at the early endpoint  (default 0.25)
#   term_frac  : fraction of cells at the terminal endpoint (default 0.25)
#   high / low : mean counts for endpoint vs. background expression
#   size / pi0 : ZINB parameters — less harsh defaults than the linear fixture
#                (size = 1, pi0 = 0.2) to keep endpoint scores well-separated
make_endpoint_fixture <- function(n_cells = 80, n_genes = 50, seed = 42,
                                  early_frac = 0.25, term_frac = 0.25,
                                  high = 8, low = 1,
                                  size = 1, pi0 = 0.2) {
  set.seed(seed)
  pt <- seq(0, 1, length.out = n_cells)

  # Background expression for all genes
  mat <- matrix(rzinb(n_genes * n_cells, mu = 2, size = size, pi0 = pi0),
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_cells))

  # Step-function means: endpoint cells get `high`, all others get `low`
  early_mu <- ifelse(pt <= early_frac,      high, low)
  term_mu  <- ifelse(pt >= 1 - term_frac,   high, low)

  for (i in 1:10)
    mat[i, ] <- rzinb(n_cells, mu = early_mu, size = size, pi0 = pi0)
  for (i in 11:20)
    mat[i, ] <- rzinb(n_cells, mu = term_mu,  size = size, pi0 = pi0)

  list(
    expr                   = mat,
    pseudotime             = setNames(pt, colnames(mat)),
    early_markers          = paste0("gene", 1:10),
    term_markers           = paste0("gene", 11:20),
    # Ground-truth endpoint cell sets for asserting GMM output directly
    true_early_endpoint    = colnames(mat)[pt <= early_frac],
    true_terminal_endpoint = colnames(mat)[pt >= 1 - term_frac],
    n_cells                = n_cells,
    n_genes                = n_genes,
    early_frac             = early_frac,
    term_frac              = term_frac,
    high                   = high,
    low                    = low,
    size                   = size,
    pi0                    = pi0
  )
}

# =============================================================================
# STRESS / ROBUSTNESS FIXTURES  (ZINB degradation tests)
# These verify that metrics degrade gracefully — assertions are directional
# (expect_gt / expect_lt), not exact.
# =============================================================================

# ── High-dropout fixture (pi0 = 0.9) ─────────────────────────────────────────
#   Mimics aggressive dropout: ~90 % of counts are zeroed.
#   D / O / E scores should degrade relative to the standard linear fixture.
make_high_dropout_fixture <- function(n_cells = 80, n_genes = 50, seed = 42) {
  make_linear_fixture(n_cells = n_cells, n_genes = n_genes, seed = seed,
                      slope = 7, size = 0.5, pi0 = 0.9)
}

# ── Low-signal fixture (slope = 0.5) ─────────────────────────────────────────
#   The mean trajectory is nearly flat (range 1 → 1.5).
#   Weak gradient tests sensitivity of metrics near the noise floor.
make_low_signal_fixture <- function(n_cells = 80, n_genes = 50, seed = 42) {
  make_linear_fixture(n_cells = n_cells, n_genes = n_genes, seed = seed,
                      slope = 0.5, size = 0.5, pi0 = 0.5)
}

# ── Small-dataset fixture (n_cells = 20, n_genes = 30) ───────────────────────
#   Tests metric behaviour under minimal sample sizes.
#   GMM may fall back to quantile/kmeans splitting; isotonic bins are reduced.
make_small_fixture <- function(seed = 42) {
  make_linear_fixture(n_cells = 20, n_genes = 30, seed = seed,
                      slope = 7, size = 0.5, pi0 = 0.5)
}

# ── High-overdispersion fixture (size = 0.2) ─────────────────────────────────
#   Very bursty counts (NB size = 0.2) create high cell-to-cell variance.
#   Tests robustness of metrics to extreme count noise.
make_high_overdispersion_fixture <- function(n_cells = 80, n_genes = 50, seed = 42) {
  make_linear_fixture(n_cells = n_cells, n_genes = n_genes, seed = seed,
                      slope = 7, size = 0.2, pi0 = 0.5)
}

# ── Stepwise-repression fixture (03 scenario — expected endpoint) ─────────────
#   Early genes are HIGH at pseudotime=0 and show a sharp step-function drop
#   (repression) partway along the trajectory.  Terminal genes show the reverse
#   (step-function activation at the terminal end).  Ground-truth endpoint cell
#   sets are retained for direct assertion.
#   Expected: E_early ≈ 1, E_term ≈ 1, D ≈ 0.8, O ≈ 0.85
make_stepwise_repression_fixture <- function(n_cells = 80, n_genes = 50, seed = 42,
                                             early_frac = 0.25, term_frac = 0.25,
                                             high = 8, low = 1, size = 1, pi0 = 0.2) {
  make_endpoint_fixture(n_cells = n_cells, n_genes = n_genes, seed = seed,
                        early_frac = early_frac, term_frac = term_frac,
                        high = high, low = low, size = size, pi0 = pi0)
}

# ── Stepwise-activation fixture (04 scenario — reverse endpoint) ──────────────
#   Same step-function data as make_stepwise_repression_fixture() but with
#   pseudotime reversed.  Now terminal cells appear at pseudotime=0 and early
#   cells at pseudotime=1 — the opposite of the correct orientation.
#   Expected: E_term ≈ 0, E_early ≈ 0 (endpoints at wrong ends of pseudotime).
make_stepwise_activation_fixture <- function(n_cells = 80, n_genes = 50, seed = 42,
                                             early_frac = 0.25, term_frac = 0.25,
                                             high = 8, low = 1, size = 1, pi0 = 0.2) {
  f <- make_endpoint_fixture(n_cells = n_cells, n_genes = n_genes, seed = seed,
                             early_frac = early_frac, term_frac = term_frac,
                             high = high, low = low, size = size, pi0 = pi0)
  f$pseudotime <- setNames(rev(f$pseudotime), names(f$pseudotime))
  f
}

# ── Gaussian branched fixture (10 scenario — ceiling branched) ────────────────
#   A → B and A → C topology with Gaussian noise (high SNR ceiling reference).
#     A cells (early progenitors): pt 0 → 0.5, genes 1-10 early (high→low)
#     B cells (terminal fate B):   pt 0.5 → 1,  genes 11-20 terminal (low→high)
#     C cells (terminal fate C):   pt 0.5 → 1,  genes 21-30 terminal (low→high)
#     genes 31-40: background noise
#   Expected per-trajectory: D≈1, O≈1, E≈1
make_gaussian_branched_fixture <- function(n_cells_per_branch = 50, seed = 42,
                                           slope = 7, sigma = 0.3) {
  set.seed(seed)
  n_total <- n_cells_per_branch * 3
  n_genes <- 40

  idx_a <- seq_len(n_cells_per_branch)
  idx_b <- seq(n_cells_per_branch + 1L,      2L * n_cells_per_branch)
  idx_c <- seq(2L * n_cells_per_branch + 1L, 3L * n_cells_per_branch)

  pt_a <- seq(0,   0.5, length.out = n_cells_per_branch)
  pt_b <- seq(0.5, 1.0, length.out = n_cells_per_branch)
  pt_c <- seq(0.5, 1.0, length.out = n_cells_per_branch)
  pt   <- c(pt_a, pt_b, pt_c)

  mat <- matrix(rnorm(n_genes * n_total, mean = 2, sd = sigma),
                nrow = n_genes, ncol = n_total)
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_total))

  # Shared early markers (genes 1-10): decrease along global pt for ALL cells
  for (i in 1:10)
    mat[i, ] <- (1 + slope) - slope * pt + rnorm(n_total, 0, sigma)

  # B-specific terminal markers (genes 11-20): increase monotonically across
  # the full AB trajectory (A cells pt 0→0.5, B cells pt 0.5→1.0) using global
  # pseudotime so that rho_term ≈ 1 on the AB subset.
  for (i in 11:20) {
    mat[i, idx_a] <- 1 + slope * pt_a + rnorm(n_cells_per_branch, 0, sigma)
    mat[i, idx_b] <- 1 + slope * pt_b + rnorm(n_cells_per_branch, 0, sigma)
    # C cells remain at background
  }

  # C-specific terminal markers (genes 21-30): increase monotonically across
  # the full AC trajectory (A cells pt 0→0.5, C cells pt 0.5→1.0) using global
  # pseudotime so that rho_term ≈ 1 on the AC subset.
  for (i in 21:30) {
    mat[i, idx_a] <- 1 + slope * pt_a + rnorm(n_cells_per_branch, 0, sigma)
    mat[i, idx_c] <- 1 + slope * pt_c + rnorm(n_cells_per_branch, 0, sigma)
    # B cells remain at background
  }

  mat[mat < 0] <- 0

  cluster_labels <- c(rep("A", n_cells_per_branch),
                       rep("B", n_cells_per_branch),
                       rep("C", n_cells_per_branch))
  names(cluster_labels) <- colnames(mat)

  list(
    expr               = mat,
    pseudotime         = setNames(pt, colnames(mat)),
    cluster_labels     = cluster_labels,
    early_markers_list = list(AB = paste0("gene", 1:10),  AC = paste0("gene", 1:10)),
    term_markers_list  = list(AB = paste0("gene", 11:20), AC = paste0("gene", 21:30)),
    branch_filters     = list(
      AB = list(include = c("A", "B")),
      AC = list(include = c("A", "C"))
    ),
    slope = slope,
    sigma = sigma
  )
}

# ── Gaussian branched reversed fixture (11 scenario) ─────────────────────────
#   Pseudotime flipped relative to the Gaussian branched fixture.
#   Expected per-branch: D≈0, O unchanged (orientation-invariant), E≈0.
make_gaussian_branched_reversed_fixture <- function(n_cells_per_branch = 50, seed = 42) {
  f <- make_gaussian_branched_fixture(n_cells_per_branch = n_cells_per_branch, seed = seed)
  f$pseudotime <- setNames(rev(f$pseudotime), names(f$pseudotime))
  f
}

# ── Branched ZINB fixture (12 scenario) ───────────────────────────────────────
#   Two-branch trajectory with ZINB noise (same as make_branched_fixture but
#   explicitly named for the manuscript scenario list).
#   Expected per-branch: D and E scores reduced relative to Gaussian baseline.
make_branched_zinb_fixture <- function(n_cells_per_branch = 50, seed = 42) {
  make_branched_fixture(n_cells_per_branch = n_cells_per_branch, seed = seed)
}

# ── Branched shuffled pseudotime fixture (negative control) ──────────────────
#   Same branched ZINB fixture, pseudotime randomly permuted across ALL cells
#   (destroys both within-branch and across-branch ordering).
#   Expected per-branch: D, O, E all collapse toward the noise floor (~0).
make_branched_shuffled_fixture <- function(n_cells_per_branch = 50, seed = 42) {
  f <- make_branched_zinb_fixture(n_cells_per_branch = n_cells_per_branch, seed = seed)
  set.seed(seed)
  f$pseudotime <- setNames(sample(f$pseudotime), names(f$pseudotime))
  f
}

# ── Branched low-signal fixture (13 scenario) ─────────────────────────────────
#   A → B and A → C topology with weak gradient (slope = 0.5) and ZINB noise.
#   Tests metric sensitivity near the noise floor in branched settings.
#   Expected per-trajectory: all metrics reduced vs. branched ZINB baseline.
make_branched_low_signal_fixture <- function(n_cells_per_branch = 50, seed = 42,
                                             slope = 0.5) {
  set.seed(seed)
  n_total <- n_cells_per_branch * 3
  n_genes <- 40

  idx_a <- seq_len(n_cells_per_branch)
  idx_b <- seq(n_cells_per_branch + 1L,      2L * n_cells_per_branch)
  idx_c <- seq(2L * n_cells_per_branch + 1L, 3L * n_cells_per_branch)

  pt_a <- seq(0,   0.5, length.out = n_cells_per_branch)
  pt_b <- seq(0.5, 1.0, length.out = n_cells_per_branch)
  pt_c <- seq(0.5, 1.0, length.out = n_cells_per_branch)
  pt   <- c(pt_a, pt_b, pt_c)

  mat <- matrix(rzinb(n_genes * n_total, mu = 2),
                nrow = n_genes, ncol = n_total)
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_total))

  # Shared early markers (genes 1-10): weak decrease along global pt
  for (i in 1:10)
    mat[i, ] <- rzinb(n_total, mu = pmax((1 + slope) - slope * pt, 1e-6))

  # B-specific terminal markers (genes 11-20): weak increase in B cells only
  pt_b_rsc <- (pt_b - 0.5) / 0.5
  for (i in 11:20)
    mat[i, idx_b] <- rzinb(n_cells_per_branch, mu = 1 + slope * pt_b_rsc)

  # C-specific terminal markers (genes 21-30): weak increase in C cells only
  pt_c_rsc <- (pt_c - 0.5) / 0.5
  for (i in 21:30)
    mat[i, idx_c] <- rzinb(n_cells_per_branch, mu = 1 + slope * pt_c_rsc)

  cluster_labels <- c(rep("A", n_cells_per_branch),
                       rep("B", n_cells_per_branch),
                       rep("C", n_cells_per_branch))
  names(cluster_labels) <- colnames(mat)

  list(
    expr               = mat,
    pseudotime         = setNames(pt, colnames(mat)),
    cluster_labels     = cluster_labels,
    early_markers_list = list(AB = paste0("gene", 1:10),  AC = paste0("gene", 1:10)),
    term_markers_list  = list(AB = paste0("gene", 11:20), AC = paste0("gene", 21:30)),
    branch_filters     = list(
      AB = list(include = c("A", "B")),
      AC = list(include = c("A", "C"))
    ),
    slope = slope
  )
}

# ── Branched stepwise repression fixture ──────────────────────────────────────
#   A → B and A → C topology with step-function markers, analogous to the
#   linear make_stepwise_repression_fixture but for branched trajectories.
#     Early markers (genes 1-5)     : HIGH at start of A (pt ≤ early_frac)
#     B-terminal markers (genes 6-10) : HIGH at end of B  (pt ≥ 1 - term_frac)
#     C-terminal markers (genes 11-15): HIGH at end of C  (pt ≥ 1 - term_frac)
#   Provides per-branch ground-truth endpoint cell sets for direct assertion.
#   Expected per branch: E_early ≈ 1, E_term ≈ 1.
make_branched_stepwise_fixture <- function(n_cells_per_branch = 50, seed = 42,
                                           early_frac = 0.25, term_frac = 0.25,
                                           high = 8, low = 1, size = 1, pi0 = 0.2) {
  set.seed(seed)
  n_total <- n_cells_per_branch * 3
  n_genes <- 40

  idx_a <- seq_len(n_cells_per_branch)
  idx_b <- seq(n_cells_per_branch + 1L,      2L * n_cells_per_branch)
  idx_c <- seq(2L * n_cells_per_branch + 1L, 3L * n_cells_per_branch)

  pt_a <- seq(0,   0.5, length.out = n_cells_per_branch)
  pt_b <- seq(0.5, 1.0, length.out = n_cells_per_branch)
  pt_c <- seq(0.5, 1.0, length.out = n_cells_per_branch)
  pt   <- c(pt_a, pt_b, pt_c)

  # Background expression
  mat <- matrix(rzinb(n_genes * n_total, mu = low, size = size, pi0 = pi0),
                nrow = n_genes, ncol = n_total)
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_total))

  # Early markers (genes 1-10): HIGH for early A cells (global pt ≤ early_frac)
  early_mu <- ifelse(pt <= early_frac, high, low)
  for (i in 1:10)
    mat[i, ] <- rzinb(n_total, mu = early_mu, size = size, pi0 = pi0)

  # B-terminal markers (genes 11-20): HIGH for late B cells (pt ≥ 1 - term_frac)
  b_term_mu <- ifelse(pt_b >= (1 - term_frac), high, low)
  for (i in 11:20)
    mat[i, idx_b] <- rzinb(n_cells_per_branch, mu = b_term_mu, size = size, pi0 = pi0)

  # C-terminal markers (genes 21-30): HIGH for late C cells (pt ≥ 1 - term_frac)
  c_term_mu <- ifelse(pt_c >= (1 - term_frac), high, low)
  for (i in 21:30)
    mat[i, idx_c] <- rzinb(n_cells_per_branch, mu = c_term_mu, size = size, pi0 = pi0)

  cluster_labels <- c(rep("A", n_cells_per_branch),
                      rep("B", n_cells_per_branch),
                      rep("C", n_cells_per_branch))
  names(cluster_labels) <- colnames(mat)
  pt_named <- setNames(pt, colnames(mat))

  list(
    expr               = mat,
    pseudotime         = pt_named,
    cluster_labels     = cluster_labels,
    early_markers_list = list(AB = paste0("gene", 1:10),  AC = paste0("gene", 1:10)),
    term_markers_list  = list(AB = paste0("gene", 11:20), AC = paste0("gene", 21:30)),
    branch_filters     = list(
      AB = list(include = c("A", "B")),
      AC = list(include = c("A", "C"))
    ),
    # Ground-truth endpoint sets (identified by expression pattern, not pseudotime,
    # so they remain valid even after pseudotime reversal)
    true_early_endpoint      = names(pt_named)[pt_named <= early_frac],
    true_terminal_endpoint_B = names(pt_named)[pt_named >= (1 - term_frac) & cluster_labels == "B"],
    true_terminal_endpoint_C = names(pt_named)[pt_named >= (1 - term_frac) & cluster_labels == "C"],
    early_frac = early_frac, term_frac = term_frac,
    high = high, low = low, size = size, pi0 = pi0
  )
}

# ── Branched stepwise reversed fixture ────────────────────────────────────────
#   Same step-function data as make_branched_stepwise_fixture but pseudotime
#   reversed, so true endpoints appear at the wrong end of the trajectory.
#   Expected per branch: E_early ≈ 0, E_term ≈ 0.
make_branched_stepwise_reversed_fixture <- function(n_cells_per_branch = 50, seed = 42,
                                                    early_frac = 0.25, term_frac = 0.25,
                                                    high = 8, low = 1, size = 1, pi0 = 0.2) {
  f <- make_branched_stepwise_fixture(
    n_cells_per_branch = n_cells_per_branch, seed = seed,
    early_frac = early_frac, term_frac = term_frac,
    high = high, low = low, size = size, pi0 = pi0
  )
  f$pseudotime <- setNames(rev(f$pseudotime), names(f$pseudotime))
  f
}
