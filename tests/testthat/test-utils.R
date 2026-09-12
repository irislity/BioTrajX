# =============================================================================
# test-utils.R  —  Unit tests for reverse_pseudotime()
# =============================================================================
library(testthat)

test_that("reverse_pseudotime: output is in [0, 1]", {
  pt  <- seq(0, 10, length.out = 50)
  rev <- reverse_pseudotime(pt)
  expect_true(all(rev >= 0 & rev <= 1))
})

test_that("reverse_pseudotime: minimum maps to 1, maximum maps to 0", {
  pt  <- c(2, 5, 8, 10)
  rev <- reverse_pseudotime(pt)
  expect_equal(rev[which.min(pt)], 1)
  expect_equal(rev[which.max(pt)], 0)
})

test_that("reverse_pseudotime: is a true reversal (rank order flipped)", {
  pt  <- c(1, 3, 5, 7, 9)
  rev <- reverse_pseudotime(pt)
  expect_equal(rank(rev), rev(rank(pt)))
})

test_that("reverse_pseudotime: preserves vector length", {
  pt <- runif(100)
  expect_length(reverse_pseudotime(pt), 100)
})

test_that("reverse_pseudotime: single unique value returns zeros with warning", {
  expect_warning(
    out <- reverse_pseudotime(rep(5, 10)),
    "identical"
  )
  expect_true(all(out == 0))
})

test_that("reverse_pseudotime: all-NA input returns all-NA", {
  out <- reverse_pseudotime(rep(NA_real_, 5))
  expect_true(all(is.na(out)))
})

test_that("reverse_pseudotime: non-numeric input throws error", {
  expect_error(reverse_pseudotime(c("a", "b")), "numeric")
})

test_that("reverse_pseudotime: already-scaled [0,1] input stays in [0,1]", {
  pt  <- seq(0, 1, length.out = 20)
  rev <- reverse_pseudotime(pt)
  expect_equal(rev[1],  1)
  expect_equal(rev[20], 0)
})

test_that("reverse_pseudotime: applying twice recovers original (up to scaling)", {
  pt      <- c(1, 4, 7, 10)
  once    <- reverse_pseudotime(pt)
  twice   <- reverse_pseudotime(once)
  # twice should be monotonically increasing (same rank order as original)
  expect_equal(rank(twice), rank(pt))
})

# ── Additional edge cases ─────────────────────────────────────────────────────

test_that("reverse_pseudotime: single-element vector returns 0", {
  expect_warning(out <- reverse_pseudotime(5), "identical")
  expect_equal(out, 0)
})

test_that("reverse_pseudotime: negative pseudotime works correctly", {
  pt  <- c(-10, -5, 0, 5, 10)
  rev <- reverse_pseudotime(pt)
  expect_equal(rev[1], 1)   # minimum maps to 1
  expect_equal(rev[5], 0)   # maximum maps to 0
  expect_true(all(rev >= 0 & rev <= 1))
})

test_that("reverse_pseudotime: named vector preserves names", {
  pt  <- c(a = 1, b = 2, c = 3)
  rev <- reverse_pseudotime(pt)
  expect_equal(names(rev), c("a", "b", "c"))
})

test_that("reverse_pseudotime: Inf input throws error or returns NA", {
  pt <- c(1, 2, Inf)
  # Either errors or produces non-finite output — should not silently return [0,1]
  out <- reverse_pseudotime(pt)
  expect_true(any(!is.finite(out)))
})
