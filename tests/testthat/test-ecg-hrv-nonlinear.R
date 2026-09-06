library(testthat)
library(PhysioECG)

# --- Poincare tests ---

test_that("ecgHRVpoincare returns SD1 near 0 for constant RR intervals", {
  n <- 50
  rr <- data.frame(
    channel = rep(1L, n),
    rr_ms = rep(800, n),
    time_sec = cumsum(rep(0.8, n))
  )

  result <- ecgHRVpoincare(rr)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(all(c("channel", "sd1", "sd2", "sd1_sd2_ratio") %in% names(result)))
  expect_equal(result$sd1, 0, tolerance = 1e-10)
  expect_equal(result$sd2, 0, tolerance = 1e-10)
})

test_that("ecgHRVpoincare returns expected values for alternating RR", {
  # Alternating pattern: 800, 900, 800, 900, ...
  n <- 100
  rr_vals <- rep(c(800, 900), length.out = n)
  rr <- data.frame(
    channel = rep(1L, n),
    rr_ms = rr_vals,
    time_sec = cumsum(rr_vals / 1000)
  )

  result <- ecgHRVpoincare(rr)

  # SD1 should reflect successive differences (all diffs are +/-100)
  expected_sd1 <- sqrt(0.5) * sd(diff(rr_vals))
  expect_equal(result$sd1, expected_sd1, tolerance = 0.01)

  # SD2 should be computed from overall variance and successive differences
  expected_sd2_sq <- 2 * var(rr_vals) - 0.5 * var(diff(rr_vals))
  if (expected_sd2_sq > 0) {
    expect_equal(result$sd2, sqrt(expected_sd2_sq), tolerance = 0.01)
  }

  # For pure alternating pattern, SD2 is ~0 so ratio is NA
  expect_true(is.na(result$sd1_sd2_ratio))
})

test_that("ecgHRVpoincare handles multiple channels", {
  rr <- data.frame(
    channel = c(rep(1L, 30), rep(2L, 30)),
    rr_ms = c(rep(800, 30), rep(c(700, 900), 15)),
    time_sec = c(cumsum(rep(0.8, 30)), cumsum(rep(0.8, 30)))
  )

  result <- ecgHRVpoincare(rr)
  expect_equal(nrow(result), 2)
  expect_equal(result$channel, c(1L, 2L))
  # Channel 1 (constant) should have SD1=0, channel 2 (alternating) should have SD1>0
  expect_equal(result$sd1[1], 0, tolerance = 1e-10)
  expect_true(result$sd1[2] > 0)
})

# --- Sample Entropy tests ---

test_that("ecgSampleEntropy: periodic signal has lower entropy than random", {
  set.seed(42)
  n <- 200

  # Periodic: repeating pattern
  periodic_rr <- rep(c(800, 850, 900, 850), length.out = n)
  rr_periodic <- data.frame(
    channel = rep(1L, n),
    rr_ms = periodic_rr,
    time_sec = cumsum(periodic_rr / 1000)
  )

  # Random: uniform noise
  random_rr <- runif(n, min = 700, max = 1000)
  rr_random <- data.frame(
    channel = rep(1L, n),
    rr_ms = random_rr,
    time_sec = cumsum(random_rr / 1000)
  )

  result_periodic <- ecgSampleEntropy(rr_periodic)
  result_random <- ecgSampleEntropy(rr_random)

  expect_s3_class(result_periodic, "data.frame")
  expect_true(all(c("channel", "sample_entropy", "m", "r") %in% names(result_periodic)))

  # Both should produce finite values
  expect_true(is.finite(result_periodic$sample_entropy))
  expect_true(is.finite(result_random$sample_entropy))

  # Periodic should have lower entropy than random
  expect_true(result_periodic$sample_entropy < result_random$sample_entropy)
})

test_that("ecgSampleEntropy returns correct m and r values", {
  n <- 100
  rr <- data.frame(
    channel = rep(1L, n),
    rr_ms = rnorm(n, mean = 850, sd = 30),
    time_sec = cumsum(rep(0.85, n))
  )

  result <- ecgSampleEntropy(rr, m = 3, r_factor = 0.15)
  expect_equal(result$m, 3L)
  expect_equal(result$r, 0.15 * sd(rr$rr_ms), tolerance = 0.001)
})

test_that("ecgSampleEntropy returns NA for constant input", {
  n <- 50
  rr <- data.frame(
    channel = rep(1L, n),
    rr_ms = rep(800, n),
    time_sec = cumsum(rep(0.8, n))
  )

  result <- ecgSampleEntropy(rr)
  expect_true(is.na(result$sample_entropy))
})

# --- DFA tests ---

test_that("ecgDFA returns valid structure with sufficient data", {
  set.seed(123)
  n <- 300
  rr_ms <- cumsum(rnorm(n, mean = 0, sd = 10)) + 850
  rr <- data.frame(
    channel = rep(1L, n),
    rr_ms = rr_ms,
    time_sec = cumsum(rr_ms / 1000)
  )

  result <- ecgDFA(rr)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(all(c("channel", "alpha1", "alpha2") %in% names(result)))
  expect_true(is.finite(result$alpha1))
})

test_that("ecgDFA alpha1 is around 1.5 for integrated white noise (Brownian)", {
  set.seed(99)
  # Brownian motion (cumulated white noise) has DFA alpha ~ 1.5
  n <- 500
  rr_ms <- cumsum(rnorm(n)) + 850
  rr <- data.frame(
    channel = rep(1L, n),
    rr_ms = rr_ms,
    time_sec = cumsum(abs(rr_ms) / 1000)
  )

  result <- ecgDFA(rr, short_range = c(4, 16), long_range = c(16, 64))
  # Alpha for Brownian motion should be approximately 1.5
  # Use wider tolerance since this is stochastic
  expect_true(result$alpha1 > 0.8 && result$alpha1 < 2.2)
})

test_that("ecgDFA handles multi-channel input", {
  set.seed(7)
  n <- 200
  rr <- data.frame(
    channel = c(rep(1L, n), rep(2L, n)),
    rr_ms = c(rnorm(n, 850, 30), rnorm(n, 900, 40)),
    time_sec = c(cumsum(rep(0.85, n)), cumsum(rep(0.9, n)))
  )

  result <- ecgDFA(rr)
  expect_equal(nrow(result), 2)
  expect_equal(result$channel, c(1L, 2L))
})

test_that("ecgDFA returns finite alpha1 AND alpha2 on a long RR series (log-spaced scales)", {
  # Regression: with LINEARLY-spaced scales, a long series (n/4 >> 64) left the
  # short range [4,16] with a single scale and the long range [16,64] empty, so
  # alpha1/alpha2 came back NA. Log-spaced scales populate both ranges.
  set.seed(2024)
  n <- 6000                                   # a realistic overnight RR length
  rr_ms <- 850 + cumsum(rnorm(n, sd = 2))
  rr <- data.frame(channel = 1L, rr_ms = rr_ms,
                   time_sec = cumsum(rr_ms / 1000))
  result <- ecgDFA(rr)
  expect_true(is.finite(result$alpha1))
  expect_true(is.finite(result$alpha2))
})

test_that("ecgDFA recovers known scaling exponents (white noise ~0.5, random walk ~1.5)", {
  set.seed(42)
  n <- 4000
  wn <- data.frame(channel = 1L, rr_ms = 850 + rnorm(n, sd = 30),
                   time_sec = cumsum(rep(0.85, n)))
  rw <- data.frame(channel = 1L, rr_ms = 850 + cumsum(rnorm(n, sd = 3)),
                   time_sec = cumsum(rep(0.85, n)))
  d_wn <- ecgDFA(wn); d_rw <- ecgDFA(rw)
  # white noise: the long-range exponent is the clean estimate (~0.5); the
  # short-range alpha1 carries the well-known small-scale positive bias.
  expect_true(d_wn$alpha2 > 0.4 && d_wn$alpha2 < 0.65)
  # random walk (integrated white noise): alpha ~ 1.5 at both ranges
  expect_true(d_rw$alpha1 > 1.3 && d_rw$alpha1 < 1.7)
  expect_true(d_rw$alpha2 > 1.3 && d_rw$alpha2 < 1.7)
})

# --- Wrapper test ---

test_that("ecgHRVnonlinear merges all metrics", {
  set.seed(10)
  n <- 200
  rr_ms <- rnorm(n, mean = 850, sd = 30)
  rr <- data.frame(
    channel = rep(1L, n),
    rr_ms = rr_ms,
    time_sec = cumsum(rr_ms / 1000)
  )

  result <- ecgHRVnonlinear(rr)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expected_cols <- c("channel", "sd1", "sd2", "sd1_sd2_ratio",
                     "sample_entropy", "m", "r", "alpha1", "alpha2")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("ecgHRVnonlinear handles multi-channel", {
  set.seed(20)
  n <- 150
  rr <- data.frame(
    channel = c(rep(1L, n), rep(2L, n)),
    rr_ms = c(rnorm(n, 850, 30), rnorm(n, 900, 40)),
    time_sec = c(cumsum(rep(0.85, n)), cumsum(rep(0.9, n)))
  )

  result <- ecgHRVnonlinear(rr)
  expect_equal(nrow(result), 2)
  expect_true(all(c(1L, 2L) %in% result$channel))
})

# --- Input validation tests ---

test_that("all nonlinear functions validate input", {
  expect_error(ecgHRVpoincare("not_df"))
  expect_error(ecgHRVpoincare(data.frame(x = 1)))
  expect_error(ecgSampleEntropy("not_df"))
  expect_error(ecgSampleEntropy(data.frame(x = 1)))
  expect_error(ecgDFA("not_df"))
  expect_error(ecgDFA(data.frame(x = 1)))
})

test_that("ecgSampleEntropy validates parameters", {
  rr <- data.frame(
    channel = rep(1L, 50),
    rr_ms = rep(800, 50),
    time_sec = cumsum(rep(0.8, 50))
  )

  expect_error(ecgSampleEntropy(rr, m = -1))
  expect_error(ecgSampleEntropy(rr, r_factor = 0))
  expect_error(ecgSampleEntropy(rr, r_factor = -0.1))
})

# --- WS5-17: rolling DFA-a1 + threshold crossings ---

test_that("rolling DFA-a1 approximates whole-series alpha1 and counts windows", {
  set.seed(21)
  n <- 1200
  x <- numeric(n); for (i in 2:n) x[i] <- 0.5 * x[i - 1] + rnorm(1)
  rr <- data.frame(channel = 1L, rr_ms = 800 + 10 * x)

  whole <- PhysioECG:::.dfa_alpha1(rr$rr_ms)
  traj <- ecgDFArolling(rr, window_beats = 300, step_beats = 60)

  expect_equal(nrow(traj), length(seq(1L, n - 300 + 1L, 60L)))   # window count
  expect_true(all(is.finite(traj$alpha1)))
  expect_lt(abs(mean(traj$alpha1) - whole) / whole, 0.10)        # rolling ~ whole
})

test_that("dfaThresholdCrossings detects AT1/AT2 on a descending ramp", {
  a1 <- seq(1.2, 0.3, length.out = 20)
  tc <- dfaThresholdCrossings(a1)
  expect_length(tc$at1_crossings, 1L)   # crosses 0.75 once
  expect_length(tc$at2_crossings, 1L)   # crosses 0.50 once
  i <- tc$at1_crossings
  expect_true(a1[i] >= 0.75 && a1[i + 1] < 0.75)
})

test_that("dfaThresholdCrossings indices map to the original (NA-containing) vector", {
  a1 <- c(NA, NA, 1.2, 1.0, 0.8, 0.6, 0.4)
  tc <- dfaThresholdCrossings(a1)
  expect_equal(tc$at1_crossings, 5L)   # 0.8 (idx5) -> 0.6 (idx6) crosses 0.75
  i <- tc$at1_crossings
  expect_true(a1[i] >= 0.75 && a1[i + 1] < 0.75)
})
