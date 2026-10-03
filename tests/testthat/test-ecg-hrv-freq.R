library(testthat)
library(PhysioECG)

test_that("ecgHRVfreq welch returns valid structure", {
  # Create RR with ~0.1 Hz modulation (LF band) and ~0.25 Hz modulation (HF band)
  n <- 300
  time_sec <- cumsum(rep(0.85, n))  # ~70 bpm base
  # Add LF (0.1 Hz) and HF (0.25 Hz) modulation
  rr_ms <- 850 + 30 * sin(2 * pi * 0.1 * time_sec) + 15 * sin(2 * pi * 0.25 * time_sec)

  rr <- data.frame(channel = rep(1L, n), rr_ms = rr_ms, time_sec = time_sec)

  result <- ecgHRVfreq(rr, method = "welch")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(all(c("channel", "vlf", "lf", "hf", "lf_hf_ratio", "total_power") %in% names(result)))
  expect_true(result$lf > 0)
  expect_true(result$hf > 0)
  expect_true(result$lf_hf_ratio > 0)
  # LF should have more power than HF due to larger amplitude modulation
  expect_true(result$lf > result$hf)
})

test_that("ecgHRVfreq lomb returns valid results", {
  n <- 200
  time_sec <- cumsum(rep(0.85, n))
  rr_ms <- 850 + 20 * sin(2 * pi * 0.1 * time_sec)

  rr <- data.frame(channel = rep(1L, n), rr_ms = rr_ms, time_sec = time_sec)

  result <- ecgHRVfreq(rr, method = "lomb")

  expect_s3_class(result, "data.frame")
  expect_true(result$lf > 0)
  expect_true(result$total_power > 0)
})

test_that("ecgHRVfreq validates input", {
  expect_error(ecgHRVfreq("not_df"))
})

# --- WS5-13: AR/Burg PSD + normalized units + peak frequencies ---

test_that("normalized units sum to 100 and AR locates the HF peak", {
  set.seed(11)
  n <- 512
  time_sec <- cumsum(rep(0.8, n))
  # balanced LF (0.10 Hz) + HF (0.25 Hz) oscillations
  rr_ms <- 800 + 30 * sin(2 * pi * 0.10 * time_sec) +
           30 * sin(2 * pi * 0.25 * time_sec) + rnorm(n, sd = 3)
  rr <- data.frame(channel = 1L, rr_ms = rr_ms, time_sec = time_sec)

  for (m in c("welch", "ar", "lomb")) {
    res <- ecgHRVfreq(rr, method = m)
    expect_equal(res$lf_nu + res$hf_nu, 100, tolerance = 1e-6)  # nu identity
    expect_equal(res$hf_peak, 0.25, tolerance = 0.01)           # 0.25 Hz resolved
    expect_true(is.finite(res$lf) && is.finite(res$hf) && res$total_power > 0)
  }

  # the two non-parametric estimators agree on the LF/HF ratio within 30%
  # (a parametric AR spectrum legitimately re-weights the bands differently)
  w <- ecgHRVfreq(rr, method = "welch")
  l <- ecgHRVfreq(rr, method = "lomb")
  expect_lt(abs(w$lf_hf_ratio - l$lf_hf_ratio) / l$lf_hf_ratio, 0.30)
})

test_that("AR order auto-selection yields a finite, positive PSD", {
  set.seed(12)
  n <- 400
  time_sec <- cumsum(rep(0.85, n))
  rr_ms <- 850 + 25 * sin(2 * pi * 0.1 * time_sec) + rnorm(n, sd = 4)
  rr <- data.frame(channel = 1L, rr_ms = rr_ms, time_sec = time_sec)

  res <- ecgHRVfreq(rr, method = "ar", ar_order = NULL)
  expect_true(all(is.finite(c(res$vlf, res$lf, res$hf, res$total_power))))
  expect_gt(res$total_power, 0)
})
