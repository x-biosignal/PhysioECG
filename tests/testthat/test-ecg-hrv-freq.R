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
