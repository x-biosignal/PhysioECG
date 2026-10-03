library(testthat)
library(PhysioECG)

make_rr <- function(n = 500, sd = 12, seed = 1) {
  set.seed(seed)
  t <- cumsum(rep(0.8, n))
  rr <- data.frame(channel = 1L,
                   rr_ms = 800 + 25 * sin(2 * pi * 0.1 * t) + rnorm(n, sd = sd))
  rr$time_sec <- cumsum(rr$rr_ms) / 1000
  rr
}

test_that("ecgHRVtimevarying returns the expected trajectory columns", {
  rr <- make_rr()
  tv <- ecgHRVtimevarying(rr, window_sec = 120, step_sec = 30)
  expect_s3_class(tv, "data.frame")
  expect_true(all(c("channel", "window", "time_start", "time_center",
                    "n_beats", "sdnn", "rmssd", "lf", "hf", "lf_hf") %in%
                    names(tv)))
  expect_true(all(tv$n_beats >= 20))
})

test_that("a single full-length window equals the static HRV functions", {
  rr <- make_rr()
  dur <- max(rr$time_sec) - min(rr$time_sec)
  tv <- ecgHRVtimevarying(rr, window_sec = dur + 10, step_sec = 30)
  expect_equal(nrow(tv), 1L)

  st_t <- ecgHRVtime(rr, rhythm_check = FALSE)
  st_f <- ecgHRVfreq(rr, method = "ar", rhythm_check = FALSE)
  expect_equal(tv$sdnn, st_t$sdnn, tolerance = 1e-8)
  expect_equal(tv$rmssd, st_t$rmssd, tolerance = 1e-8)
  expect_equal(tv$lf, st_f$lf, tolerance = 1e-8)
  expect_equal(tv$hf, st_f$hf, tolerance = 1e-8)
  expect_equal(tv$lf_hf, st_f$lf_hf_ratio, tolerance = 1e-8)
})

test_that("the number of windows matches (duration - window) / step", {
  rr <- make_rr()
  dur <- max(rr$time_sec) - min(rr$time_sec)
  for (W in c(90, 120, 180)) {
    for (S in c(15, 30, 45)) {
      tv <- ecgHRVtimevarying(rr, window_sec = W, step_sec = S)
      expect_equal(nrow(tv), floor((dur - W) / S) + 1)
    }
  }
})

test_that("the HF trajectory rises at a mid-recording HF increase", {
  set.seed(2)
  n <- 700
  t <- cumsum(rep(0.8, n))
  mid <- t[n] / 2
  hf_amp <- ifelse(t > mid, 35, 2)             # respiratory HF turns on midway
  rr <- data.frame(channel = 1L,
    rr_ms = 800 + 20 * sin(2 * pi * 0.1 * t) +
      hf_amp * sin(2 * pi * 0.25 * t) + rnorm(n, sd = 6))
  rr$time_sec <- cumsum(rr$rr_ms) / 1000

  tv <- ecgHRVtimevarying(rr, window_sec = 100, step_sec = 20)
  hf_before <- mean(tv$hf[tv$time_center < mid * 0.8], na.rm = TRUE)
  hf_after <- mean(tv$hf[tv$time_center > mid * 1.2], na.rm = TRUE)
  expect_gt(hf_after, hf_before * 2)           # HF clearly rises

  # the rise happens near the transition, not at the start
  rise_window <- tv$time_center[which(tv$hf > hf_before * 3)[1]]
  expect_gt(rise_window, mid * 0.5)
})

test_that("windows with too few beats yield NA and multi-channel works", {
  rr <- make_rr(n = 300)
  # tiny window -> some windows have < min_beats
  tv <- ecgHRVtimevarying(rr, window_sec = 12, step_sec = 12, min_beats = 20)
  expect_true(any(is.na(tv$sdnn)))
  expect_true(all(tv$n_beats[is.na(tv$sdnn)] < 20))

  rr2 <- rbind(make_rr(n = 300, seed = 1),
               transform(make_rr(n = 300, seed = 2), channel = 2L))
  tv2 <- ecgHRVtimevarying(rr2, window_sec = 100, step_sec = 40)
  expect_setequal(unique(tv2$channel), c(1L, 2L))
})
