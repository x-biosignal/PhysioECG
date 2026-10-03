bp <- function(sig, lo, hi, fs) {
  n <- length(sig)
  p <- Mod(stats::fft(sig - mean(sig)))^2 / n
  f <- (seq_len(n) - 1) * fs / n
  half <- seq_len(n %/% 2)
  sum(p[half][f[half] >= lo & f[half] < hi])
}

test_that("smoothnessPriorsDetrend removes slow trend, keeps HF, near-zero mean", {
  fs <- 4
  t <- seq(0, 300, by = 1 / fs)                 # 5 min at 4 Hz
  hf <- sin(2 * pi * 0.25 * t)                  # 0.25 Hz HF oscillation
  trend <- 0.02 * t + 3 * sin(2 * pi * 0.005 * t)  # linear + deep-VLF
  x <- trend + hf

  xd <- smoothnessPriorsDetrend(x, lambda = 500, fs = fs)
  expect_lt(abs(mean(xd)), 1e-6)

  vlf_ratio <- bp(xd, 0.003, 0.04, fs) / bp(x, 0.003, 0.04, fs)
  hf_ratio  <- bp(xd, 0.15,  0.40, fs) / bp(x, 0.15,  0.40, fs)
  expect_lt(vlf_ratio, 0.10)   # VLF power reduced by >90%
  expect_gt(hf_ratio, 0.90)    # HF power retained >90%
})

test_that("higher lambda lowers the cutoff frequency (monotone)", {
  cut <- function(l) attr(smoothnessPriorsDetrend(rnorm(200), lambda = l, fs = 4),
                          "cutoff_hz")
  expect_gt(cut(100), cut(500))
  expect_gt(cut(500), cut(2000))
})

test_that("ecgHRVfreq detrend=TRUE reduces VLF relative to detrend=FALSE", {
  set.seed(5)
  fs_rr <- 1
  n <- 600
  time_sec <- cumsum(rep(1, n))
  rr_ms <- 800 + 60 * sin(2 * pi * 0.004 * time_sec) +   # strong VLF drift
           20 * sin(2 * pi * 0.25 * time_sec) + rnorm(n, sd = 3)
  rr <- data.frame(channel = 1L, rr_ms = rr_ms, time_sec = time_sec)

  plain <- ecgHRVfreq(rr, method = "welch", detrend = FALSE)
  det   <- ecgHRVfreq(rr, method = "welch", detrend = TRUE)
  expect_lt(det$vlf, plain$vlf)
})

test_that("cutoff_hz is the -3 dB (half-power) point of the high-pass filter", {
  lambda <- 500; fs <- 4
  cut <- attr(smoothnessPriorsDetrend(rnorm(50), lambda = lambda, fs = fs), "cutoff_hz")
  w <- 2 * pi * cut / fs
  s <- 2 - 2 * cos(w)
  L <- 1 / (1 + lambda^2 * s^2)   # low-pass (trend) amplitude gain
  expect_equal((1 - L)^2, 0.5, tolerance = 1e-6)   # high-pass power = 1/2
})

test_that("smoothnessPriorsDetrend attaches attributes even for short input", {
  x <- smoothnessPriorsDetrend(c(1, 2), lambda = 500, fs = 4)
  expect_equal(attr(x, "lambda"), 500)
  expect_true("cutoff_hz" %in% names(attributes(x)))
})
