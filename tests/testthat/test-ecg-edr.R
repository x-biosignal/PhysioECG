library(testthat)
library(PhysioECG)

# Synthetic ECG whose R-wave amplitude is modulated by respiration at f_resp.
make_resp_ecg <- function(sr = 250, dur = 140, hr = 75, f_resp = 0.25,
                          depth = 0.35, seed = 1) {
  set.seed(seed)
  n <- sr * dur
  RR <- as.integer(sr * 60 / hr)
  sig <- rnorm(n, sd = 0.02)
  qw <- as.integer(round(0.04 * sr))
  for (p in seq(RR, n - qw, by = RR)) {
    amp <- 1.5 * (1 + depth * sin(2 * pi * f_resp * p / sr))
    idx <- seq(max(1L, p - qw), min(n, p + qw))
    sig[idx] <- sig[idx] + amp * exp(-((idx - p)^2) / (2 * (qw / 3)^2))
  }
  PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                   colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
                   samplingRate = sr)
}

test_that("EDR recovers the imposed respiratory frequency within 0.02 Hz", {
  for (fr in c(0.20, 0.25, 0.30)) {
    pe <- make_resp_ecg(f_resp = fr)
    pk <- ecgDetectRpeaks(pe)
    for (m in c("amplitude", "area", "slope")) {
      edr <- ecgEDR(pe, pk, method = m)
      expect_true(all(c("method", "beats", "edr", "resp_rate") %in%
                        names(edr)))
      expect_equal(edr$resp_rate$resp_rate_hz, fr, tolerance = 0.02)
      # respiratory rate within 1 breath/min
      expect_equal(edr$resp_rate$resp_rate_bpm, fr * 60, tolerance = 1)
    }
  }
})

test_that("EDR returns per-beat and resampled surrogate series", {
  pe <- make_resp_ecg(f_resp = 0.25)
  pk <- ecgDetectRpeaks(pe)
  edr <- ecgEDR(pe, pk, method = "area")
  expect_s3_class(edr$beats, "data.frame")
  expect_true(all(c("channel", "time_sec", "feature") %in% names(edr$beats)))
  expect_true(all(c("channel", "time_sec", "edr") %in% names(edr$edr)))
  expect_equal(nrow(edr$beats), nrow(pk))
  expect_gt(nrow(edr$edr), nrow(pk))          # upsampled onto a uniform grid
})

test_that("RSA peak-valley amplitude increases with HF modulation depth", {
  depths <- c(5, 15, 30, 50)
  rsa <- vapply(depths, function(d) {
    set.seed(9)
    t <- cumsum(rep(0.8, 400))
    rr <- data.frame(rr_ms = 800 + d * sin(2 * pi * 0.25 * t) +
                       rnorm(400, sd = 4), time_sec = t)
    ecgRSA(rr, method = "peak_valley")$rsa
  }, numeric(1))
  expect_true(all(diff(rsa) > 0))             # monotone increasing
})

test_that("RSA Porges-Bohrer variance also increases monotonically", {
  depths <- c(5, 15, 30, 50)
  rsa <- vapply(depths, function(d) {
    set.seed(9)
    t <- cumsum(rep(0.8, 400))
    rr <- data.frame(rr_ms = 800 + d * sin(2 * pi * 0.25 * t) +
                       rnorm(400, sd = 4), time_sec = t)
    ecgRSA(rr, method = "porges_bohrer")$rsa
  }, numeric(1))
  expect_true(all(diff(rsa) > 0))
})

test_that("ecgRSA accepts a numeric RR vector and reports cycle count", {
  set.seed(1)
  v <- 800 + 40 * sin(2 * pi * 0.25 * cumsum(rep(0.8, 300))) + rnorm(300, 5)
  # a numeric vector derives beat times as cumsum(rr)/1000; supplying the same
  # time_sec to the data.frame form must give the identical RSA
  res_df <- ecgRSA(data.frame(rr_ms = v, time_sec = cumsum(v) / 1000))
  res_vec <- ecgRSA(v)
  expect_equal(res_df$rsa, res_vec$rsa, tolerance = 1e-6)
  expect_gt(res_df$n_cycles, 0)
})
