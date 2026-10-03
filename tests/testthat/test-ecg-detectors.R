add_noise <- function(pe, sd, seed) {
  set.seed(seed)
  a <- SummarizedExperiment::assay(pe, "raw")
  SummarizedExperiment::assay(pe, "raw") <- a + matrix(rnorm(length(a), sd = sd),
                                                       nrow = nrow(a))
  pe
}

test_that("all detectors agree on beat count and handle inverted leads", {
  pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)  # ~11-12 beats
  ref <- nrow(ecgDetectRpeaks(pe, method = "pan_tompkins"))
  for (m in c("hamilton", "elgendi", "christov")) {
    n <- nrow(ecgDetectRpeaks(pe, method = m))
    expect_lte(abs(n - ref), 1L)          # within +/-1 beat of Pan-Tompkins
  }

  peA <- pe
  SummarizedExperiment::assay(peA, "raw") <- -SummarizedExperiment::assay(pe, "raw")
  for (m in c("hamilton", "elgendi", "christov")) {
    p <- ecgDetectRpeaks(peA, method = m)
    expect_true(nrow(p) >= 8 && nrow(p) <= 15)
    expect_true(all(p$amplitude < 0))     # amplitude reported from inverted lead
  }
})

test_that("bSQI is ~1 on clean ECG and drops on noisy ECG", {
  set.seed(1)
  pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
  expect_gt(ecgBeatSQI(pe)$bSQI, 0.95)                       # clean -> agreement
  expect_lt(ecgBeatSQI(add_noise(pe, 1.0, 11))$bSQI, 0.7)    # noisy -> disagreement
})

test_that("ecgBeatSQIgate accepts clean and rejects noisy channels", {
  pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
  g <- ecgBeatSQIgate(pe, threshold = 0.8)
  expect_true("accept" %in% names(g))
  expect_true(g$accept[1])
  expect_false(ecgBeatSQIgate(add_noise(pe, 1.0, 2), threshold = 0.8)$accept[1])
})

test_that("bandpass keeps its passband: gain peaks at the band centre, not 2x it", {
  fb <- PhysioECG:::.fir_bandpass
  sr <- 500; n <- 5000; t <- (0:(n - 1)) / sr
  gain <- function(f) {
    y <- fb(sin(2 * pi * f * t), sr, 5, 15, order = 65)
    sqrt(2 * mean(y[500:(n - 500)]^2))            # ~peak amplitude of a sinusoid
  }
  expect_gt(gain(10), 0.9)                          # centre of 5-15 Hz -> unit gain
  expect_gt(gain(10), gain(20))                     # not doubled to a 10-30 band
  expect_lt(gain(30), 0.2)                          # well outside the true band
})

test_that("R-peak fiducial stays on R for a deep-S (biphasic) lead", {
  sr <- 500; n <- 5000
  sig <- numeric(n)
  locs <- seq(400, n - 200, by = 400)
  for (L in locs) { sig[L] <- 1.0; sig[L + 30] <- -1.6 }   # positive R, deeper S
  sig <- as.numeric(stats::filter(sig, rep(1 / 5, 5), sides = 2))
  sig[is.na(sig)] <- 0
  pe <- PhysioCore::PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                                     samplingRate = sr)
  for (m in c("pan_tompkins", "hamilton", "elgendi", "christov")) {
    pk <- ecgDetectRpeaks(pe, method = m)
    expect_true(all(pk$amplitude > 0),
                info = paste(m, "should report the positive R amplitude"))
    nearest <- vapply(pk$sample, function(s) min(abs(s - locs)), numeric(1))
    expect_true(all(nearest <= 5),
                info = paste(m, "fiducial should land within 5 samples of R"))
  }
})

test_that(".find_local_maxima handles a signal with a single local maximum", {
  flm <- PhysioECG:::.find_local_maxima
  expect_equal(flm(c(0, 0, 1, 0, 0), min_distance = 2L), 3L)   # was a crash
  expect_equal(flm(c(0, 0, 0, 0), min_distance = 2L), integer(0))
})

test_that(".ecg_match_beats matches one-to-one within the window", {
  match <- PhysioECG:::.ecg_match_beats
  expect_equal(match(c(100, 200, 300), c(105, 205, 295), 10), 3L)
  expect_equal(match(c(100, 200), 150, 10), 0L)              # both too far
  expect_equal(match(integer(0), c(1, 2), 10), 0L)           # empty a
})
