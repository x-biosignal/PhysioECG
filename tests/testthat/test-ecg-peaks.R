library(testthat)
library(PhysioECG)

test_that("ecgDetectRpeaks detects peaks in synthetic ECG", {
  pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
  peaks <- ecgDetectRpeaks(pe)

  expect_s3_class(peaks, "data.frame")
  expect_true(all(c("channel", "sample", "time_sec", "amplitude") %in% names(peaks)))
  # At 72 bpm over 10 sec, expect ~11-12 peaks

  expect_true(nrow(peaks) >= 8 && nrow(peaks) <= 15)
  # Peaks should be positive
  expect_true(all(peaks$amplitude > 0))
})

test_that("ecgDetectRpeaks respects refractory period", {
  pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
  peaks <- ecgDetectRpeaks(pe, refractory_ms = 200)

  if (nrow(peaks) >= 2) {
    diffs <- diff(peaks$sample)
    min_refractory <- as.integer(round(200 / 1000 * 500))
    expect_true(all(diffs >= min_refractory))
  }
})

test_that("ecgDetectRpeaks works with multiple channels", {
  pe <- make_ecg(n_time = 5000, n_channels = 3, sr = 500, heart_rate = 72)
  peaks <- ecgDetectRpeaks(pe)

  expect_s3_class(peaks, "data.frame")
  expect_true(all(unique(peaks$channel) %in% 1:3))
  # Each channel should have detected peaks
  for (ch in 1:3) {
    ch_peaks <- peaks[peaks$channel == ch, ]
    expect_true(nrow(ch_peaks) >= 5)
  }
})

test_that("ecgRRintervals computes intervals correctly", {
  pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 60)
  peaks <- ecgDetectRpeaks(pe)
  rr <- ecgRRintervals(pe, peaks)

  expect_s3_class(rr, "data.frame")
  expect_true(all(c("channel", "rr_ms", "time_sec") %in% names(rr)))
  # At 60 bpm, RR should be ~1000ms
  expect_true(all(rr$rr_ms > 800 & rr$rr_ms < 1200))
})

test_that("ecgRRintervals returns empty data.frame for single peak", {
  pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 60)
  single_peak <- data.frame(channel = 1L, sample = 500L, time_sec = 1.0,
                            amplitude = 1.5)
  rr <- ecgRRintervals(pe, single_peak)

  expect_s3_class(rr, "data.frame")
  expect_equal(nrow(rr), 0)
})

test_that("ecgDetectRpeaks validates input", {
  expect_error(ecgDetectRpeaks("not_pe"))
})

test_that("ecgRRintervals validates input", {
  expect_error(ecgRRintervals("not_pe", data.frame()))
  pe <- make_ecg(n_time = 1000, sr = 500)
  expect_error(ecgRRintervals(pe, data.frame(a = 1)))
})

test_that("ecgDetectRpeaks handles amplitude-varying ECG (adaptive threshold)", {
  pe <- make_ecg_amplitude_varying(n_time = 10000, sr = 500, heart_rate = 72)
  peaks <- ecgDetectRpeaks(pe)

  expect_s3_class(peaks, "data.frame")
  # At 72 bpm over 20 sec, expect ~23-24 beats; adaptive threshold should

  # catch most of them including the low-amplitude middle section.
  # With 10000 samples at 500 Hz = 20 sec, ~24 beats expected.
  expect_true(nrow(peaks) >= 18,
              label = sprintf("Expected >=18 peaks for amplitude-varying ECG, got %d", nrow(peaks)))
  # Peaks should span the full signal, not just the high-amplitude ends
  signal_duration <- 10000 / 500
  expect_true(max(peaks$time_sec) > signal_duration * 0.7,
              label = "Peaks should be detected in the latter portion of the signal")
  expect_true(min(peaks$time_sec) < signal_duration * 0.3,
              label = "Peaks should be detected in the early portion of the signal")
})

test_that("ecgDetectRpeaks handles inverted ECG signal", {
  pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
  # Get reference peaks from normal signal
  peaks_normal <- ecgDetectRpeaks(pe)

  # Invert the signal
  inv_data <- -SummarizedExperiment::assay(pe, "raw")
  SummarizedExperiment::assay(pe, "raw") <- inv_data
  peaks_inv <- ecgDetectRpeaks(pe)

  expect_s3_class(peaks_inv, "data.frame")
  # Should detect a similar number of peaks as the normal signal
  expect_true(nrow(peaks_inv) >= 8 && nrow(peaks_inv) <= 15,
              label = sprintf("Inverted signal: expected 8-15 peaks, got %d", nrow(peaks_inv)))
  # Amplitudes from the inverted signal should be negative (original data is negative)
  expect_true(all(peaks_inv$amplitude < 0),
              label = "Amplitudes should be negative for inverted signal")
  # Peak locations should be close to the normal signal's peak locations
  expect_equal(nrow(peaks_inv), nrow(peaks_normal), tolerance = 2)
})
