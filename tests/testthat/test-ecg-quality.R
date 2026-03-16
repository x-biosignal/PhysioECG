library(testthat)
library(PhysioECG)

test_that("ecgQualityCheck detects ectopic beats", {
  rr <- data.frame(
    channel = rep(1L, 10),
    rr_ms = c(800, 850, 820, 500, 1200, 830, 840, 810, 850, 820),
    time_sec = cumsum(c(0, c(800, 850, 820, 500, 1200, 830, 840, 810, 850) / 1000))
  )

  result <- ecgQualityCheck(rr, threshold_ms = 200)

  expect_true("is_ectopic" %in% names(result))
  # The 500ms and 1200ms intervals should be marked as ectopic
  expect_true(result$is_ectopic[4])
  expect_true(result$is_ectopic[5])
  # Normal intervals should not be ectopic
  expect_false(result$is_ectopic[1])
  expect_false(result$is_ectopic[6])
})

test_that("ecgRRcorrect interpolate replaces ectopic values", {
  rr <- data.frame(
    channel = rep(1L, 8),
    rr_ms = c(800, 850, 500, 1200, 830, 840, 810, 850),
    time_sec = cumsum(c(0, c(800, 850, 500, 1200, 830, 840, 810) / 1000)),
    is_ectopic = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  )

  corrected <- ecgRRcorrect(rr, method = "interpolate")

  expect_false("is_ectopic" %in% names(corrected))
  expect_equal(nrow(corrected), 8)
  # Corrected values should be closer to normal range
  expect_true(all(corrected$rr_ms > 700 & corrected$rr_ms < 1000))
})

test_that("ecgRRcorrect remove drops ectopic rows", {
  rr <- data.frame(
    channel = rep(1L, 6),
    rr_ms = c(800, 500, 850, 1200, 830, 840),
    time_sec = cumsum(c(0, c(800, 500, 850, 1200, 830) / 1000)),
    is_ectopic = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE)
  )

  corrected <- ecgRRcorrect(rr, method = "remove")

  expect_equal(nrow(corrected), 4)
  expect_false("is_ectopic" %in% names(corrected))
})

test_that("ecgQualityCheck validates input", {
  expect_error(ecgQualityCheck("not_df"))
})

test_that("ecgRRcorrect validates input", {
  rr <- data.frame(channel = 1L, rr_ms = 800, time_sec = 0)
  expect_error(ecgRRcorrect(rr))  # missing is_ectopic
})


# -----------------------------------------------------------------------
# Tests for ecgSignalQuality
# -----------------------------------------------------------------------

test_that("ecgSignalQuality returns correct structure", {
  pe <- make_ecg(n_time = 5000, n_channels = 2, sr = 500)
  result <- ecgSignalQuality(pe)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(c("channel", "snr_db", "baseline_wander",
                     "saturation_ratio", "quality_score") %in% names(result)))
})

test_that("ecgSignalQuality SNR distinguishes clean vs noisy signal", {
  sr <- 500
  n <- 5000

  # Clean ECG-like signal (low noise)
  pe_clean <- make_ecg(n_time = n, n_channels = 1, sr = sr)
  q_clean <- ecgSignalQuality(pe_clean)

  # Noisy version: add large amount of noise
  data_noisy <- SummarizedExperiment::assay(pe_clean, "raw") +
    matrix(rnorm(n, sd = 2.0), ncol = 1)
  pe_noisy <- PhysioExperiment(
    assays = list(raw = data_noisy),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )
  q_noisy <- ecgSignalQuality(pe_noisy)

  # Clean signal should have higher SNR than noisy one
  expect_true(q_clean$snr_db > q_noisy$snr_db)
  # Both should be finite
  expect_true(is.finite(q_clean$snr_db))
  expect_true(is.finite(q_noisy$snr_db))
})

test_that("ecgSignalQuality SNR with peaks uses QRS-based estimation", {
  pe <- make_ecg(n_time = 5000, n_channels = 1, sr = 500)
  peaks <- ecgDetectRpeaks(pe)

  result_with <- ecgSignalQuality(pe, peaks = peaks)
  result_without <- ecgSignalQuality(pe)

  # Both should return valid SNR values but they may differ

  expect_true(is.finite(result_with$snr_db))
  expect_true(is.finite(result_without$snr_db))
})

test_that("ecgSignalQuality detects baseline wander", {
  sr <- 500
  n <- 5000
  t <- seq(0, (n - 1) / sr, length.out = n)

  # Clean signal (no wander)
  clean <- rnorm(n, sd = 0.05)
  data_clean <- matrix(clean, ncol = 1)
  pe_clean <- PhysioExperiment(
    assays = list(raw = data_clean),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )

  # Signal with 0.5 Hz baseline drift
  drifty <- clean + 2 * sin(2 * pi * 0.5 * t)
  data_drift <- matrix(drifty, ncol = 1)
  pe_drift <- PhysioExperiment(
    assays = list(raw = data_drift),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )

  q_clean <- ecgSignalQuality(pe_clean)
  q_drift <- ecgSignalQuality(pe_drift)

  # Drifty signal should have more baseline wander
  expect_true(q_drift$baseline_wander > q_clean$baseline_wander)
})

test_that("ecgSignalQuality detects saturation (clipping)", {
  sr <- 500
  n <- 5000
  sig <- rnorm(n, sd = 1)

  # Unclipped
  data_normal <- matrix(sig, ncol = 1)
  pe_normal <- PhysioExperiment(
    assays = list(raw = data_normal),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )

  # Clipped at +/- 1
  sig_clipped <- pmax(pmin(sig, 1), -1)
  data_clipped <- matrix(sig_clipped, ncol = 1)
  pe_clipped <- PhysioExperiment(
    assays = list(raw = data_clipped),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )

  q_normal <- ecgSignalQuality(pe_normal)
  q_clipped <- ecgSignalQuality(pe_clipped)

  # Clipped signal should have higher saturation ratio
  expect_true(q_clipped$saturation_ratio > q_normal$saturation_ratio)
})

test_that("ecgSignalQuality quality_score is in [0, 1]", {
  pe <- make_ecg(n_time = 5000, n_channels = 3, sr = 500)
  result <- ecgSignalQuality(pe)

  expect_true(all(result$quality_score >= 0))
  expect_true(all(result$quality_score <= 1))
})

test_that("ecgSignalQuality validates input", {
  expect_error(ecgSignalQuality("not_pe"))
  pe <- make_ecg(n_time = 1000, sr = 500)
  expect_error(ecgSignalQuality(pe, peaks = "not_df"))
})


# -----------------------------------------------------------------------
# Tests for ecgBaselineCorrect
# -----------------------------------------------------------------------

test_that("ecgBaselineCorrect highpass removes low-frequency drift", {
  sr <- 500
  n <- 10000
  t <- seq(0, (n - 1) / sr, length.out = n)

  # ECG-like signal with large slow baseline drift
  set.seed(42)
  ecg_sig <- rnorm(n, sd = 0.05)
  drift <- 3 * sin(2 * pi * 0.1 * t)  # 0.1 Hz drift (well below 0.5 Hz cutoff)
  data <- matrix(ecg_sig + drift, ncol = 1)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )

  pe_corr <- ecgBaselineCorrect(pe, method = "highpass", cutoff = 0.5)

  expect_true("baseline_corrected" %in% SummarizedExperiment::assayNames(pe_corr))

  corrected_sig <- SummarizedExperiment::assay(pe_corr, "baseline_corrected")[, 1]
  # RMS of corrected signal should be much less than RMS of drift
  rms_drift <- sqrt(mean(drift^2))
  rms_corrected <- sqrt(mean(corrected_sig^2))
  expect_true(rms_corrected < rms_drift * 0.5)
})

test_that("ecgBaselineCorrect median method works", {
  sr <- 500
  n <- 10000
  t <- seq(0, (n - 1) / sr, length.out = n)

  set.seed(42)
  ecg_sig <- rnorm(n, sd = 0.05)
  drift <- 2 * sin(2 * pi * 0.1 * t)  # 0.1 Hz drift
  data <- matrix(ecg_sig + drift, ncol = 1)
  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )

  pe_corr <- ecgBaselineCorrect(pe, method = "median", cutoff = 0.5)

  expect_true("baseline_corrected" %in% SummarizedExperiment::assayNames(pe_corr))

  corrected_sig <- SummarizedExperiment::assay(pe_corr, "baseline_corrected")[, 1]
  rms_drift <- sqrt(mean(drift^2))
  rms_corrected <- sqrt(mean(corrected_sig^2))
  expect_true(rms_corrected < rms_drift * 0.5)
})

test_that("ecgBaselineCorrect preserves original assay", {
  pe <- make_ecg(n_time = 3000, sr = 500)
  original <- SummarizedExperiment::assay(pe, "raw")

  pe_corr <- ecgBaselineCorrect(pe, method = "highpass")

  # Original raw assay should be unchanged
  expect_equal(SummarizedExperiment::assay(pe_corr, "raw"), original)
  # New assay should exist
  expect_true("baseline_corrected" %in% SummarizedExperiment::assayNames(pe_corr))
})

test_that("ecgBaselineCorrect custom output_assay name", {
  pe <- make_ecg(n_time = 3000, sr = 500)
  pe_corr <- ecgBaselineCorrect(pe, output_assay = "my_corrected")

  expect_true("my_corrected" %in% SummarizedExperiment::assayNames(pe_corr))
})

test_that("ecgBaselineCorrect validates input", {
  expect_error(ecgBaselineCorrect("not_pe"))
  pe <- make_ecg(n_time = 1000, sr = 500)
  expect_error(ecgBaselineCorrect(pe, method = "invalid"))
  expect_error(ecgBaselineCorrect(pe, cutoff = -1))
  expect_error(ecgBaselineCorrect(pe, cutoff = "abc"))
})

test_that("ecgBaselineCorrect works with multiple channels", {
  pe <- make_ecg(n_time = 3000, n_channels = 3, sr = 500)
  pe_corr <- ecgBaselineCorrect(pe, method = "highpass")

  corrected <- SummarizedExperiment::assay(pe_corr, "baseline_corrected")
  expect_equal(ncol(corrected), 3)
  expect_equal(nrow(corrected), 3000)
})
