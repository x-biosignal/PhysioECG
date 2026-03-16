library(testthat)
library(PhysioECG)

# --- ecgDelineate tests ---

test_that("ecgDelineate returns correct structure", {
  result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72)
  pe <- result$pe
  fiducials <- result$fiducials

  peaks <- data.frame(
    channel = rep(1L, nrow(fiducials)),
    sample = fiducials$r_peak,
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)

  expect_s3_class(delin, "data.frame")
  expected_cols <- c("channel", "beat", "r_peak", "qrs_onset", "qrs_offset",
                     "qrs_duration_ms", "p_peak", "t_peak")
  expect_true(all(expected_cols %in% names(delin)))
  expect_equal(nrow(delin), nrow(fiducials))
})

test_that("ecgDelineate detects QRS boundaries within tolerance", {
  set.seed(123)
  result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72,
                            noise_sd = 0.01)
  pe <- result$pe
  fiducials <- result$fiducials

  peaks <- data.frame(
    channel = rep(1L, nrow(fiducials)),
    sample = fiducials$r_peak,
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)

  # QRS onset should be within +/-10 samples of known values
  for (i in seq_len(nrow(fiducials))) {
    expect_true(
      abs(delin$qrs_onset[i] - fiducials$qrs_onset[i]) <= 10,
      info = sprintf("Beat %d: QRS onset %d vs expected %d (diff=%d)",
                     i, delin$qrs_onset[i], fiducials$qrs_onset[i],
                     delin$qrs_onset[i] - fiducials$qrs_onset[i])
    )
  }

  # QRS offset should be within +/-10 samples of known values
  for (i in seq_len(nrow(fiducials))) {
    expect_true(
      abs(delin$qrs_offset[i] - fiducials$qrs_offset[i]) <= 10,
      info = sprintf("Beat %d: QRS offset %d vs expected %d (diff=%d)",
                     i, delin$qrs_offset[i], fiducials$qrs_offset[i],
                     delin$qrs_offset[i] - fiducials$qrs_offset[i])
    )
  }
})

test_that("ecgDelineate detects P and T wave peaks", {
  set.seed(123)
  result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72,
                            noise_sd = 0.01)
  pe <- result$pe
  fiducials <- result$fiducials

  peaks <- data.frame(
    channel = rep(1L, nrow(fiducials)),
    sample = fiducials$r_peak,
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)

  # P and T peaks should not be NA for the PQRST data
  expect_true(all(!is.na(delin$p_peak)))
  expect_true(all(!is.na(delin$t_peak)))

  # P peak should be within +/-10 samples of known location
  for (i in seq_len(nrow(fiducials))) {
    expect_true(
      abs(delin$p_peak[i] - fiducials$p_peak[i]) <= 10,
      info = sprintf("Beat %d: P peak %d vs expected %d",
                     i, delin$p_peak[i], fiducials$p_peak[i])
    )
  }

  # T peak should be within +/-10 samples of known location
  for (i in seq_len(nrow(fiducials))) {
    expect_true(
      abs(delin$t_peak[i] - fiducials$t_peak[i]) <= 10,
      info = sprintf("Beat %d: T peak %d vs expected %d",
                     i, delin$t_peak[i], fiducials$t_peak[i])
    )
  }
})

test_that("ecgDelineate QRS duration is physiologically reasonable", {
  result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72)
  pe <- result$pe
  fiducials <- result$fiducials

  peaks <- data.frame(
    channel = rep(1L, nrow(fiducials)),
    sample = fiducials$r_peak,
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)

  # QRS duration should be between 20ms and 200ms (physiological range)
  expect_true(all(delin$qrs_duration_ms > 20))
  expect_true(all(delin$qrs_duration_ms < 200))
})

test_that("ecgDelineate works with multiple channels", {
  result <- make_ecg_pqrst(n_time = 10000, n_channels = 2, sr = 500,
                            heart_rate = 72)
  pe <- result$pe
  fiducials <- result$fiducials

  peaks <- data.frame(
    channel = c(rep(1L, nrow(fiducials)), rep(2L, nrow(fiducials))),
    sample = rep(fiducials$r_peak, 2),
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)

  expect_true(all(c(1L, 2L) %in% delin$channel))
  ch1 <- delin[delin$channel == 1, ]
  ch2 <- delin[delin$channel == 2, ]
  expect_equal(nrow(ch1), nrow(fiducials))
  expect_equal(nrow(ch2), nrow(fiducials))
})

test_that("ecgDelineate returns empty data.frame for empty peaks", {
  result <- make_ecg_pqrst(n_time = 10000, sr = 500)
  pe <- result$pe

  empty_peaks <- data.frame(channel = integer(0), sample = integer(0))
  delin <- ecgDelineate(pe, empty_peaks)

  expect_s3_class(delin, "data.frame")
  expect_equal(nrow(delin), 0)
})

test_that("ecgDelineate validates input", {
  expect_error(ecgDelineate("not_pe", data.frame(channel = 1, sample = 100)))

  result <- make_ecg_pqrst(n_time = 10000, sr = 500)
  pe <- result$pe
  expect_error(ecgDelineate(pe, data.frame(a = 1, b = 2)),
               "Missing required columns")
})

# --- ecgIntervals tests ---

test_that("ecgIntervals returns correct structure", {
  result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72)
  pe <- result$pe
  fiducials <- result$fiducials

  peaks <- data.frame(
    channel = rep(1L, nrow(fiducials)),
    sample = fiducials$r_peak,
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)
  intervals <- ecgIntervals(delin, samplingRate(pe))

  expect_s3_class(intervals, "data.frame")
  expected_cols <- c("channel", "beat", "pr_ms", "qt_ms", "qtc_ms",
                     "qrs_ms", "rr_ms")
  expect_true(all(expected_cols %in% names(intervals)))
  expect_equal(nrow(intervals), nrow(delin))
})

test_that("ecgIntervals computes QRS duration correctly", {
  result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72)
  pe <- result$pe
  fiducials <- result$fiducials
  sr <- samplingRate(pe)

  peaks <- data.frame(
    channel = rep(1L, nrow(fiducials)),
    sample = fiducials$r_peak,
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)
  intervals <- ecgIntervals(delin, sr)

  # QRS duration from intervals should match delineation
  expect_equal(intervals$qrs_ms, delin$qrs_duration_ms)
})

test_that("ecgIntervals computes RR intervals for 72 bpm", {
  result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72)
  pe <- result$pe
  fiducials <- result$fiducials
  sr <- samplingRate(pe)

  peaks <- data.frame(
    channel = rep(1L, nrow(fiducials)),
    sample = fiducials$r_peak,
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)
  intervals <- ecgIntervals(delin, sr)

  # At 72 bpm, RR should be ~833ms
  rr_valid <- intervals$rr_ms[!is.na(intervals$rr_ms)]
  expect_true(length(rr_valid) > 0)
  expect_true(all(abs(rr_valid - 833.3) < 50))

  # Last beat should have NA RR
  expect_true(is.na(intervals$rr_ms[nrow(intervals)]))
})

test_that("ecgIntervals computes QTc with Bazett formula", {
  result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72)
  pe <- result$pe
  fiducials <- result$fiducials
  sr <- samplingRate(pe)

  peaks <- data.frame(
    channel = rep(1L, nrow(fiducials)),
    sample = fiducials$r_peak,
    stringsAsFactors = FALSE
  )

  delin <- ecgDelineate(pe, peaks)
  intervals <- ecgIntervals(delin, sr)

  # Verify Bazett formula: QTc = QT / sqrt(RR_sec)
  for (i in seq_len(nrow(intervals))) {
    if (!is.na(intervals$qtc_ms[i]) && !is.na(intervals$qt_ms[i]) &&
        !is.na(intervals$rr_ms[i])) {
      rr_sec <- intervals$rr_ms[i] / 1000
      expected_qtc <- intervals$qt_ms[i] / sqrt(rr_sec)
      expect_equal(intervals$qtc_ms[i], expected_qtc, tolerance = 0.01)
    }
  }
})

test_that("ecgIntervals handles missing P and T waves", {
  # Create delineation with NA p_peak, t_peak, and t_end
  delin <- data.frame(
    channel = c(1L, 1L),
    beat = c(1L, 2L),
    r_peak = c(500L, 1000L),
    qrs_onset = c(480L, 980L),
    qrs_offset = c(520L, 1020L),
    qrs_duration_ms = c(80, 80),
    p_peak = c(NA_integer_, NA_integer_),
    t_peak = c(NA_integer_, NA_integer_),
    t_end = c(NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )

  intervals <- ecgIntervals(delin, sr = 500)

  expect_true(all(is.na(intervals$pr_ms)))
  expect_true(all(is.na(intervals$qt_ms)))
  # QTc should also be NA when QT is NA
  expect_true(all(is.na(intervals$qtc_ms)))
  # QRS and RR should still be valid
  expect_equal(intervals$qrs_ms, c(80, 80))
  expect_equal(intervals$rr_ms[1], 1000)
  expect_true(is.na(intervals$rr_ms[2]))
})

test_that("ecgIntervals returns empty data.frame for empty delineation", {
  empty_delin <- data.frame(
    channel = integer(0), beat = integer(0), r_peak = integer(0),
    qrs_onset = integer(0), qrs_offset = integer(0),
    qrs_duration_ms = numeric(0), p_peak = integer(0),
    t_peak = integer(0), t_end = integer(0)
  )

  intervals <- ecgIntervals(empty_delin, sr = 500)

  expect_s3_class(intervals, "data.frame")
  expect_equal(nrow(intervals), 0)
})

test_that("ecgIntervals validates input", {
  expect_error(ecgIntervals(data.frame(a = 1), sr = 500),
               "Missing required columns")
  expect_error(ecgIntervals(data.frame(channel = 1, beat = 1, r_peak = 100,
                                        qrs_onset = 90, qrs_offset = 110,
                                        p_peak = 70, t_peak = 200,
                                        t_end = 250),
                            sr = -1))
})
