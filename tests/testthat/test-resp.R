library(testthat)


test_that("breathingRate recovers an annotated sinusoid by both methods", {
  experiment <- make_rip(f = 0.25)
  peaks <- breathingRate(
    experiment, channel = "RC", method = "peaks"
  )
  spectral <- breathingRate(
    experiment, channel = "RC", method = "spectral"
  )

  expect_equal(peaks$rate_bpm, 15, tolerance = 1)
  expect_equal(peaks$n_breaths, 14L)
  expect_equal(peaks$sd_bpm, 0, tolerance = 1e-12)
  expect_equal(spectral$rate_bpm, 15, tolerance = 1e-6)
  expect_true(is.na(spectral$n_breaths))
  expect_true(is.na(spectral$sd_bpm))
})


test_that("respiratoryOnsets preserve channels and analytic timing", {
  experiment <- make_rip(f = 0.25, phase_ab = pi / 2)
  onsets <- respiratoryOnsets(experiment)

  expect_identical(sort(unique(onsets$channel)), c(1L, 2L))
  expect_identical(as.integer(table(onsets$channel)), c(14L, 14L))
  rc <- onsets[onsets$channel == 1L, ]
  expect_equal(diff(rc$onset_sample), rep(40, 13))
  expect_equal(rc$offset_time - rc$onset_time, rep(4, 14))
  expect_true(all(
    rc$onset_sample < rc$peak_sample &
      rc$peak_sample < rc$offset_sample
  ))
  expect_equal(rc$insp_amplitude, rep(2, 14), tolerance = 1e-10)
})


test_that("tidalMetrics timing and minute ventilation are analytic", {
  experiment <- make_rip(f = 0.25, a_rc = 1)
  metrics <- tidalMetrics(experiment, channel = "RC")

  expect_equal(mean(metrics$breaths$Ti), 2, tolerance = 0.1)
  expect_equal(mean(metrics$breaths$Te), 2, tolerance = 0.1)
  expect_equal(mean(metrics$breaths$duty_cycle), 0.5, tolerance = 0.02)
  expect_equal(
    mean(metrics$breaths$tidal_volume), 2, tolerance = 0.05
  )
  expect_equal(
    metrics$summary$minute_ventilation, 30, tolerance = 0.5
  )
  expect_equal(metrics$summary$n_breaths, 14L)

  scaled <- tidalMetrics(
    experiment, channel = "RC", calibration = 0.25
  )
  expect_equal(scaled$summary$tidal_volume_mean, 0.5, tolerance = 0.02)
  expect_equal(scaled$summary$minute_ventilation, 7.5, tolerance = 0.1)
})


test_that("ripCalibrate LSQ recovers Konno-Mead coefficients", {
  sr <- 10
  time <- (0:599) / sr
  rc <- sin(2 * pi * 0.20 * time)
  ab <- sin(2 * pi * 0.30 * time)
  reference <- 0.6 * rc + 0.4 * ab
  calibration <- ripCalibrate(
    rc, ab,
    reference = reference,
    sr = sr,
    method = "lsq"
  )

  expect_s3_class(calibration, "rip_calibration")
  expect_equal(unname(calibration$coef["rc"]), 0.6, tolerance = 1e-8)
  expect_equal(unname(calibration$coef["ab"]), 0.4, tolerance = 1e-8)
  expect_equal(calibration$k, 2 / 3, tolerance = 1e-8)
  expect_equal(calibration$r_squared, 1, tolerance = 1e-8)
  expect_equal(calibration$volume, reference, tolerance = 1e-10)
  expect_output(print(calibration), "method=lsq")
})


test_that("ripCalibrate QDC recovers SD ratio and absolute scale", {
  sr <- 10
  time <- (0:599) / sr
  rc <- 2 * sin(2 * pi * 0.25 * time)
  ab <- sin(2 * pi * 0.25 * time)
  calibration <- ripCalibrate(rc, ab, sr = sr, method = "qdc")

  expect_equal(calibration$k, 2, tolerance = 1e-8)
  expect_equal(unname(calibration$coef["rc"]), 1)
  expect_equal(unname(calibration$coef["ab"]), 2, tolerance = 1e-8)
  expect_true(is.na(calibration$r_squared))

  absolute <- ripCalibrate(
    rc, ab,
    sr = sr,
    method = "qdc",
    known_volume = 0.8
  )
  expect_equal(
    mean(respiratoryOnsets(absolute)$insp_amplitude),
    0.8,
    tolerance = 1e-8
  )
})


test_that("ripCalibrate resolves PhysioExperiment belt labels", {
  experiment <- make_rip(a_rc = 2, a_ab = 1)
  calibration <- ripCalibrate(experiment, ab = "AB", method = "qdc")

  expect_equal(calibration$sr, 10)
  expect_equal(calibration$k, 2, tolerance = 1e-8)
  expect_equal(length(calibration$volume), 600)
  expect_error(
    ripCalibrate(experiment, ab = "RC", method = "qdc"),
    "must be different"
  )
})


test_that("respiration calibration rejects degenerate inputs", {
  expect_error(
    ripCalibrate(rep(1, 100), rep(0, 100), sr = 10, method = "qdc"),
    "zero variance"
  )
  expect_error(
    ripCalibrate(rep(0, 100), sin(2 * pi * (0:99) / 50),
                 sr = 10, method = "qdc"),
    "ribcage belt has zero variance"
  )
  time <- (0:99) / 10
  belt <- sin(2 * pi * 0.2 * time)
  expect_error(
    ripCalibrate(
      belt, belt,
      reference = belt,
      sr = 10,
      method = "lsq"
    ),
    "collinear"
  )
  expect_error(
    ripCalibrate(belt, belt[-1], sr = 10),
    "equal-length"
  )
  expect_error(
    ripCalibrate(belt, rep(0, 100), sr = 10, window = c(0, 10)),
    "window"
  )
})


test_that("flat respiration returns typed empty results", {
  experiment <- make_rip(a_rc = 0)
  onsets <- respiratoryOnsets(experiment, channel = "RC")
  rate_peaks <- breathingRate(experiment, channel = "RC")
  rate_spectral <- breathingRate(
    experiment, channel = "RC", method = "spectral"
  )
  metrics <- tidalMetrics(experiment, channel = "RC", onsets = onsets)

  expect_identical(names(onsets), c(
    "channel", "breath", "onset_sample", "peak_sample", "offset_sample",
    "onset_time", "peak_time", "offset_time", "insp_amplitude"
  ))
  expect_equal(nrow(onsets), 0L)
  expect_true(is.na(rate_peaks$rate_bpm))
  expect_equal(rate_peaks$n_breaths, 0L)
  expect_true(is.na(rate_spectral$rate_bpm))
  expect_equal(nrow(metrics$breaths), 0L)
  expect_equal(metrics$summary$n_breaths, 0L)
  expect_true(is.na(metrics$summary$minute_ventilation))
})


test_that("respiration input resolution validates labels and sampling", {
  experiment <- make_rip()
  expect_error(
    respiratoryOnsets(experiment, channel = "missing"),
    "Unknown"
  )
  expect_error(
    respiratoryOnsets(sin(1:100), sr = NULL),
    "sr"
  )
  expect_error(
    breathingRate(sin(1:100), sr = 10, band = c(0.1, 6)),
    "Nyquist"
  )
  expect_error(
    tidalMetrics(experiment, calibration = 0),
    "calibration"
  )
})


test_that("ventilatoryThreshold finds a planted VT1", {
  vo2 <- seq(1, 3, by = 0.05)
  knot <- 2
  vco2 <- ifelse(
    vo2 <= knot,
    0.9 * vo2,
    0.9 * knot + 1.3 * (vo2 - knot)
  )
  ve <- 25 * vco2
  threshold <- ventilatoryThreshold(vo2, vco2, ve)

  expect_s3_class(threshold, "ventilatory_threshold")
  expect_equal(threshold$vt1$vo2, 2, tolerance = 1e-6)
  expect_equal(threshold$vt1$slope1, 0.9, tolerance = 1e-8)
  expect_equal(threshold$vt1$slope2, 1.3, tolerance = 1e-8)
  expect_true(is.na(threshold$vt2$vo2))
  expect_output(print(threshold), "VT1")
})


test_that("ventilatoryThreshold finds planted respiratory compensation", {
  vo2 <- seq(1, 4, by = 0.05)
  vco2 <- ifelse(
    vo2 <= 2,
    0.9 * vo2,
    0.9 * 2 + 1.3 * (vo2 - 2)
  )
  rcp_vco2 <- vco2[which.min(abs(vo2 - 3))]
  ve <- ifelse(
    vco2 <= rcp_vco2,
    22 * vco2,
    22 * rcp_vco2 + 34 * (vco2 - rcp_vco2)
  )
  threshold <- ventilatoryThreshold(
    data.frame(
      time = seq_along(vo2) * 10,
      vo2 = rev(vo2),
      vco2 = rev(vco2),
      ve = rev(ve)
    )
  )

  expect_equal(threshold$vt1$vo2, 2, tolerance = 1e-6)
  expect_equal(threshold$vt2$vo2, 3, tolerance = 0.1)
  expect_equal(threshold$vt2$vco2, rcp_vco2, tolerance = 0.1)
  expect_gt(threshold$vt2$slope2, threshold$vt2$slope1)
  expect_true(is.finite(threshold$vt2$time))
})


test_that("ventilatoryThreshold is stable to noise and block averaging", {
  set.seed(1)
  vo2 <- seq(1, 3, by = 0.05)
  vco2 <- ifelse(
    vo2 <= 2,
    0.9 * vo2,
    1.8 + 1.3 * (vo2 - 2)
  ) + rnorm(length(vo2), 0, 0.02)
  ve <- 25 * vco2

  noisy <- ventilatoryThreshold(vo2, vco2, ve)
  averaged <- ventilatoryThreshold(
    vo2, vco2, ve, average_n = 2, min_segment = 3
  )
  expect_equal(noisy$vt1$vo2, 2, tolerance = 0.15)
  expect_equal(averaged$vt1$vo2, 2, tolerance = 0.2)
  expect_equal(nrow(averaged$data), ceiling(length(vo2) / 2))
})


test_that("ventilatoryThreshold reports absent and invalid thresholds", {
  vo2 <- seq(1, 2, length.out = 20)
  vco2 <- 0.9 * vo2
  ve <- 20 * vco2
  absent <- ventilatoryThreshold(vo2, vco2, ve)
  short <- ventilatoryThreshold(vo2[1:5], vco2[1:5], ve[1:5])

  expect_true(is.na(absent$vt1$vo2))
  expect_true(is.na(absent$vt2$vo2))
  expect_true(is.na(short$vt1$vo2))
  expect_error(
    ventilatoryThreshold(vo2, vco2[-1], ve),
    "equal-length"
  )
  expect_error(
    ventilatoryThreshold(vo2, vco2, ve, min_segment = 1),
    "at least 2"
  )
})
