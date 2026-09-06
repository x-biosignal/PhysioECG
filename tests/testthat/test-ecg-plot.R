library(testthat)
library(PhysioECG)

# Force a headless render so build/device errors surface without a device.
built <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  invisible(p)
}

rr_fixture <- function(n_time = 40000, sr = 500) {
  pe <- make_ecg(n_time = n_time, sr = sr)
  ecgQualityCheck(ecgRRintervals(pe, ecgDetectRpeaks(pe)))
}

test_that("plotTachogram renders and accepts an ecgProcess object", {
  rr <- rr_fixture()
  built(plotTachogram(rr))
  built(plotTachogram(ecgProcess(make_ecg(n_time = 30000, sr = 500))))
})

test_that("plotPoincare renders and its ellipse SD1/SD2 equal ecgHRVpoincare", {
  rr <- rr_fixture()
  p <- built(plotPoincare(rr))
  pc <- ecgHRVpoincare(rr)
  expect_equal(unname(attr(p, "sd1")[1]), pc$sd1, tolerance = 1e-8)
  expect_equal(unname(attr(p, "sd2")[1]), pc$sd2, tolerance = 1e-8)
})

test_that("plotHRVpsd band shading integrates to the ecgHRVfreq band powers", {
  rr <- rr_fixture()
  p <- built(plotHRVpsd(rr, method = "welch"))
  bp <- attr(p, "band_power")
  fr <- ecgHRVfreq(rr, method = "welch", rhythm_check = FALSE)

  get <- function(b) bp$power[bp$band == b]
  expect_lt(abs(get("LF") - fr$lf) / fr$lf, 0.05)
  expect_lt(abs(get("HF") - fr$hf) / fr$hf, 0.05)
  expect_lt(abs(get("VLF") - fr$vlf) / fr$vlf, 0.05)
})

test_that("plotEcgWave renders with P/QRS/T marks and editable peak indices", {
  pe <- make_ecg(n_time = 10000, sr = 500)
  pk <- ecgDetectRpeaks(pe)
  p <- built(plotEcgWave(pe, pk))
  # editable peak indices exposed for beat editing
  expect_equal(attr(p, "peaks"), pk$sample[pk$channel == 1])
  expect_s3_class(attr(p, "delineation"), "data.frame")
  expect_true("r_peak" %in% names(attr(p, "delineation")))
  # window_sec limits the displayed span
  built(plotEcgWave(pe, pk, window_sec = 4))
})

test_that("the ECG plots render on multiple make_ecg* fixtures", {
  for (pe in list(make_ecg(n_time = 30000, sr = 500),
                  make_ecg_noisy(n_time = 30000, sr = 500),
                  make_ecg_pqrst(n_time = 30000, sr = 500)$pe)) {
    pk <- ecgDetectRpeaks(pe)
    rr <- ecgQualityCheck(ecgRRintervals(pe, pk))
    built(plotTachogram(rr))
    built(plotPoincare(rr))
    built(plotHRVpsd(rr))
    built(plotEcgWave(pe, pk))
  }
})
