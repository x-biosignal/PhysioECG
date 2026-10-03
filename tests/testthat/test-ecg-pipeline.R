library(testthat)
library(PhysioECG)

test_that("ecgProcess records peaks, corrected RR, bSQI and provenance", {
  pe <- make_ecg_noisy(n_time = 40000, sr = 500)
  proc <- ecgProcess(pe)

  expect_s4_class(proc, "PhysioExperiment")
  md <- S4Vectors::metadata(proc)$ecg
  expect_false(is.null(md))
  expect_true(all(c("peaks", "rr", "rr_corrected", "bsqi", "detect_assay") %in%
                    names(md)))
  expect_gt(nrow(md$peaks), 0)
  expect_gt(nrow(md$rr_corrected), 0)
  expect_true("bSQI" %in% names(md$bsqi))
  expect_true("is_ectopic" %in% names(md$rr))

  # provenance recorded via PhysioCore
  prov <- PhysioCore::provenance(proc)
  expect_s3_class(prov, "data.frame")
  expect_true(any(prov$step == "ecgProcess"))
})

test_that("ecgAnalyze columns equal the individual functions on the same RR", {
  pe <- make_ecg(n_time = 40000, sr = 500)
  proc <- ecgProcess(pe)
  an <- ecgAnalyze(proc)
  rc <- S4Vectors::metadata(proc)$ecg$rr_corrected

  tm <- ecgHRVtime(rc, rhythm_check = FALSE)
  fr <- ecgHRVfreq(rc, method = "ar", rhythm_check = FALSE)
  nl <- ecgHRVnonlinear(rc, rhythm_check = FALSE)

  expect_equal(an$sdnn, tm$sdnn, tolerance = 1e-6)
  expect_equal(an$rmssd, tm$rmssd, tolerance = 1e-6)
  expect_equal(an$pnn50, tm$pnn50, tolerance = 1e-6)
  expect_equal(an$lf, fr$lf, tolerance = 1e-6)
  expect_equal(an$hf, fr$hf, tolerance = 1e-6)
  expect_equal(an$lf_nu, fr$lf_nu, tolerance = 1e-6)
  expect_equal(an$hf_nu, fr$hf_nu, tolerance = 1e-6)
  expect_equal(an$sd1, nl$sd1, tolerance = 1e-6)
  expect_equal(an$sd2, nl$sd2, tolerance = 1e-6)
  expect_equal(an$alpha1, nl$alpha1, tolerance = 1e-6)
})

test_that("ecgAnalyze returns one row per channel with the full metric panel", {
  proc <- ecgProcess(make_ecg(n_time = 40000, sr = 500))
  an <- ecgAnalyze(proc)
  expect_s3_class(an, "data.frame")
  expect_equal(nrow(an), 1L)
  expect_true(all(c("channel", "sdnn", "rmssd", "lf", "hf", "lf_nu",
                    "HRV_triangular_index", "sd1", "alpha1",
                    "qtc_bazett", "qrs_ms") %in% names(an)))

  expect_error(ecgAnalyze(make_ecg(sr = 500)), "run ecgProcess")
})

test_that("the pipeline runs end-to-end on every make_ecg* fixture", {
  fixtures <- list(
    make_ecg(n_time = 30000, sr = 500),
    make_ecg_noisy(n_time = 30000, sr = 500),
    make_ecg_irregular(n_time = 20000, sr = 500),
    make_ecg_pqrst(n_time = 30000, sr = 500)$pe)
  for (pe in fixtures) {
    proc <- ecgProcess(pe)
    an <- expect_no_error(ecgAnalyze(proc))
    expect_s3_class(an, "data.frame")
    expect_gt(ncol(an), 15)
  }
})

test_that("ecgProcess options (no baseline, no gate, alternative correction) work", {
  pe <- make_ecg(n_time = 30000, sr = 500)
  p1 <- ecgProcess(pe, baseline = FALSE, sqi_gate = FALSE,
                   correct = "interpolate")
  md <- S4Vectors::metadata(p1)$ecg
  expect_equal(md$detect_assay, "raw")
  expect_null(md$bsqi)
  expect_gt(nrow(md$rr_corrected), 0)
})
