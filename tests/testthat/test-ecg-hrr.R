library(testthat)
library(PhysioECG)

test_that("HRR equals the constructed HR drop within 1 bpm", {
  t <- seq(0, 120, by = 1)
  hr <- 170 - (170 - 120) * t / 120           # linear: HR(60) = 145
  res <- ecgHRR(hr, time_sec = t, peak_time = 0)

  expect_type(res, "list")
  expect_equal(res$peak_hr, 170)
  expect_equal(res$hrr1, 25, tolerance = 1)    # 170 - 145
  expect_equal(res$hrr2, 50, tolerance = 1)    # 170 - 120
  expect_equal(res$hrr$offset_sec, c(60, 120))
  expect_lt(res$slope_bpm_per_min, 0)          # HR falls during recovery
  expect_false(res$abnormal)
})

test_that("the abnormal-recovery flag triggers at the threshold", {
  t <- seq(0, 120, by = 1)
  slow <- 170 - (170 - 150) * t / 120          # HR(60) = 160 -> HRR60 = 10
  fast <- 170 - (170 - 120) * t / 120          # HRR60 = 25

  expect_true(ecgHRR(slow, time_sec = t, peak_time = 0)$abnormal)   # 10 < 12
  expect_false(ecgHRR(fast, time_sec = t, peak_time = 0)$abnormal)  # 25 >= 12
  # threshold is configurable
  expect_true(ecgHRR(fast, time_sec = t, peak_time = 0,
                     abnormal_bpm = 30)$abnormal)
})

test_that("HR and RR inputs give the same HRR", {
  t <- seq(0, 120, by = 1)
  hr <- 170 - (170 - 120) * t / 120
  rr <- 60000 / hr                             # equivalent RR (ms)

  r_hr <- ecgHRR(hr, time_sec = t, peak_time = 0, input = "hr")
  r_rr <- ecgHRR(rr, time_sec = t, peak_time = 0, input = "rr")
  r_auto <- ecgHRR(rr, time_sec = t, peak_time = 0)   # auto -> rr by magnitude

  expect_equal(r_rr$hrr1, r_hr$hrr1, tolerance = 1e-6)
  expect_equal(r_auto$hrr1, r_hr$hrr1, tolerance = 1e-6)
})

test_that("data.frame input and auto peak detection work", {
  t <- seq(0, 120, by = 1)
  hr <- 170 - (170 - 120) * t / 120
  res <- ecgHRR(data.frame(time_sec = t, hr = hr))   # auto peak = max HR
  expect_equal(res$peak_time, 0)
  expect_equal(res$peak_hr, 170)
  expect_equal(res$hrr1, 25, tolerance = 1)

  # rr_ms column is recognized
  res_rr <- ecgHRR(data.frame(time_sec = t, rr_ms = 60000 / hr), peak_time = 0)
  expect_equal(res_rr$hrr1, 25, tolerance = 1)
})

test_that("auto peak detection finds the peak of an exercise-recovery ramp", {
  # HR ramps 90 -> 170 (exercise) then recovers 170 -> 120
  t <- seq(0, 180, by = 1)
  hr <- ifelse(t <= 60, 90 + (170 - 90) * t / 60,
               170 - (170 - 120) * (t - 60) / 120)
  res <- ecgHRR(hr, time_sec = t)              # peak auto-detected at t = 60
  expect_equal(res$peak_time, 60)
  expect_equal(res$peak_hr, 170)
  expect_equal(res$hrr1, 25, tolerance = 1)    # HR(120) = 145
})

test_that("offsets beyond the recording yield NA", {
  t <- seq(0, 30, by = 1)                      # only 30 s of recovery
  hr <- 170 - t
  res <- ecgHRR(hr, time_sec = t, peak_time = 0,
                recovery_offsets_sec = c(60, 120))
  expect_true(all(is.na(res$hrr$hrr_bpm)))     # 60 s / 120 s out of range
  expect_error(ecgHRR(1, time_sec = 1), "at least two")
})
