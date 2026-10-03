library(testthat)
library(PhysioECG)

# Deceleration-dominant RR: a sharp rise (few beats) then a gradual fall (many).
decel_dominant_rr <- function(cycles = 30, seed = 1) {
  set.seed(seed)
  base <- as.numeric(sapply(seq_len(cycles), function(i)
    c(seq(800, 860, length.out = 4), seq(855, 805, length.out = 16))))
  base + rnorm(length(base), sd = 1.5)
}

test_that("HRT gives negative TO and positive TS on the canonical pattern", {
  # beats 1..5 sinus(800), beat 7 is a PVC: rr[6]=coupling(500), rr[7]=pause(1100)
  # post-PVC: early acceleration (780, 785) then steady deceleration
  rr <- c(rep(800, 5), 500, 1100, 780, 785, 790 + 5 * (0:12))
  h <- ecgHRT(rr, pvc_index = 7)

  expect_type(h, "list")
  expect_lt(h$to, 0)                            # acceleration onset
  expect_gt(h$ts, 0)                            # deceleration slope
  expect_equal(h$n_pvc, 1L)
  expect_length(h$tachogram, 15L)
})

test_that("HRT averages multiple PVC tachograms", {
  block <- c(rep(800, 5), 500, 1100, 780, 785, 790 + 5 * (0:12))
  rr <- c(block, block)
  p2 <- length(block) + 7L
  h <- ecgHRT(rr, pvc_index = c(7L, p2))
  expect_equal(h$n_pvc, 2L)
  expect_lt(h$to, 0)
  expect_gt(h$ts, 0)

  # PVCs too close to the edges are skipped
  expect_equal(ecgHRT(rr, pvc_index = 2L)$n_pvc, 0L)
})

test_that("PRSA gives DC > 0 and AC < 0 with dominant decelerations", {
  rr <- decel_dominant_rr()
  p <- ecgPRSA(rr)

  expect_gt(p$dc, 0)
  expect_lt(p$ac, 0)
  expect_gt(abs(p$dc), abs(p$ac))              # deceleration-dominant
  expect_gt(p$n_dc, 0)
  expect_gt(p$n_ac, 0)
})

test_that("PRSA sign flips when the series is time-reversed", {
  rr <- decel_dominant_rr()
  p <- ecgPRSA(rr)
  pr <- ecgPRSA(rev(rr))

  # decelerations become accelerations under time reversal
  expect_gt(pr$dc, 0)
  expect_lt(pr$ac, 0)
  expect_lt(abs(pr$dc), abs(pr$ac))            # now acceleration-dominant
  # DC(reversed) ~ -AC(original) and AC(reversed) ~ -DC(original)
  expect_equal(pr$dc, -p$ac, tolerance = 0.05)
  expect_equal(pr$ac, -p$dc, tolerance = 0.05)
})

test_that("PRSA anchor selection is reproducible and accepts a data.frame", {
  rr <- decel_dominant_rr()
  a <- ecgPRSA(rr)
  b <- ecgPRSA(rr)
  expect_identical(a$dc, b$dc)
  expect_identical(a$ac, b$ac)
  expect_identical(a$n_dc, b$n_dc)
  expect_identical(a$n_ac, b$n_ac)

  df <- data.frame(rr_ms = rr, time_sec = cumsum(rr) / 1000)
  expect_equal(ecgPRSA(df)$dc, a$dc, tolerance = 1e-8)

  # larger wavelet scale still yields DC > 0 / AC < 0
  s4 <- ecgPRSA(rr, s = 4)
  expect_gt(s4$dc, 0)
  expect_lt(s4$ac, 0)
})
