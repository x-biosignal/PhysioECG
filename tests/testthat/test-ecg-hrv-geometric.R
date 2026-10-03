test_that("HRV triangular index equals N / modal count (hand-checked histogram)", {
  # bin_ms = 10: values fall in bins (780,790],(790,800],(800,810]
  #   790 x10 -> bin1, 800 x60 -> bin2 (mode), 810 x30 -> bin3
  rr <- c(rep(800, 60), rep(810, 30), rep(790, 10))
  g <- ecgHRVgeometric(rr, bin_ms = 10)
  expect_equal(g$HRV_triangular_index, length(rr) / 60, tolerance = 1e-6)
  expect_gte(g$TINN, 0)
})

test_that("TINN is non-negative and scales with NN dispersion", {
  set.seed(1)
  narrow <- 800 + rnorm(2000, sd = 10)
  wide   <- 800 + rnorm(2000, sd = 60)
  tn <- ecgHRVgeometric(narrow)$TINN
  tw <- ecgHRVgeometric(wide)$TINN
  expect_gte(tn, 0)
  expect_gt(tw, tn)
})

test_that("ecgHRVgeometric handles degenerate input", {
  expect_true(is.na(ecgHRVgeometric(numeric(0))$TINN))
  expect_true(is.na(ecgHRVgeometric(c(800))$HRV_triangular_index))
})

test_that("ecgHRVsegmented: SDANN=0 for equal segment means, SDNN_index~SDNN when stationary", {
  set.seed(2)
  seg <- function() 1000 + rnorm(300, sd = 30)  # ~5 min at 1000 ms/beat
  rr <- c(seg(), seg(), seg())
  res <- ecgHRVsegmented(rr, segment_sec = 300)
  expect_equal(res$n_segments, 3L)              # short trailing segment dropped
  # stationary series: mean of per-segment SDNN approximates the whole SDNN
  expect_equal(res$SDNN_index, stats::sd(rr), tolerance = 0.25 * stats::sd(rr))

  # every complete 5-min segment has identical mean -> SDANN == 0
  const_seg <- rep(1000, 300)
  res2 <- ecgHRVsegmented(c(const_seg, const_seg, const_seg), segment_sec = 300)
  expect_equal(res2$SDANN, 0, tolerance = 1e-9)

  # off-multiple recording: a short tail must NOT inflate SDANN
  rr3 <- c(rep(1000, 600), rep(600, 50))        # 600 s of 1000 ms + 30 s tail
  res3 <- ecgHRVsegmented(rr3, segment_sec = 300)
  expect_equal(res3$SDANN, 0, tolerance = 1e-9)
  expect_equal(res3$n_segments, 2L)
})
