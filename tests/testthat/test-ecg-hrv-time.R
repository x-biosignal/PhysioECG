library(testthat)
library(PhysioECG)

test_that("ecgHRVtime computes correct metrics from known RR intervals", {
  # Known RR intervals: 800, 850, 900, 850, 800, 850, 900, 850, 800, 850 ms
  rr <- data.frame(
    channel = rep(1L, 10),
    rr_ms = c(800, 850, 900, 850, 800, 850, 900, 850, 800, 850),
    time_sec = cumsum(c(0, rep(0.85, 9)))
  )

  result <- ecgHRVtime(rr)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$mean_rr, mean(rr$rr_ms))
  expect_equal(result$sdnn, sd(rr$rr_ms), tolerance = 0.01)

  diffs <- diff(rr$rr_ms)
  expected_rmssd <- sqrt(mean(diffs^2))
  expect_equal(result$rmssd, expected_rmssd, tolerance = 0.01)

  expected_pnn50 <- 100 * sum(abs(diffs) > 50) / length(diffs)
  expect_equal(result$pnn50, expected_pnn50)

  expect_equal(result$mean_hr, 60000 / mean(rr$rr_ms), tolerance = 0.01)
})

test_that("ecgHRVtime handles multiple channels", {
  rr <- data.frame(
    channel = c(rep(1L, 5), rep(2L, 5)),
    rr_ms = c(800, 850, 900, 850, 800, 1000, 1050, 1000, 950, 1000),
    time_sec = c(cumsum(c(0, rep(0.85, 4))), cumsum(c(0, rep(1.0, 4))))
  )

  result <- ecgHRVtime(rr)
  expect_equal(nrow(result), 2)
  expect_equal(result$channel, c(1L, 2L))
})

test_that("ecgHRVtime validates input", {
  expect_error(ecgHRVtime("not_df"))
  expect_error(ecgHRVtime(data.frame(x = 1)))
})
