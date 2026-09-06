# =============================================================================
# Golden REGRESSION tests for PhysioECG
# =============================================================================
# Each test rebuilds the SAME deterministic input used by data-raw/golden.R,
# computes the metric via the PACKAGE function, and compares against the stored
# INDEPENDENT-reference golden (skips if the fixture is absent).
#
# Golden sources (see data-raw/golden.R):
#   hrv_time_domain -> base-R sd()/sqrt(mean(diff^2))/pNN50
#   dfa_alpha       -> inline standard Peng et al. (1994) DFA (base stats::lm)
# =============================================================================

library(testthat)
library(PhysioECG)

# ---- Deterministic input builders (MUST match data-raw/golden.R) ----

build_rr_hrv <- function() {
  set.seed(2024)
  n <- 256L
  base <- 800 + cumsum(rnorm(n, 0, 4))
  base[seq(10L, n, by = 20L)] <- base[seq(10L, n, by = 20L)] + 70
  rr_ms <- pmax(base, 400)
  data.frame(channel = 1L, rr_ms = rr_ms,
             time_sec = cumsum(rr_ms) / 1000, stringsAsFactors = FALSE)
}

build_rr_dfa <- function() {
  set.seed(99)
  n <- 512L
  rr_ms <- 800 + cumsum(rnorm(n, 0, 3))
  rr_ms <- pmax(rr_ms, 400)
  data.frame(channel = 1L, rr_ms = rr_ms,
             time_sec = cumsum(rr_ms) / 1000, stringsAsFactors = FALSE)
}


test_that("ecgHRVtime matches independent base-R reference (golden)", {
  rr <- build_rr_hrv()
  # rhythm_check = FALSE to test the pure time-domain formulas independently of
  # the AF/ectopy gate (the reference does not depend on it).
  res <- ecgHRVtime(rr, rhythm_check = FALSE)

  actual <- data.frame(
    channel = as.integer(res$channel),
    mean_rr = res$mean_rr,
    sdnn    = res$sdnn,
    rmssd   = res$rmssd,
    pnn50   = res$pnn50,
    mean_hr = res$mean_hr,
    stringsAsFactors = FALSE
  )

  expect_equal_golden(actual, "hrv_time_domain", tol = 1e-8)
})


test_that("ecgDFA alpha1/alpha2 match independent standard-DFA reference (golden)", {
  rr <- build_rr_dfa()
  res <- ecgDFA(rr, short_range = c(4, 16), long_range = c(16, 64))

  actual <- data.frame(
    channel = as.integer(res$channel),
    alpha1  = res$alpha1,
    alpha2  = res$alpha2,
    stringsAsFactors = FALSE
  )

  expect_equal_golden(actual, "dfa_alpha", tol = 1e-8)
})
