make_af <- function(seed, n = 300) {
  set.seed(seed)
  data.frame(channel = 1L, rr_ms = runif(n, 480, 920),
             time_sec = cumsum(runif(n, 0.48, 0.92)))
}
make_sinus <- function(seed, n = 300) {
  set.seed(seed)
  data.frame(channel = 1L, rr_ms = 800 + rnorm(n, sd = 25),
             time_sec = cumsum(rep(0.8, n)))
}

test_that("ecgDetectAF distinguishes AF from sinus", {
  res <- ecgDetectAF(make_af(1))
  expect_true(all(c("rhythm", "ectopic_frac", "rmssd_ratio", "drr_entropy",
                    "n_beats") %in% names(res)))
  expect_equal(res$rhythm, "AF")
  expect_gte(res$rmssd_ratio, 0.15)
  expect_gte(res$drr_entropy, 0.70)

  expect_equal(ecgDetectAF(make_sinus(1))$rhythm, "sinus")
})

test_that("ecgHRVtime gates AF: metrics NA, hrv_valid FALSE, one warning", {
  expect_warning(res <- ecgHRVtime(make_af(2)), "atrial fibrillation")
  expect_equal(res$rhythm, "AF")
  expect_false(res$hrv_valid)
  expect_true(is.na(res$sdnn) && is.na(res$rmssd))
})

test_that("rhythm_check leaves sinus HRV values byte-identical (backward compatible)", {
  sinus <- make_sinus(3)
  gated <- ecgHRVtime(sinus, rhythm_check = TRUE)
  plain <- ecgHRVtime(sinus, rhythm_check = FALSE)
  expect_true(gated$hrv_valid)
  expect_equal(gated$sdnn, plain$sdnn)
  expect_equal(gated$rmssd, plain$rmssd)
  expect_false("rhythm" %in% names(plain))     # gate off -> no extra columns
})

test_that("AF gate also applies to ecgHRVfreq and ecgHRVnonlinear", {
  af <- make_af(4)
  expect_warning(f <- ecgHRVfreq(af), "atrial fibrillation")
  expect_false(f$hrv_valid)
  expect_true(is.na(f$lf))
  expect_warning(n <- ecgHRVnonlinear(af), "atrial fibrillation")
  expect_false(n$hrv_valid)
})

test_that("ecgDetectAF honors its documented channel + rr_ms contract (no time_sec)", {
  set.seed(1)
  rr <- data.frame(channel = 1L, rr_ms = 800 + rnorm(50, sd = 20))
  expect_silent(res <- ecgDetectAF(rr))          # was: error on missing time_sec
  expect_equal(res$rhythm, "sinus")
})

test_that("ecgDetectAF tolerates NA in rr_ms rather than crashing", {
  set.seed(1)
  rr <- data.frame(channel = 1L,
                   rr_ms = c(800 + rnorm(49, sd = 20), NA_real_),
                   time_sec = cumsum(rep(0.8, 50)))
  expect_error(ecgDetectAF(rr), NA)              # completes without error
})

test_that(".drr_shannon_entropy rejects n_bins < 2 instead of returning NaN", {
  ent <- PhysioECG:::.drr_shannon_entropy
  expect_error(ent(rnorm(50, sd = 100), n_bins = 1L), "n_bins")
  expect_false(is.nan(ent(rnorm(50, sd = 100), n_bins = 16L)))
})

test_that("undetermined rhythm (too few beats) is hrv_valid = NA, metrics retained", {
  set.seed(5)
  short <- data.frame(channel = 1L, rr_ms = 800 + rnorm(12, sd = 20),
                      time_sec = cumsum(rep(0.8, 12)))          # < min_beats (20)
  res <- ecgHRVtime(short)
  expect_true(is.na(res$rhythm))
  expect_true(is.na(res$hrv_valid))              # unknown, not silently FALSE
  expect_false(is.na(res$sdnn))                  # metrics kept
})


# --- PhysioNet afdb golden: real labelled AF-vs-sinus RR (WSCB-05 / WSCB-09) ---

test_that("ecgDetectAF discriminates AF from sinus on a real PhysioNet record", {
  fx <- readRDS(test_path("fixtures", "afdb_af_rr.rds"))
  rr <- fx$rr_ms; lab <- fx$label
  window <- 60L
  nw <- length(rr) %/% window
  # Evaluate a deterministic, evenly-spaced subset of windows: ~250 is ample for
  # a stable AUC and keeps the test quick.
  stride <- max(1L, nw %/% 250L)
  ks <- seq(1L, nw, by = stride)

  wins <- lapply(ks, function(k) {
    ix <- ((k - 1L) * window + 1L):(k * window)
    list(rr = rr[ix],
         af_prob = ecgDetectAF(data.frame(channel = 1L, rr_ms = rr[ix]))$af_prob,
         label = names(sort(table(lab[ix]), decreasing = TRUE))[1])
  })
  wins <- Filter(function(w) is.finite(w$af_prob) && !is.na(w$label), wins)
  ap <- vapply(wins, `[[`, numeric(1), "af_prob")
  wl <- vapply(wins, `[[`, character(1), "label")
  af <- ap[wl == "AF"]; sn <- ap[wl == "sinus"]
  expect_gt(length(af), 20)
  expect_gt(length(sn), 20)

  # the continuous AF probability separates the two rhythms (AUC via rank sum)
  auc <- as.numeric(stats::wilcox.test(af, sn)$statistic) / (length(af) * length(sn))
  expect_gt(auc, 0.90)
  expect_gt(mean(af) - mean(sn), 0.15)

  # the binary gate fires correctly on the clearest AF window (taken from the
  # window itself, not a reconstructed index)
  clearest_af <- wins[wl == "AF"][[which.max(af)]]
  expect_equal(ecgDetectAF(data.frame(channel = 1L,
                                      rr_ms = clearest_af$rr))$rhythm, "AF")
})
