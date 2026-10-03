library(testthat)
library(PhysioECG)

ecg1 <- function(sig, sr = 500) {
  PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                   colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
                   samplingRate = sr)
}
add_qrs <- function(sig, pos, sigma, amp) {
  n <- length(sig)
  tw <- seq(-4 * sigma, 4 * sigma)
  idx <- round(pos) + tw
  ok <- idx >= 1 & idx <= n
  sig[idx[ok]] <- sig[idx[ok]] + amp * exp(-(tw[ok])^2 / (2 * sigma^2))
  sig
}

# Regular sinus + injected wide/premature PVCs with a full compensatory pause.
make_pvc_pe <- function(sr = 500, n_beats = 25, rr_ms = 800,
                        pvc_at = c(6, 12, 18), seed = 1) {
  set.seed(seed)
  RR <- as.integer(rr_ms / 1000 * sr)
  margin <- RR
  g <- (1:n_beats) * RR + margin
  n <- max(g) + RR
  sig <- rnorm(n, sd = 0.02)
  truth <- integer(0)
  for (k in seq_len(n_beats)) {
    if (k %in% pvc_at && k > 1) {
      pvc <- g[k - 1] + round(0.6 * RR)          # premature; g[k] blocked
      sig <- add_qrs(sig, pvc, round(0.035 * sr), -2.2)  # wide, inverted
      truth <- c(truth, pvc)
    } else {
      sig <- add_qrs(sig, g[k], round(0.007 * sr), 1.5)  # narrow, normal
    }
  }
  list(pe = ecg1(sig, sr), truth = truth)
}

test_that("AF detector separates AF from sinus (ROC AUC > 0.9)", {
  set.seed(1)
  mk <- function(v) data.frame(channel = 1L, rr_ms = v)
  af_s <- sinus_s <- numeric(30)
  for (i in seq_len(30)) {
    af <- runif(80, 400, 1200)                        # irregular-irregular
    resp <- 800 + 30 * sin(2 * pi * seq_len(80) / 8)  # respiratory sinus
    sinus <- resp + rnorm(80, sd = 8)
    af_s[i] <- ecgDetectAF(mk(af))$af_prob
    sinus_s[i] <- ecgDetectAF(mk(sinus))$af_prob
  }
  gt <- outer(af_s, sinus_s, ">")
  auc <- (sum(gt) + 0.5 * sum(outer(af_s, sinus_s, "=="))) / (30 * 30)
  expect_gt(auc, 0.9)
  expect_gt(mean(af_s), mean(sinus_s))
})

test_that("ecgDetectAF keeps the per-channel verdict and adds af_prob", {
  set.seed(2)
  af <- data.frame(channel = 1L, rr_ms = runif(80, 400, 1200))
  sinus <- data.frame(channel = 1L,
                      rr_ms = 800 + 20 * sin(2 * pi * seq_len(80) / 8))
  res_af <- ecgDetectAF(af)
  expect_true(all(c("rhythm", "ectopic_frac", "rmssd_ratio", "drr_entropy",
                    "af_prob", "n_beats") %in% names(res_af)))
  expect_equal(res_af$rhythm, "AF")
  expect_gt(res_af$af_prob, ecgDetectAF(sinus)$af_prob)
})

test_that("ecgDetectAF window mode returns per-window AF probability", {
  set.seed(3)
  rr <- data.frame(channel = 1L,
                   rr_ms = c(800 + rnorm(60, sd = 12), runif(60, 400, 1200)))
  w <- ecgDetectAF(rr, window_beats = 30)
  expect_true(all(c("channel", "window", "af_prob", "sample_entropy",
                    "sd1_sd2") %in% names(w)))
  expect_equal(nrow(w), 4L)
  # the AF half scores higher than the sinus half
  expect_gt(mean(w$af_prob[3:4]), mean(w$af_prob[1:2]))
})

test_that("PVC detector flags injected beats (recall > 0.9, FP < 0.1)", {
  for (seed in 1:5) {
    f <- make_pvc_pe(seed = seed)
    pk <- ecgDetectRpeaks(f$pe)
    pv <- ecgDetectPVC(f$pe, pk, ecgRRintervals(f$pe, pk))

    pvc_beats <- vapply(f$truth, function(t) which.min(abs(pv$sample - t)),
                        integer(1))
    recall <- mean(pv$is_pvc[pvc_beats])
    n_normal <- nrow(pv) - length(pvc_beats)
    fp <- sum(pv$is_pvc & !(seq_len(nrow(pv)) %in% pvc_beats)) / n_normal
    expect_gt(recall, 0.9)
    expect_lt(fp, 0.1)
    # PVC beats are wide (low morphology corr) and premature
    expect_true(all(pv$morph_corr[pvc_beats] < 0.9))
    expect_true(all(pv$prematurity[pvc_beats] < 0.9))
  }
})

test_that("compensatory-pause logic distinguishes a PVC from a non-compensatory APB", {
  # PVC: full compensatory pause (comp ~ 1.0)
  f <- make_pvc_pe(seed = 1)
  pk <- ecgDetectRpeaks(f$pe)
  pv <- ecgDetectPVC(f$pe, pk, ecgRRintervals(f$pe, pk))
  pvc_beats <- vapply(f$truth, function(t) which.min(abs(pv$sample - t)),
                      integer(1))
  expect_true(all(pv$compensatory[pvc_beats] > 0.95))

  # APB: premature NORMAL beat, SA node resets -> non-compensatory pause
  sr <- 500; RR <- 400; margin <- 400; nb <- 25; apb_at <- 12
  set.seed(1)
  n <- (nb + 3) * RR
  sig <- rnorm(n, sd = 0.02)
  pos <- margin + RR
  pks <- integer(0); apb_truth <- NA_integer_
  for (k in seq_len(nb)) {
    if (k == apb_at) {
      apb <- pos - round(0.4 * RR); pks <- c(pks, apb); apb_truth <- apb
      pos <- apb + RR                                   # rhythm reset
    } else {
      pks <- c(pks, pos); pos <- pos + RR
    }
  }
  for (p in pks) sig <- add_qrs(sig, p, 3, 1.5)         # all normal morphology
  pe <- ecg1(sig, sr)
  pv2 <- ecgDetectPVC(pe, ecgDetectRpeaks(pe),
                      ecgRRintervals(pe, ecgDetectRpeaks(pe)))
  ai <- which.min(abs(pv2$sample - apb_truth))
  expect_lt(pv2$prematurity[ai], 0.85)                  # premature
  expect_lt(pv2$compensatory[ai], 0.9)                  # NOT compensatory
  expect_gt(pv2$morph_corr[ai], 0.9)                    # normal morphology
  expect_false(pv2$is_pvc[ai])                          # not a PVC
  expect_equal(sum(pv2$is_pvc), 0L)
})
