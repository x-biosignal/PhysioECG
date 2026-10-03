library(testthat)
library(PhysioECG)

# VAL-05: validate R-peak detection and HRV against REAL MIT-BIH data (mitdb),
# replacing the synthetic-RR substitute in the ECG "real data" validation path.
# Fixture: raw ECG (MLII, 360 Hz) + cardiologist beat annotations for a clean
# record (100) and the noisiest record (203), plus a 5-minute annotation-derived
# reference RR series (record 100). Provenance/generation:
# data-raw/ecg_mitbih_reference.{py,R}.

fx_path <- test_path("fixtures", "mitbih-reference.rds")

# precision/recall/F1 by greedy nearest matching within `tol` samples (AAMI EC57)
.beat_prf <- function(detected, reference, tol) {
  detected <- sort(unique(detected))
  matched <- logical(length(reference))
  tp <- 0
  for (d in detected) {
    i <- which(!matched & abs(reference - d) <= tol)
    if (length(i)) { matched[i[which.min(abs(reference[i] - d))]] <- TRUE; tp <- tp + 1 }
  }
  fp <- length(detected) - tp
  fn <- sum(!matched)
  list(precision = tp / (tp + fp), recall = tp / (tp + fn),
       f1 = 2 * tp / (2 * tp + fp + fn), tp = tp, fp = fp, fn = fn)
}

test_that("ecgDetectRpeaks matches MIT-BIH beat annotations (precision/recall/F1)", {
  skip_if(!file.exists(fx_path), "MIT-BIH reference fixture not bundled")
  fx <- readRDS(fx_path)
  # per-record F1 floors: clean record (100) near-perfect; noisiest record (203)
  # is the hardest in mitdb (muscle artefact, frequent ectopy).
  floors <- c("mitdb/100" = 0.99, "mitdb/203" = 0.90)
  for (d in fx$detection) {
    sig <- matrix(d$signal, ncol = 1)
    colnames(sig) <- d$channel
    pe <- PhysioExperiment(assays = list(raw = sig), samplingRate = d$fs)
    det <- ecgDetectRpeaks(pe, method = "pan_tompkins")$sample
    m <- .beat_prf(det, d$ref_peaks, round(0.15 * d$fs))
    expect_gt(m$f1, floors[[d$record]], label = sprintf("%s F1", d$record))
    expect_gt(m$precision, 0.9, label = sprintf("%s precision", d$record))
    expect_gt(m$recall, 0.9, label = sprintf("%s recall", d$record))
  }
})

test_that("ecgHRVtime reproduces the standard HRV definitions on real MIT-BIH RR", {
  skip_if(!file.exists(fx_path), "MIT-BIH reference fixture not bundled")
  fx <- readRDS(fx_path)
  rr <- fx$hrv_rr
  rr_df <- data.frame(channel = 1L, rr_ms = rr, time_sec = cumsum(rr) / 1000)
  h <- ecgHRVtime(rr_df, rhythm_check = FALSE)
  expect_equal(h$sdnn, sd(rr), tolerance = 1e-6)
  expect_equal(h$rmssd, sqrt(mean(diff(rr)^2)), tolerance = 1e-6)
  expect_equal(h$pnn50, 100 * mean(abs(diff(rr)) > 50), tolerance = 1e-6)
})

test_that("ecgHRVtime agrees with RHRV on real MIT-BIH RR", {
  skip_if(!file.exists(fx_path), "MIT-BIH reference fixture not bundled")
  skip_if_not_installed("RHRV")
  fx <- readRDS(fx_path)
  rr <- fx$hrv_rr
  rr_df <- data.frame(channel = 1L, rr_ms = rr, time_sec = cumsum(rr) / 1000)
  h <- ecgHRVtime(rr_df, rhythm_check = FALSE)
  hd <- RHRV::CreateHRVData()
  hd$Beat <- data.frame(Time = cumsum(c(0, rr)) / 1000)
  hd <- RHRV::BuildNIHR(hd)
  hd <- suppressWarnings(RHRV::CreateTimeAnalysis(hd))
  ta <- hd$TimeAnalysis[[1]]
  # RHRV drops the first NIHR interval, so allow ~1% (SDNN/RMSSD) tolerance.
  expect_lt(abs(h$sdnn - ta$SDNN) / ta$SDNN, 0.01)
  expect_lt(abs(h$rmssd - ta$rMSSD) / ta$rMSSD, 0.01)
})
