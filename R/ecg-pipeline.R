#' One-Call ECG Processing Pipeline
#'
#' Runs the standard ECG preprocessing chain end to end -- baseline correction,
#' R-peak detection, RR extraction, ectopic-beat correction and (optionally) a
#' bSQI signal-quality gate -- and stores the peaks, RR series and provenance in
#' the object's \code{metadata()} (after the neurokit2 \code{ecg_process()}
#' pattern).
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param detector R-peak detector passed to \code{\link{ecgDetectRpeaks}}
#'   (default "pan_tompkins").
#' @param correct Ectopic-beat correction passed to \code{\link{ecgRRcorrect}}:
#'   "cubic_spline" (default), "interpolate" or "remove".
#' @param sqi_gate If TRUE (default), also compute the bSQI quality gate
#'   (\code{\link{ecgBeatSQIgate}}).
#' @param baseline If TRUE (default), baseline-correct with
#'   \code{\link{ecgBaselineCorrect}} and detect on the corrected assay.
#' @param assay_name Input assay name (default: first assay).
#' @return The (baseline-corrected) PhysioExperiment with an \code{ecg} entry in
#'   \code{metadata()} holding \code{peaks}, \code{rr} (with ectopic flags),
#'   \code{rr_corrected}, \code{bsqi} and the detect-assay name, plus a
#'   provenance record appended via \code{\link[PhysioCore]{withProvenance}}.
#' @seealso \code{\link{ecgAnalyze}}, \code{\link{ecgDetectRpeaks}},
#'   \code{\link{ecgRRcorrect}}, \code{\link{ecgBeatSQIgate}}
#' @export
#' @examples
#' pe <- make_ecg(n_time = 20000, sr = 500)
#' proc <- ecgProcess(pe)
#' names(S4Vectors::metadata(proc)$ecg)
ecgProcess <- function(x, detector = "pan_tompkins", correct = "cubic_spline",
                       sqi_gate = TRUE, baseline = TRUE, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)

  proc <- x
  detect_assay <- assay_name
  if (baseline) {
    proc <- ecgBaselineCorrect(x, assay_name = assay_name)
    detect_assay <- "baseline_corrected"
  }

  peaks <- ecgDetectRpeaks(proc, method = detector, assay_name = detect_assay)
  rr <- ecgRRintervals(proc, peaks)
  rr <- ecgQualityCheck(rr)
  rr_corrected <- ecgRRcorrect(rr, method = correct)
  bsqi <- if (sqi_gate) ecgBeatSQIgate(proc) else NULL

  md <- S4Vectors::metadata(proc)
  md$ecg <- list(detect_assay = detect_assay, detector = detector,
                 correct = correct, peaks = peaks, rr = rr,
                 rr_corrected = rr_corrected, bsqi = bsqi)
  S4Vectors::metadata(proc) <- md

  proc <- PhysioCore::withProvenance(
    x, proc, step = "ecgProcess",
    params = list(detector = detector, correct = correct,
                  baseline = baseline, sqi_gate = sqi_gate))
  proc
}

# Build a one-row-per-channel NA frame with the given numeric columns.
.ecg_na_frame <- function(channels, cols) {
  df <- data.frame(channel = channels)
  for (c in cols) df[[c]] <- NA_real_
  df
}

#' Aggregate ECG HRV and Morphology Metrics
#'
#' Aggregates the full HRV panel -- time-domain, frequency-domain (AR spectrum
#' with normalized units), geometric and nonlinear -- together with a QTc
#' morphology summary into a single per-channel data.frame, from an object
#' produced by \code{\link{ecgProcess}} (after the neurokit2
#' \code{ecg_analyze()} pattern).
#'
#' @param processed A PhysioExperiment returned by \code{\link{ecgProcess}}.
#' @return A data.frame with one row per channel and columns for the time-domain
#'   (\code{sdnn}, \code{rmssd}, \code{pnn50}, \code{mean_hr}), frequency-domain
#'   (\code{lf}, \code{hf}, \code{lf_hf_ratio}, \code{lf_nu}, \code{hf_nu},
#'   \code{total_power}), geometric (\code{HRV_triangular_index}, \code{TINN}),
#'   nonlinear (\code{sd1}, \code{sd2}, \code{sample_entropy}, \code{alpha1},
#'   \code{alpha2}) and morphology (\code{qt_ms}, \code{qtc_bazett},
#'   \code{qrs_ms}) metrics. Metrics that cannot be computed (e.g. DFA on a short
#'   recording) are \code{NA}.
#' @seealso \code{\link{ecgProcess}}, \code{\link{ecgHRVtime}},
#'   \code{\link{ecgHRVfreq}}, \code{\link{ecgHRVgeometric}},
#'   \code{\link{ecgHRVnonlinear}}, \code{\link{ecgIntervals}}
#' @export
#' @examples
#' pe <- make_ecg(n_time = 40000, sr = 500)
#' ecgAnalyze(ecgProcess(pe))
ecgAnalyze <- function(processed) {
  stopifnot(inherits(processed, "PhysioExperiment"))
  ecg <- S4Vectors::metadata(processed)$ecg
  if (is.null(ecg)) {
    stop("'processed' has no ecg metadata; run ecgProcess() first",
         call. = FALSE)
  }
  rr <- ecg$rr_corrected
  peaks <- ecg$peaks
  sr <- samplingRate(processed)
  channels <- sort(unique(rr$channel))

  tm <- ecgHRVtime(rr, rhythm_check = FALSE)[,
    c("channel", "sdnn", "rmssd", "pnn50", "mean_hr")]
  fr <- ecgHRVfreq(rr, method = "ar", rhythm_check = FALSE)[,
    c("channel", "lf", "hf", "lf_hf_ratio", "lf_nu", "hf_nu", "total_power")]

  geo <- do.call(rbind, lapply(channels, function(ch) {
    g <- tryCatch(ecgHRVgeometric(rr$rr_ms[rr$channel == ch]),
                  error = function(e) NULL)
    data.frame(channel = ch,
               HRV_triangular_index = if (is.null(g)) NA_real_ else
                 g$HRV_triangular_index,
               TINN = if (is.null(g)) NA_real_ else g$TINN)
  }))

  nl <- tryCatch(
    ecgHRVnonlinear(rr, rhythm_check = FALSE)[,
      c("channel", "sd1", "sd2", "sample_entropy", "alpha1", "alpha2")],
    error = function(e) .ecg_na_frame(channels,
      c("sd1", "sd2", "sample_entropy", "alpha1", "alpha2")))

  morph <- tryCatch({
    del <- ecgDelineate(processed, peaks, assay_name = ecg$detect_assay)
    iv <- ecgIntervals(del, sr)
    do.call(rbind, lapply(channels, function(ch) {
      d <- iv[iv$channel == ch, , drop = FALSE]
      data.frame(channel = ch,
                 qt_ms = stats::median(d$qt_ms, na.rm = TRUE),
                 qtc_bazett = stats::median(d$qtc_bazett, na.rm = TRUE),
                 qrs_ms = stats::median(d$qrs_ms, na.rm = TRUE))
    }))
  }, error = function(e)
    .ecg_na_frame(channels, c("qt_ms", "qtc_bazett", "qrs_ms")))

  Reduce(function(a, b) merge(a, b, by = "channel", sort = TRUE),
         list(tm, fr, geo, nl, morph))
}
