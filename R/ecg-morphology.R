#' Delineate ECG Waveform Morphology
#'
#' For each detected R-peak, identifies QRS complex boundaries (onset and
#' offset) and P-wave and T-wave peaks. QRS boundaries are detected using a
#' gradient-based search from the R-peak, while P and T waves are found as
#' local maxima in physiologically plausible time windows.
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param peaks A data.frame of detected R-peaks as returned by
#'   \code{\link{ecgDetectRpeaks}}, with columns \code{channel} and
#'   \code{sample}.
#' @param assay_name Name of the assay to use. If \code{NULL}, the default
#'   assay is used.
#' @return A data.frame with one row per beat and the following columns:
#'   \describe{
#'     \item{channel}{Integer channel index (1-based).}
#'     \item{beat}{Integer beat number within the channel (1-based).}
#'     \item{r_peak}{Sample index of the R-peak.}
#'     \item{qrs_onset}{Sample index of QRS complex onset.}
#'     \item{qrs_offset}{Sample index of QRS complex offset (J-point).}
#'     \item{qrs_duration_ms}{QRS complex duration in milliseconds.}
#'     \item{p_peak}{Sample index of P-wave peak, or \code{NA} if not found
#'       in the search window (300--80 ms before R-peak).}
#'     \item{t_peak}{Sample index of T-wave peak, or \code{NA} if not found
#'       in the search window (80--500 ms after R-peak).}
#'     \item{t_end}{Sample index of T-wave end estimated by
#'       tangent-intercept method, or \code{NA} if T-wave not found.}
#'   }
#'   Returns a zero-row data.frame with the same column structure if no beats
#'   are delineated.
#'
#' @references Goldberger, A.L., et al. (2000). "PhysioBank, PhysioToolkit,
#'   and PhysioNet: Components of a new research resource for complex
#'   physiologic signals." \emph{Circulation}, 101(23), e215--e220.
#'   \doi{10.1161/01.CIR.101.23.e215}
#'
#' @seealso \code{\link{ecgDetectRpeaks}} for R-peak detection,
#'   \code{\link{ecgIntervals}} for computing clinical ECG intervals from
#'   delineation results, \code{\link{ecgSignalQuality}} for signal quality
#'   assessment.
#'
#' @export
#' @examples
#' \dontrun{
#' pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
#' peaks <- ecgDetectRpeaks(pe)
#' delin <- ecgDelineate(pe, peaks)
#' head(delin)
#' }
ecgDelineate <- function(x, peaks, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.data.frame(peaks))
  required <- c("channel", "sample")
  missing_cols <- setdiff(required, names(peaks))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in peaks: ",
         paste(missing_cols, collapse = ", "))
  }

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)

  # Window parameters in samples
  p_window_start_ms <- 300   # search P wave from -300ms before R

  p_window_end_ms <- 80      # to -80ms before R
  t_window_start_ms <- 80    # search T wave from +80ms after R
  t_window_end_ms <- 500     # to +500ms after R

  p_win_start <- as.integer(round(p_window_start_ms / 1000 * sr))
  p_win_end <- as.integer(round(p_window_end_ms / 1000 * sr))
  t_win_start <- as.integer(round(t_window_start_ms / 1000 * sr))
  t_win_end <- as.integer(round(t_window_end_ms / 1000 * sr))

  channels <- unique(peaks$channel)
  results <- list()

  for (ch in channels) {
    ch_peaks <- peaks[peaks$channel == ch, ]
    ch_peaks <- ch_peaks[order(ch_peaks$sample), ]
    sig <- data[, ch]

    for (beat_idx in seq_len(nrow(ch_peaks))) {
      r_pos <- ch_peaks$sample[beat_idx]

      # --- QRS onset: search backward from R-peak ---
      qrs_onset <- .find_qrs_onset(sig, r_pos, sr, n_time)

      # --- QRS offset: search forward from R-peak past S wave ---
      qrs_offset <- .find_qrs_offset(sig, r_pos, sr, n_time)

      qrs_duration_ms <- (qrs_offset - qrs_onset) / sr * 1000

      # --- P wave: local max in [-300ms, -80ms] window before R ---
      p_peak <- .find_wave_peak(sig, r_pos, -p_win_start, -p_win_end,
                                n_time, find_max = TRUE)

      # --- T wave: local max in [+80ms, +500ms] window after R ---
      t_peak <- .find_wave_peak(sig, r_pos, t_win_start, t_win_end,
                                n_time, find_max = TRUE)

      # --- T wave end: search forward from T peak for baseline return ---
      t_end <- if (is.na(t_peak)) {
        NA_integer_
      } else {
        .find_t_end(sig, t_peak, sr, n_time)
      }

      results[[length(results) + 1]] <- data.frame(
        channel = as.integer(ch),
        beat = as.integer(beat_idx),
        r_peak = as.integer(r_pos),
        qrs_onset = as.integer(qrs_onset),
        qrs_offset = as.integer(qrs_offset),
        qrs_duration_ms = qrs_duration_ms,
        p_peak = if (is.na(p_peak)) NA_integer_ else as.integer(p_peak),
        t_peak = if (is.na(t_peak)) NA_integer_ else as.integer(t_peak),
        t_end = if (is.na(t_end)) NA_integer_ else as.integer(t_end),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(
      channel = integer(0), beat = integer(0), r_peak = integer(0),
      qrs_onset = integer(0), qrs_offset = integer(0),
      qrs_duration_ms = numeric(0), p_peak = integer(0),
      t_peak = integer(0), t_end = integer(0)
    )
  }
}


#' Compute ECG Intervals from Delineation
#'
#' Calculates standard clinical ECG intervals from the waveform delineation
#' produced by \code{\link{ecgDelineate}}: PR interval, QT interval, QTc
#' (Bazett correction), QRS duration, and RR interval.
#'
#' @param delineation A data.frame as returned by \code{\link{ecgDelineate}},
#'   with columns \code{channel}, \code{beat}, \code{r_peak}, \code{qrs_onset},
#'   \code{qrs_offset}, \code{p_peak}, \code{t_peak}, \code{t_end}.
#' @param sr Sampling rate in Hz.
#' @return A data.frame with one row per beat and the following columns:
#'   \describe{
#'     \item{channel}{Integer channel index (1-based).}
#'     \item{beat}{Integer beat number within the channel.}
#'     \item{pr_ms}{PR interval in milliseconds (P-wave peak to QRS onset),
#'       or \code{NA} if the P wave was not detected.}
#'     \item{qt_ms}{QT interval in milliseconds (QRS onset to T-wave end),
#'       or \code{NA} if the T wave was not detected.}
#'     \item{qtc_ms}{Corrected QT interval using Bazett's formula
#'       (\code{QT / sqrt(RR_sec)}), or \code{NA} if QT or RR is
#'       unavailable.}
#'     \item{qrs_ms}{QRS complex duration in milliseconds.}
#'     \item{rr_ms}{RR interval in milliseconds to the next beat, or
#'       \code{NA} for the last beat in each channel.}
#'   }
#'   Returns a zero-row data.frame with the same column structure if no beats
#'   are present.
#'
#' @references Goldberger, A.L., et al. (2000). "PhysioBank, PhysioToolkit,
#'   and PhysioNet: Components of a new research resource for complex
#'   physiologic signals." \emph{Circulation}, 101(23), e215--e220.
#'   \doi{10.1161/01.CIR.101.23.e215}
#'
#' @seealso \code{\link{ecgDelineate}} for waveform delineation,
#'   \code{\link{ecgDetectRpeaks}} for R-peak detection,
#'   \code{\link{ecgRRintervals}} for RR interval computation.
#'
#' @export
#' @examples
#' \dontrun{
#' pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
#' peaks <- ecgDetectRpeaks(pe)
#' delin <- ecgDelineate(pe, peaks)
#' intervals <- ecgIntervals(delin, samplingRate(pe))
#' head(intervals)
#' }
ecgIntervals <- function(delineation, sr) {
  stopifnot(is.data.frame(delineation))
  required <- c("channel", "beat", "r_peak", "qrs_onset", "qrs_offset",
                 "p_peak", "t_peak", "t_end")
  missing_cols <- setdiff(required, names(delineation))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in delineation: ",
         paste(missing_cols, collapse = ", "))
  }
  stopifnot(is.numeric(sr) && length(sr) == 1 && sr > 0)

  channels <- unique(delineation$channel)
  results <- list()

  for (ch in channels) {
    ch_del <- delineation[delineation$channel == ch, ]
    ch_del <- ch_del[order(ch_del$beat), ]
    n_beats <- nrow(ch_del)

    for (i in seq_len(n_beats)) {
      row <- ch_del[i, ]

      # QRS duration
      qrs_ms <- (row$qrs_offset - row$qrs_onset) / sr * 1000

      # PR interval: P peak to QRS onset
      pr_ms <- if (is.na(row$p_peak)) {
        NA_real_
      } else {
        (row$qrs_onset - row$p_peak) / sr * 1000
      }

      # QT interval: QRS onset to T-wave end (clinical standard)
      qt_ms <- if (is.na(row$t_end)) {
        NA_real_
      } else {
        (row$t_end - row$qrs_onset) / sr * 1000
      }

      # RR interval to next beat
      rr_ms <- if (i < n_beats) {
        (ch_del$r_peak[i + 1] - row$r_peak) / sr * 1000
      } else {
        NA_real_
      }

      # QTc (Bazett): QT / sqrt(RR in seconds)
      qtc_ms <- if (!is.na(qt_ms) && !is.na(rr_ms) && rr_ms > 0) {
        rr_sec <- rr_ms / 1000
        qt_ms / sqrt(rr_sec)
      } else {
        NA_real_
      }

      results[[length(results) + 1]] <- data.frame(
        channel = as.integer(ch),
        beat = row$beat,
        pr_ms = pr_ms,
        qt_ms = qt_ms,
        qtc_ms = qtc_ms,
        qrs_ms = qrs_ms,
        rr_ms = rr_ms,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(
      channel = integer(0), beat = integer(0),
      pr_ms = numeric(0), qt_ms = numeric(0), qtc_ms = numeric(0),
      qrs_ms = numeric(0), rr_ms = numeric(0)
    )
  }
}


# --- Internal helpers --------------------------------------------------------

#' Find QRS onset by searching backward from Q-wave trough
#'
#' Strategy: first locate the Q-wave trough (local minimum before R), then
#' estimate baseline from far before Q, and walk backward from Q until
#' the signal amplitude is within a small fraction of the Q deflection from
#' baseline.
#' @keywords internal
.find_qrs_onset <- function(sig, r_pos, sr, n_time) {
  # Step 1: Find Q trough in [-60ms, 0ms] window before R
  q_search <- as.integer(round(0.060 * sr))
  q_start <- max(1L, r_pos - q_search)
  q_segment <- sig[q_start:r_pos]
  q_rel <- which.min(q_segment)
  q_pos <- q_start + q_rel - 1L

  # Step 2: Search backward from Q for onset up to 80ms before Q
  onset_search <- as.integer(round(0.080 * sr))
  search_start <- max(1L, q_pos - onset_search)

  if (search_start >= q_pos) return(q_pos)

  segment <- sig[search_start:q_pos]
  n_seg <- length(segment)

  if (n_seg < 3) return(search_start)

  # Estimate baseline from the first third of search window
  baseline_n <- max(1L, as.integer(n_seg / 3))
  baseline <- stats::median(segment[seq_len(baseline_n)])

  # Q deflection amplitude from baseline
  q_amp <- abs(sig[q_pos] - baseline)
  if (q_amp < 1e-10) return(q_pos)

  # Walk backward from Q: find where signal is within 10% of Q amplitude
  # from baseline (i.e., the QRS deflection has essentially ended)
  amp_thr <- 0.10 * q_amp
  onset_rel <- 1L
  for (k in seq(n_seg, 1L)) {
    if (abs(segment[k] - baseline) < amp_thr) {
      onset_rel <- k
      break
    }
  }

  search_start + onset_rel - 1L
}


#' Find QRS offset by searching forward from S-wave trough
#'
#' Strategy: first locate the S-wave trough (local minimum after R), then
#' estimate baseline from far after S, and walk forward from S until
#' the signal returns to near baseline.
#' @keywords internal
.find_qrs_offset <- function(sig, r_pos, sr, n_time) {
  # Step 1: Find S trough in [0ms, +60ms] window after R
  s_search <- as.integer(round(0.060 * sr))
  s_end <- min(n_time, r_pos + s_search)
  s_segment <- sig[r_pos:s_end]
  s_rel <- which.min(s_segment)
  s_pos <- r_pos + s_rel - 1L

  # Step 2: Search forward from S for offset up to 80ms after S
  offset_search <- as.integer(round(0.080 * sr))
  search_end <- min(n_time, s_pos + offset_search)

  if (search_end <= s_pos) return(s_pos)

  segment <- sig[s_pos:search_end]
  n_seg <- length(segment)

  if (n_seg < 3) return(s_pos)

  # Estimate baseline from the last third of search window
  baseline_n <- max(1L, as.integer(n_seg / 3))
  baseline <- stats::median(segment[seq(n_seg - baseline_n + 1L, n_seg)])

  # S deflection amplitude from baseline
  s_amp <- abs(sig[s_pos] - baseline)
  if (s_amp < 1e-10) return(s_pos)

  # Walk forward from S: find where signal returns to within 10% of S
  # amplitude from baseline
  amp_thr <- 0.10 * s_amp
  offset_rel <- n_seg
  for (k in seq_len(n_seg)) {
    if (abs(segment[k] - baseline) < amp_thr) {
      offset_rel <- k
      break
    }
  }

  s_pos + offset_rel - 1L
}


#' Find T-wave end by searching forward from T-peak for baseline return
#'
#' Uses a tangent-intercept method: draws a tangent at the steepest downslope
#' of the T wave and finds where it crosses the baseline level.
#' Falls back to amplitude threshold if tangent method fails.
#' @param sig Signal vector.
#' @param t_peak Sample index of the T-wave peak.
#' @param sr Sampling rate in Hz.
#' @param n_time Total number of samples.
#' @return Sample index of T-wave end, or NA_integer_.
#' @keywords internal
.find_t_end <- function(sig, t_peak, sr, n_time) {
  # Search window: T-peak to +200ms after T-peak
  search_ms <- 200
  search_samples <- as.integer(round(search_ms / 1000 * sr))
  search_end <- min(n_time, t_peak + search_samples)

  if (search_end <= t_peak + 2L) return(NA_integer_)

  segment <- sig[t_peak:search_end]
  n_seg <- length(segment)

  if (n_seg < 5) return(NA_integer_)

  # Estimate baseline from the last quarter of search window
  baseline_n <- max(2L, as.integer(n_seg / 4))
  baseline <- stats::median(segment[seq(n_seg - baseline_n + 1L, n_seg)])

  # T-peak amplitude above baseline
  t_amp <- segment[1] - baseline
  if (abs(t_amp) < 1e-10) return(t_peak)

  # Find steepest downslope (most negative first derivative)
  deriv <- diff(segment)
  steepest_idx <- which.min(deriv)

  if (length(steepest_idx) == 0 || steepest_idx < 1) {
    # Fallback: amplitude threshold (10% of T amplitude above baseline)
    thr <- baseline + 0.1 * t_amp
    for (k in seq_len(n_seg)) {
      if ((t_amp > 0 && segment[k] <= thr) ||
          (t_amp < 0 && segment[k] >= thr)) {
        return(t_peak + k - 1L)
      }
    }
    return(search_end)
  }

  # Tangent line at steepest point: y = slope * (x - x0) + y0
  slope <- deriv[steepest_idx]
  x0 <- steepest_idx
  y0 <- segment[steepest_idx]

  # Find where tangent crosses baseline: baseline = slope * (x - x0) + y0
  if (abs(slope) < 1e-10) {
    # Nearly flat slope — use amplitude threshold fallback
    thr <- baseline + 0.1 * t_amp
    for (k in seq_len(n_seg)) {
      if ((t_amp > 0 && segment[k] <= thr) ||
          (t_amp < 0 && segment[k] >= thr)) {
        return(t_peak + k - 1L)
      }
    }
    return(search_end)
  }

  x_intercept <- x0 + (baseline - y0) / slope
  t_end_rel <- as.integer(round(x_intercept))

  # Clamp to valid range
  t_end_rel <- max(1L, min(n_seg, t_end_rel))
  t_peak + t_end_rel - 1L
}


#' Find a wave peak (P or T) in a time window relative to R-peak
#' @keywords internal
.find_wave_peak <- function(sig, r_pos, win_start_samples, win_end_samples,
                            n_time, find_max = TRUE) {
  lo <- r_pos + win_start_samples
  hi <- r_pos + win_end_samples

  # Ensure proper ordering
  if (lo > hi) {
    tmp <- lo
    lo <- hi
    hi <- tmp
  }

  lo <- max(1L, lo)
  hi <- min(n_time, hi)

  if (lo > hi) return(NA_integer_)

  segment <- sig[lo:hi]

  if (find_max) {
    peak_rel <- which.max(segment)
  } else {
    peak_rel <- which.min(segment)
  }

  lo + peak_rel - 1L
}
