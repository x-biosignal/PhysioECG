#' Detect R-Peaks in ECG Signal
#'
#' Identifies R-peaks in ECG data using one of four QRS detectors. All share a
#' bandpass front end and automatic inverted-lead handling but differ in their
#' feature transform and threshold rule: \code{"pan_tompkins"} (5-15 Hz,
#' derivative-square-integrate, adaptive dual threshold), \code{"hamilton"}
#' (8-16 Hz, rectified derivative, 80 ms integration, running mean of the last
#' eight QRS/noise peaks), \code{"elgendi"} (8-20 Hz, squared, two moving
#' averages with an event-block threshold), and \code{"christov"} (9-30 Hz,
#' rectified derivative smoothed over 30 ms, adaptive steep-slope M-threshold
#' with decay).
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param method Detection method: one of \code{"pan_tompkins"} (default),
#'   \code{"hamilton"}, \code{"elgendi"} or \code{"christov"}.
#' @param threshold_factor Retained for backward compatibility; the
#'   Pan-Tompkins path uses its built-in fractional threshold (default: 0.6).
#' @param refractory_ms Refractory period in milliseconds. No two peaks can
#'   be closer than this (default: 200).
#' @param assay_name Name of the assay to use. If \code{NULL}, the default
#'   assay is used.
#' @param beta Offset factor for the Elgendi event threshold
#'   \code{MA_beat + beta * mean(squared signal)} (default: 0.08); ignored by
#'   the other methods.
#' @return A data.frame with one row per detected R-peak and the following
#'   columns:
#'   \describe{
#'     \item{channel}{Integer channel index (1-based).}
#'     \item{sample}{Integer sample index of the R-peak within the assay
#'       matrix.}
#'     \item{time_sec}{Time of the R-peak in seconds from signal onset.}
#'     \item{amplitude}{Amplitude of the raw signal at the R-peak location
#'       (in original units, not inverted).}
#'   }
#'   Returns a zero-row data.frame with the same column structure if no peaks
#'   are detected.
#'
#' @references Pan, J. & Tompkins, W.J. (1985). "A real-time QRS detection
#'   algorithm." \emph{IEEE Transactions on Biomedical Engineering}, 32(3),
#'   230--236. \doi{10.1109/TBME.1985.325532}
#'
#'   Hamilton, P. (2002). "Open source ECG analysis." \emph{Computers in
#'   Cardiology}, 29, 101--104.
#'
#'   Elgendi, M. (2013). "Fast QRS detection with an optimized knowledge-based
#'   method." \emph{PLoS ONE}, 8(9), e73557. \doi{10.1371/journal.pone.0073557}
#'
#'   Christov, I. (2004). "Real time electrocardiogram QRS detection using
#'   combined adaptive threshold." \emph{BioMedical Engineering OnLine}, 3, 28.
#'   \doi{10.1186/1475-925X-3-28}
#'
#' @seealso \code{\link{ecgRRintervals}} for computing RR intervals from
#'   detected peaks, \code{\link{ecgBeatSQI}} for detector-agreement quality,
#'   \code{\link{ecgDelineate}} for full waveform morphology
#'   analysis, \code{\link{ecgSignalQuality}} for signal quality assessment.
#'
#' @export
#' @examples
#' \dontrun{
#' pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
#' peaks <- ecgDetectRpeaks(pe)
#' head(peaks)
#' }
ecgDetectRpeaks <- function(x, method = c("pan_tompkins", "hamilton",
                                          "elgendi", "christov"),
                            threshold_factor = 0.6,
                            refractory_ms = 200,
                            assay_name = NULL,
                            beta = 0.08) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)

  results <- list()
  for (ch in seq_len(ncol(data))) {
    sig <- data[, ch]
    refined <- .ecg_run_detector(sig, sr, method,
                                 threshold_factor = threshold_factor,
                                 refractory_ms = refractory_ms, beta = beta)
    if (length(refined) > 0) {
      orig_sig <- data[, ch]
      results[[length(results) + 1]] <- data.frame(
        channel = rep(as.integer(ch), length(refined)),
        sample = refined,
        time_sec = (refined - 1) / sr,
        amplitude = orig_sig[refined],
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(channel = integer(0), sample = integer(0),
               time_sec = numeric(0), amplitude = numeric(0))
  }
}

# Dispatch a single-channel signal to the requested detector; returns refined
# R-peak sample indices (integer vector).
.ecg_run_detector <- function(sig, sr, method, threshold_factor = 0.6,
                              refractory_ms = 200, beta = 0.08) {
  switch(method,
    pan_tompkins = .ecg_detect_pan_tompkins(sig, sr, threshold_factor, refractory_ms),
    hamilton     = .ecg_detect_hamilton(sig, sr, refractory_ms),
    elgendi      = .ecg_detect_elgendi(sig, sr, beta, refractory_ms),
    christov     = .ecg_detect_christov(sig, sr, refractory_ms),
    stop("unknown detector method: ", method)
  )
}

# Snap candidate feature peaks to the true R-peak. Detector feature streams lag
# the QRS, so the search window is biased backward; this aligns every detector's
# output on the same R-peak (essential for cross-detector agreement / bSQI).
# Shifting all peaks by a common delay leaves RR intervals unchanged.
#
# `sig` must be the ORIGINAL (un-flipped) lead: polarity is decided here, per
# beat-set, biased toward the positive R. A normal lead with a deep S wave (e.g.
# precordial V1/V2, where the global |min| exceeds |max|) is NOT treated as
# inverted, so the fiducial and reported amplitude land on the R-peak rather
# than the S trough. A genuinely inverted lead (positive excursion < half the
# negative one) snaps to the dominant negative deflection.
.ecg_refine_peaks <- function(sig, peak_indices, sr) {
  if (length(peak_indices) == 0) return(integer(0))
  n_time <- length(sig)
  back <- max(1L, as.integer(round(0.22 * sr)))
  fwd  <- max(1L, as.integer(round(0.06 * sr)))
  lo <- pmax(1L, peak_indices - back)
  hi <- pmin(n_time, peak_indices + fwd)
  base <- stats::median(sig)
  posv <- negv <- numeric(length(peak_indices))
  for (j in seq_along(peak_indices)) {
    seg <- sig[lo[j]:hi[j]] - base
    posv[j] <- max(seg)
    negv[j] <- abs(min(seg))
  }
  inverted <- stats::median(posv) < 0.5 * stats::median(negv)
  refined <- integer(length(peak_indices))
  for (j in seq_along(peak_indices)) {
    seg <- sig[lo[j]:hi[j]]
    refined[j] <- if (inverted) lo[j] + which.min(seg) - 1L
                  else          lo[j] + which.max(seg) - 1L
  }
  sort(unique(refined))
}

# Keep peaks in ascending order with a minimum separation (drop the later of
# any pair closer than refractory_samples).
.ecg_apply_refractory <- function(peaks, refractory_samples) {
  peaks <- sort(unique(peaks))
  if (length(peaks) <= 1) return(peaks)
  keep <- integer(0)
  last <- -Inf
  for (p in peaks) {
    if (p - last >= refractory_samples) { keep <- c(keep, p); last <- p }
  }
  keep
}

# --- Pan-Tompkins (1985): derivative-square-integrate + dual EMA threshold ---
.ecg_detect_pan_tompkins <- function(sig, sr, threshold_factor = 0.6,
                                     refractory_ms = 200) {
  refractory_samples <- max(1L, as.integer(round(refractory_ms / 1000 * sr)))
  int_window <- max(1L, as.integer(round(0.150 * sr)))
  orig_sig <- sig
  if (abs(min(sig)) > abs(max(sig))) sig <- -sig

  bp_sig <- .fir_bandpass(sig, sr, low = 5, high = 15, order = 65)
  dsig <- c(diff(bp_sig), 0)
  sq_sig <- dsig^2
  int_sig <- as.numeric(stats::filter(sq_sig, rep(1 / int_window, int_window),
                                      sides = 1))
  int_sig[is.na(int_sig)] <- 0

  local_max_idx <- .find_local_maxima(int_sig, refractory_samples)
  if (length(local_max_idx) == 0) return(integer(0))

  init_max <- max(int_sig[seq_len(min(length(int_sig),
                                      2L * as.integer(round(sr))))])
  signal_peak <- 0.5 * init_max
  noise_peak <- 0.1 * init_max
  alpha <- 0.125
  peak_indices <- integer(0)
  last_peak_sample <- -refractory_samples
  for (idx in local_max_idx) {
    if ((idx - last_peak_sample) < refractory_samples) next
    threshold1 <- noise_peak + 0.25 * (signal_peak - noise_peak)
    val <- int_sig[idx]
    if (val > threshold1) {
      peak_indices <- c(peak_indices, idx)
      last_peak_sample <- idx
      signal_peak <- alpha * val + (1 - alpha) * signal_peak
    } else {
      noise_peak <- alpha * val + (1 - alpha) * noise_peak
    }
  }
  .ecg_refine_peaks(orig_sig, peak_indices, sr)
}

# --- Hamilton (2002): rectified derivative, 80 ms integration, running mean-8 ---
.ecg_detect_hamilton <- function(sig, sr, refractory_ms = 200) {
  refractory_samples <- max(1L, as.integer(round(refractory_ms / 1000 * sr)))
  int_window <- max(1L, as.integer(round(0.080 * sr)))
  orig_sig <- sig
  if (abs(min(sig)) > abs(max(sig))) sig <- -sig

  bp_sig <- .fir_bandpass(sig, sr, low = 8, high = 16, order = 65)
  rect <- abs(c(diff(bp_sig), 0))                 # rectify (not square)
  int_sig <- as.numeric(stats::filter(rect, rep(1 / int_window, int_window),
                                      sides = 1))
  int_sig[is.na(int_sig)] <- 0

  local_max_idx <- .find_local_maxima(int_sig, refractory_samples)
  if (length(local_max_idx) == 0) return(integer(0))

  init_max <- max(int_sig[seq_len(min(length(int_sig),
                                      2L * as.integer(round(sr))))])
  qrs_buf <- rep(0.5 * init_max, 8L)
  noise_buf <- rep(0.1 * init_max, 8L)
  qrs_mean <- mean(qrs_buf); noise_mean <- mean(noise_buf)
  tf <- 0.3125                                     # Hamilton detection fraction
  peaks <- integer(0); last <- -refractory_samples
  for (idx in local_max_idx) {
    val <- int_sig[idx]
    dt <- noise_mean + tf * (qrs_mean - noise_mean)
    if (val > dt && (idx - last) >= refractory_samples) {
      peaks <- c(peaks, idx); last <- idx
      qrs_buf <- c(qrs_buf[-1L], val); qrs_mean <- mean(qrs_buf)
    } else {
      noise_buf <- c(noise_buf[-1L], val); noise_mean <- mean(noise_buf)
    }
  }
  .ecg_refine_peaks(orig_sig, peaks, sr)
}

# --- Elgendi (2013): two moving averages with an event-block threshold ---
.ecg_detect_elgendi <- function(sig, sr, beta = 0.08, refractory_ms = 200) {
  refractory_samples <- max(1L, as.integer(round(refractory_ms / 1000 * sr)))
  orig_sig <- sig
  if (abs(min(sig)) > abs(max(sig))) sig <- -sig

  bp_sig <- .fir_bandpass(sig, sr, low = 8, high = 20, order = 65)
  y <- bp_sig^2
  w1 <- max(1L, as.integer(round(0.097 * sr)))    # QRS window
  w2 <- max(1L, as.integer(round(0.611 * sr)))    # beat window
  ma_qrs <- as.numeric(stats::filter(y, rep(1 / w1, w1), sides = 2))
  ma_beat <- as.numeric(stats::filter(y, rep(1 / w2, w2), sides = 2))
  ma_qrs[is.na(ma_qrs)] <- 0
  ma_beat[is.na(ma_beat)] <- 0

  thr1 <- ma_beat + beta * mean(y)
  blocks <- ma_qrs > thr1
  r <- rle(blocks)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  peaks <- integer(0)
  for (k in seq_along(r$values)) {
    if (isTRUE(r$values[k]) && r$lengths[k] >= w1) {  # accept wide blocks only
      lo <- starts[k]; hi <- ends[k]
      peaks <- c(peaks, lo + which.max(ma_qrs[lo:hi]) - 1L)   # feature peak
    }
  }
  peaks <- .ecg_apply_refractory(peaks, refractory_samples)
  .ecg_refine_peaks(orig_sig, peaks, sr)   # align to raw R-peak like the others
}

# --- Christov (2004): steep-slope adaptive M-threshold with decay ---
# Implements the dominant M (steep-slope) term of Christov's combined adaptive
# threshold with its 200 ms -> 1200 ms decay; the slower baseline (F) and rate
# (R) terms are omitted, which is sufficient for clean-signal detection.
.ecg_detect_christov <- function(sig, sr, refractory_ms = 200) {
  refractory_samples <- max(1L, as.integer(round(refractory_ms / 1000 * sr)))
  orig_sig <- sig
  if (abs(min(sig)) > abs(max(sig))) sig <- -sig

  bp_sig <- .fir_bandpass(sig, sr, low = 9, high = 30, order = 65)
  dY <- abs(c(diff(bp_sig), 0))
  g <- max(1L, as.integer(round(0.030 * sr)))     # 30 ms smoother
  y <- as.numeric(stats::filter(dY, rep(1 / g, g), sides = 1))
  y[is.na(y)] <- 0

  init_win <- min(length(y), as.integer(round(5 * sr)))
  m <- 0.6 * max(y[seq_len(init_win)])
  cand <- .find_local_maxima(y, refractory_samples)
  if (length(cand) == 0) return(integer(0))

  peaks <- integer(0); last <- -refractory_samples
  for (idx in cand) {
    dt <- (idx - last) / sr
    if (dt <= 0.200) next                          # absolute refractory
    m_th <- if (dt <= 1.200) m * (1 - 0.4 * (dt - 0.200) / 1.000) else 0.6 * m
    if (y[idx] > m_th) {
      peaks <- c(peaks, idx); last <- idx
      m <- 0.6 * m + 0.4 * (0.6 * y[idx])          # EMA update of M
    }
  }
  .ecg_refine_peaks(orig_sig, peaks, sr)
}


#' Find local maxima in a signal with minimum separation
#' @keywords internal
.find_local_maxima <- function(sig, min_distance) {
  n <- length(sig)
  if (n < 3) return(integer(0))

  # Find all local maxima (strictly greater than both neighbours)
  maxima <- integer(0)
  for (i in 2:(n - 1)) {
    if (sig[i] > sig[i - 1] && sig[i] >= sig[i + 1] && sig[i] > 0) {
      maxima <- c(maxima, i)
    }
  }

  if (length(maxima) == 0) return(integer(0))

  # Sort by amplitude (descending) and greedily keep peaks with min_distance
  ord <- order(sig[maxima], decreasing = TRUE)
  maxima <- maxima[ord]
  kept <- logical(length(maxima))
  kept[1] <- TRUE

  # seq_len(n)[-1] is empty for a single maximum (avoids the 2:1 reverse loop
  # that would dereference the non-existent maxima[2]).
  for (i in seq_len(length(maxima))[-1]) {
    if (all(abs(maxima[i] - maxima[which(kept)]) >= min_distance)) {
      kept[i] <- TRUE
    }
  }

  sort(maxima[kept])
}


#' Compute RR Intervals from Detected R-Peaks
#'
#' Calculates the time intervals between consecutive R-peaks for each channel.
#' The resulting RR interval series is the standard input for all HRV analysis
#' functions in this package.
#'
#' @param x A PhysioExperiment object.
#' @param peaks A data.frame of detected peaks as returned by
#'   \code{\link{ecgDetectRpeaks}}, with columns \code{channel},
#'   \code{sample}, and \code{time_sec}.
#' @return A data.frame with one row per consecutive beat pair and the
#'   following columns:
#'   \describe{
#'     \item{channel}{Integer channel index (1-based).}
#'     \item{rr_ms}{RR interval in milliseconds (time between successive
#'       R-peaks).}
#'     \item{time_sec}{Time of the first beat in each pair (seconds from
#'       signal onset).}
#'   }
#'   Returns a zero-row data.frame with the same column structure if fewer
#'   than two peaks are available.
#'
#' @references Pan, J. & Tompkins, W.J. (1985). "A real-time QRS detection
#'   algorithm." \emph{IEEE Transactions on Biomedical Engineering}, 32(3),
#'   230--236. \doi{10.1109/TBME.1985.325532}
#'
#' @seealso \code{\link{ecgDetectRpeaks}} for R-peak detection,
#'   \code{\link{ecgHRVtime}} for time-domain HRV metrics,
#'   \code{\link{ecgHRVfreq}} for frequency-domain HRV analysis,
#'   \code{\link{ecgQualityCheck}} for ectopic beat detection.
#'
#' @export
#' @examples
#' \dontrun{
#' pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 60)
#' peaks <- ecgDetectRpeaks(pe)
#' rr <- ecgRRintervals(pe, peaks)
#' head(rr)
#' }
ecgRRintervals <- function(x, peaks) {
  stopifnot(inherits(x, "PhysioExperiment"))
  stopifnot(is.data.frame(peaks))
  stopifnot(all(c("channel", "sample", "time_sec") %in% names(peaks)))

  sr <- samplingRate(x)
  channels <- unique(peaks$channel)
  results <- list()

  for (ch in channels) {
    ch_peaks <- peaks[peaks$channel == ch, ]
    ch_peaks <- ch_peaks[order(ch_peaks$sample), ]

    if (nrow(ch_peaks) < 2) next

    diffs <- diff(ch_peaks$sample)
    rr_ms <- diffs / sr * 1000

    results[[length(results) + 1]] <- data.frame(
      channel = rep(as.integer(ch), length(rr_ms)),
      rr_ms = rr_ms,
      time_sec = ch_peaks$time_sec[-nrow(ch_peaks)],
      stringsAsFactors = FALSE
    )
  }

  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame(channel = integer(0), rr_ms = numeric(0),
               time_sec = numeric(0))
  }
}


# --- Internal helper: windowed-sinc FIR bandpass filter ---

#' @keywords internal
.fir_bandpass <- function(sig, sr, low, high, order = 65) {
  n <- length(sig)
  # Ensure odd order for symmetric filter
  if (order %% 2 == 0) order <- order + 1L
  half <- (order - 1L) %/% 2L

  # Normalised cutoffs in cycles/sample (fraction of the sampling rate) -- the
  # convention required by the sin(2*pi*fc*m)/(pi*m) windowed-sinc kernel below.
  # (Using fraction-of-Nyquist here would double every cutoff.)
  fc_low <- low / sr
  fc_high <- high / sr

  # Sinc-based bandpass filter coefficients
  m <- seq(-half, half)
  # Avoid division by zero at centre tap
  hp <- ifelse(m == 0, 2 * fc_high, sin(2 * pi * fc_high * m) / (pi * m))
  lp <- ifelse(m == 0, 2 * fc_low, sin(2 * pi * fc_low * m) / (pi * m))
  bp <- hp - lp

  # Apply Hamming window
  w <- 0.54 - 0.46 * cos(2 * pi * seq(0, order - 1) / (order - 1))
  bp <- bp * w

  # Normalise to unit gain at centre frequency
  centre_freq <- (low + high) / 2
  centre_response <- sum(bp * cos(2 * pi * centre_freq / sr * m))
  if (abs(centre_response) > 1e-10) bp <- bp / centre_response

  # Apply filter (causal via stats::filter, then compensate for group delay)
  filtered <- as.numeric(stats::filter(sig, bp, sides = 1))
  # Replace leading NAs from filter delay with zeros
  filtered[is.na(filtered)] <- 0
  filtered
}

# Greedy one-to-one nearest matching of two sorted index vectors within a
# window; returns the number of matched beats.
.ecg_match_beats <- function(a, b, win_samples) {
  if (length(a) == 0 || length(b) == 0) return(0L)
  a <- sort(a); b <- sort(b)
  used_b <- logical(length(b))
  matched <- 0L
  for (ai in a) {
    d <- abs(b - ai)
    d[used_b] <- Inf
    j <- which.min(d)
    if (is.finite(d[j]) && d[j] <= win_samples) {
      used_b[j] <- TRUE
      matched <- matched + 1L
    }
  }
  matched
}

#' Beat-level Signal Quality Index (bSQI)
#'
#' Quantifies ECG quality per channel by the agreement between two independent
#' QRS detectors. A clean signal makes both detectors agree on nearly every
#' beat (bSQI near 1); noise makes them disagree (bSQI drops).
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param detector_a,detector_b Detector methods to compare (see
#'   \code{\link{ecgDetectRpeaks}}); defaults \code{"pan_tompkins"} and
#'   \code{"elgendi"}.
#' @param match_window_ms Two beats match when within this window (default 150).
#' @param threshold_factor,refractory_ms,beta Passed to the detectors.
#' @param assay_name Assay to use (default: the default assay).
#' @return A data.frame with one row per channel: \code{n_a}, \code{n_b}
#'   (beat counts), \code{n_matched}, and \code{bSQI = matched / (n_a + n_b -
#'   matched)} (Jaccard agreement in \code{[0, 1]}; \code{NA} for a silent
#'   channel).
#' @references Li, Q., Mark, R.G. & Clifford, G.D. (2008). Robust heart rate
#'   estimation from multiple asynchronous noisy sources. \emph{Physiological
#'   Measurement}, 29(1), 15-32.
#' @seealso \code{\link{ecgDetectRpeaks}}, \code{\link{ecgBeatSQIgate}}
#' @export
#' @examples
#' pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
#' ecgBeatSQI(pe)
ecgBeatSQI <- function(x, detector_a = "pan_tompkins", detector_b = "elgendi",
                       match_window_ms = 150, threshold_factor = 0.6,
                       refractory_ms = 200, beta = 0.08, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  win_samples <- max(1L, as.integer(round(match_window_ms / 1000 * sr)))

  results <- list()
  for (ch in seq_len(ncol(data))) {
    sig <- data[, ch]
    a <- .ecg_run_detector(sig, sr, detector_a, threshold_factor, refractory_ms, beta)
    b <- .ecg_run_detector(sig, sr, detector_b, threshold_factor, refractory_ms, beta)
    n_a <- length(a); n_b <- length(b)
    nm <- .ecg_match_beats(a, b, win_samples)
    bsqi <- if (n_a + n_b == 0) NA_real_ else nm / (n_a + n_b - nm)
    results[[length(results) + 1]] <- data.frame(
      channel = as.integer(ch), n_a = n_a, n_b = n_b,
      n_matched = nm, bSQI = bsqi, stringsAsFactors = FALSE
    )
  }
  do.call(rbind, results)
}

#' Accept/reject ECG channels by bSQI
#'
#' Thin wrapper over \code{\link{ecgBeatSQI}} adding a logical \code{accept}
#' column for channels whose bSQI meets a threshold.
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param threshold Minimum bSQI to accept a channel (default 0.8).
#' @param ... Passed to \code{\link{ecgBeatSQI}}.
#' @return The \code{\link{ecgBeatSQI}} data.frame with an added logical
#'   \code{accept} column.
#' @seealso \code{\link{ecgBeatSQI}}
#' @export
#' @examples
#' pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
#' ecgBeatSQIgate(pe, threshold = 0.8)
ecgBeatSQIgate <- function(x, threshold = 0.8, ...) {
  sqi <- ecgBeatSQI(x, ...)
  sqi$accept <- !is.na(sqi$bSQI) & sqi$bSQI >= threshold
  sqi
}
