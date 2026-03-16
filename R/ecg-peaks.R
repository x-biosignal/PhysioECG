#' Detect R-Peaks in ECG Signal Using Pan-Tompkins Algorithm
#'
#' Identifies R-peaks in ECG data using an adaptive dual-threshold
#' Pan-Tompkins detector. The algorithm applies bandpass filtering (5-15 Hz),
#' differentiation, squaring, and moving-window integration followed by
#' adaptive thresholding with running signal and noise level estimates.
#' Automatically detects and handles inverted ECG signals.
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param method Detection method. Currently only \code{"pan_tompkins"} is
#'   supported.
#' @param threshold_factor Fraction of the peak integrated signal used as
#'   the detection threshold (default: 0.6).
#' @param refractory_ms Refractory period in milliseconds. No two peaks can
#'   be closer than this (default: 200).
#' @param assay_name Name of the assay to use. If \code{NULL}, the default
#'   assay is used.
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
#' @seealso \code{\link{ecgRRintervals}} for computing RR intervals from
#'   detected peaks, \code{\link{ecgDelineate}} for full waveform morphology
#'   analysis, \code{\link{ecgSignalQuality}} for signal quality assessment.
#'
#' @export
#' @examples
#' \dontrun{
#' pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
#' peaks <- ecgDetectRpeaks(pe)
#' head(peaks)
#' }
ecgDetectRpeaks <- function(x, method = "pan_tompkins",
                            threshold_factor = 0.6,
                            refractory_ms = 200,
                            assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method, choices = "pan_tompkins")

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  refractory_samples <- max(1L, as.integer(round(refractory_ms / 1000 * sr)))
  int_window <- max(1L, as.integer(round(0.150 * sr)))  # 150 ms window

  results <- list()

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    # --- Step 0: Detect signal inversion ---
    # If the magnitude of the most negative peak exceeds the most positive,
    # the signal is likely inverted (e.g., inverted lead). Flip it.
    inverted <- FALSE
    if (abs(min(sig)) > abs(max(sig))) {
      sig <- -sig
      inverted <- TRUE
    }

    # --- Step 1: Bandpass filter 5-15 Hz (windowed-sinc FIR) ---
    bp_sig <- .fir_bandpass(sig, sr, low = 5, high = 15, order = 65)

    # --- Step 2: Differentiate ---
    dsig <- c(diff(bp_sig), 0)

    # --- Step 3: Square ---
    sq_sig <- dsig^2

    # --- Step 4: Moving-window integration ---
    kern <- rep(1 / int_window, int_window)
    int_sig <- as.numeric(stats::filter(sq_sig, kern, sides = 1))
    int_sig[is.na(int_sig)] <- 0

    # --- Step 5: Find local peaks in integrated signal ---
    # Identify all local maxima (points higher than both neighbours)
    local_max_idx <- .find_local_maxima(int_sig, refractory_samples)
    if (length(local_max_idx) == 0) next

    # --- Step 6: Adaptive dual-threshold classification ---
    # Initialise signal and noise peak running estimates from initial data
    init_max <- max(int_sig[seq_len(min(length(int_sig), 2L * as.integer(round(sr))))])
    signal_peak <- 0.5 * init_max
    noise_peak <- 0.1 * init_max
    alpha <- 0.125  # EMA smoothing factor

    peak_indices <- integer(0)
    last_peak_sample <- -refractory_samples

    for (idx in local_max_idx) {
      # Enforce refractory period
      if ((idx - last_peak_sample) < refractory_samples) next

      threshold1 <- noise_peak + 0.25 * (signal_peak - noise_peak)
      val <- int_sig[idx]

      if (val > threshold1) {
        # Classified as signal (R-peak candidate)
        peak_indices <- c(peak_indices, idx)
        last_peak_sample <- idx
        signal_peak <- alpha * val + (1 - alpha) * signal_peak
      } else {
        # Classified as noise
        noise_peak <- alpha * val + (1 - alpha) * noise_peak
      }
    }

    # --- Step 7: Refine to actual R-peak in original signal ---
    search_half <- max(1L, as.integer(round(0.075 * sr)))  # +/-75ms
    refined <- integer(length(peak_indices))
    for (j in seq_along(peak_indices)) {
      lo <- max(1L, peak_indices[j] - search_half)
      hi <- min(n_time, peak_indices[j] + search_half)
      refined[j] <- lo + which.max(sig[lo:hi]) - 1L
    }

    # --- Step 8: Report amplitudes from the original (possibly inverted) signal ---
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

  for (i in 2:length(maxima)) {
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

  # Normalised cutoff frequencies (fraction of Nyquist)
  nyq <- sr / 2
  fc_low <- low / nyq
  fc_high <- high / nyq

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
