#' HRV Frequency-Domain Analysis
#'
#' Computes heart rate variability frequency-domain metrics from RR interval
#' data using Welch's method or Lomb-Scargle periodogram.
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @param method Spectral estimation method: \code{"welch"} (default) for
#'   uniformly resampled FFT-based PSD, or \code{"lomb"} for Lomb-Scargle
#'   periodogram on unevenly sampled data.
#' @param vlf_band Numeric vector of length 2 defining VLF band in Hz
#'   (default: c(0.003, 0.04)).
#' @param lf_band Numeric vector of length 2 defining LF band in Hz
#'   (default: c(0.04, 0.15)).
#' @param hf_band Numeric vector of length 2 defining HF band in Hz
#'   (default: c(0.15, 0.4)).
#' @return A data.frame with one row per channel and the following columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{vlf}{Absolute power in the very-low-frequency band (ms^2).}
#'     \item{lf}{Absolute power in the low-frequency band (ms^2), associated
#'       with sympathetic and parasympathetic modulation.}
#'     \item{hf}{Absolute power in the high-frequency band (ms^2), associated
#'       with parasympathetic (vagal) modulation.}
#'     \item{lf_hf_ratio}{Ratio of LF to HF power, or \code{NA} if HF power
#'       is zero.}
#'     \item{total_power}{Sum of VLF, LF, and HF power (ms^2).}
#'   }
#'
#' @references Task Force of the European Society of Cardiology and the North
#'   American Society of Pacing and Electrophysiology (1996). "Heart rate
#'   variability: Standards of measurement, physiological interpretation and
#'   clinical use." \emph{Circulation}, 93(5), 1043--1065.
#'
#' @seealso \code{\link{ecgRRintervals}} for computing RR intervals,
#'   \code{\link{ecgHRVtime}} for time-domain HRV metrics,
#'   \code{\link{ecgHRVnonlinear}} for nonlinear HRV analysis,
#'   \code{\link{ecgRRcorrect}} for ectopic beat correction before analysis.
#'
#' @export
#'
#' @examples
#' n <- 300
#' time_sec <- cumsum(rep(0.85, n))
#' rr_ms <- 850 + 30 * sin(2 * pi * 0.1 * time_sec)
#' rr <- data.frame(channel = rep(1L, n), rr_ms = rr_ms, time_sec = time_sec)
#' result <- ecgHRVfreq(rr, method = "welch")
ecgHRVfreq <- function(rr,
                        method = c("welch", "lomb"),
                        vlf_band = c(0.003, 0.04),
                        lf_band = c(0.04, 0.15),
                        hf_band = c(0.15, 0.4)) {
  stopifnot(is.data.frame(rr))
  required_cols <- c("channel", "rr_ms", "time_sec")
  missing_cols <- setdiff(required_cols, names(rr))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  method <- match.arg(method)

  channels <- unique(rr$channel)
  results <- list()

  for (ch in channels) {
    ch_data <- rr[rr$channel == ch, ]
    t_sec <- ch_data$time_sec
    rr_ms <- ch_data$rr_ms

    if (method == "welch") {
      psd_result <- .hrv_welch(t_sec, rr_ms)
    } else {
      psd_result <- .hrv_lomb(t_sec, rr_ms)
    }

    freqs <- psd_result$freqs
    psd <- psd_result$psd

    vlf <- .band_power(freqs, psd, vlf_band)
    lf <- .band_power(freqs, psd, lf_band)
    hf <- .band_power(freqs, psd, hf_band)
    total_power <- vlf + lf + hf
    lf_hf_ratio <- if (hf > 0) lf / hf else NA_real_

    results[[length(results) + 1]] <- data.frame(
      channel = ch,
      vlf = vlf,
      lf = lf,
      hf = hf,
      lf_hf_ratio = lf_hf_ratio,
      total_power = total_power,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

#' Welch's method for HRV PSD estimation
#' @param t_sec Time stamps in seconds.
#' @param rr_ms RR intervals in milliseconds.
#' @return List with freqs and psd vectors.
#' @keywords internal
.hrv_welch <- function(t_sec, rr_ms) {
  # Resample to uniform 4 Hz grid
  fs <- 4
  t_uniform <- seq(min(t_sec), max(t_sec), by = 1 / fs)
  rr_uniform <- stats::approx(t_sec, rr_ms, xout = t_uniform)$y

  # Detrend: subtract mean
  rr_uniform <- rr_uniform - mean(rr_uniform, na.rm = TRUE)

  n_total <- length(rr_uniform)
  seg_len <- min(256L, n_total)
  overlap_samples <- as.integer(seg_len %/% 2)
  step <- seg_len - overlap_samples

  # Hanning window
  hanning <- 0.5 * (1 - cos(2 * pi * seq(0, seg_len - 1) / (seg_len - 1)))
  win_power <- sum(hanning^2)

  # Frequency axis for one-sided PSD
  n_fft <- seg_len
  n_freq <- n_fft %/% 2 + 1
  freqs <- seq(0, fs / 2, length.out = n_freq)

  # Accumulate PSDs across segments
  psd_sum <- rep(0, n_freq)
  n_segments <- 0L

  start <- 1L
  while (start + seg_len - 1 <= n_total) {
    segment <- rr_uniform[start:(start + seg_len - 1)]
    segment <- segment * hanning
    ft <- fft(segment)
    power <- (Mod(ft)^2) / win_power
    psd_half <- power[seq_len(n_freq)]

    # Scale for one-sided spectrum (double non-DC, non-Nyquist bins)
    if (n_freq > 2) {
      psd_half[2:(n_freq - 1)] <- 2 * psd_half[2:(n_freq - 1)]
    }

    psd_sum <- psd_sum + psd_half
    n_segments <- n_segments + 1L
    start <- start + step
  }

  # If no complete segment, use the whole signal

  if (n_segments == 0L) {
    segment <- rr_uniform * hanning[seq_len(n_total)]
    ft <- fft(segment)
    n_freq <- n_total %/% 2 + 1
    freqs <- seq(0, fs / 2, length.out = n_freq)
    power <- (Mod(ft)^2) / sum(hanning[seq_len(n_total)]^2)
    psd_sum <- power[seq_len(n_freq)]
    if (n_freq > 2) {
      psd_sum[2:(n_freq - 1)] <- 2 * psd_sum[2:(n_freq - 1)]
    }
    n_segments <- 1L
  }

  psd_avg <- psd_sum / (n_segments * fs)

  # Remove DC
  psd_avg[1] <- 0

  list(freqs = freqs, psd = psd_avg)
}

#' Lomb-Scargle periodogram for HRV PSD estimation
#' @param t_sec Time stamps in seconds.
#' @param rr_ms RR intervals in milliseconds.
#' @return List with freqs and psd vectors.
#' @keywords internal
.hrv_lomb <- function(t_sec, rr_ms) {
  # Detrend: subtract mean
  x <- rr_ms - mean(rr_ms, na.rm = TRUE)
  t <- t_sec

  # Test frequencies: 512 evenly spaced from 0.003 to 0.4 Hz
  n_freq <- 512L
  freqs <- seq(0.003, 0.4, length.out = n_freq)
  psd <- numeric(n_freq)

  two_pi <- 2 * pi

  for (i in seq_len(n_freq)) {
    f <- freqs[i]
    omega <- two_pi * f

    # Compute tau
    sin2 <- sum(sin(2 * omega * t))
    cos2 <- sum(cos(2 * omega * t))
    tau <- atan2(sin2, cos2) / (2 * omega)

    # Compute Lomb-Scargle power
    cos_term <- cos(omega * (t - tau))
    sin_term <- sin(omega * (t - tau))

    num_cos <- sum(x * cos_term)^2
    den_cos <- sum(cos_term^2)

    num_sin <- sum(x * sin_term)^2
    den_sin <- sum(sin_term^2)

    psd[i] <- (num_cos / den_cos + num_sin / den_sin) / 2
  }

  list(freqs = freqs, psd = psd)
}

#' Integrate PSD within a frequency band
#' @param freqs Frequency vector.
#' @param psd PSD vector.
#' @param band Numeric vector of length 2 (low, high) in Hz.
#' @return Numeric scalar: integrated power in the band.
#' @keywords internal
.band_power <- function(freqs, psd, band) {
  idx <- which(freqs >= band[1] & freqs <= band[2])
  if (length(idx) < 2) return(0)
  df <- freqs[idx[2]] - freqs[idx[1]]
  sum(psd[idx]) * df
}
