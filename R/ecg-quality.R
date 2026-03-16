#' Detect Ectopic Beats in RR Interval Data
#'
#' Identifies ectopic (abnormal) beats by comparing each RR interval to the
#' local median computed over a sliding window of 5 beats. Beats with deviation
#' exceeding the threshold are marked as ectopic.
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @param threshold_ms Maximum allowed deviation from local median in
#'   milliseconds (default: 300).
#' @return The input data.frame with an additional logical column
#'   \code{is_ectopic}.
#'
#' @references Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006).
#'   \emph{Advanced Methods and Tools for ECG Data Analysis}. Artech House.
#'
#'   Task Force of the European Society of Cardiology and the North American
#'   Society of Pacing and Electrophysiology (1996). "Heart rate variability:
#'   Standards of measurement, physiological interpretation and clinical use."
#'   \emph{Circulation}, 93(5), 1043--1065.
#'
#' @seealso \code{\link{ecgRRcorrect}} for correcting detected ectopic beats,
#'   \code{\link{ecgRRintervals}} for computing RR intervals,
#'   \code{\link{ecgHRVtime}} for time-domain HRV analysis,
#'   \code{\link{ecgSignalQuality}} for signal quality assessment.
#'
#' @export
ecgQualityCheck <- function(rr, threshold_ms = 300) {
  stopifnot(is.data.frame(rr))
  required <- c("channel", "rr_ms", "time_sec")
  missing_cols <- setdiff(required, names(rr))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  rr$is_ectopic <- FALSE

  channels <- unique(rr$channel)

  for (ch in channels) {
    idx <- which(rr$channel == ch)
    rr_vals <- rr$rr_ms[idx]
    n <- length(rr_vals)

    for (i in seq_len(n)) {
      win_start <- max(1L, i - 2L)
      win_end <- min(n, i + 2L)
      # Exclude the current beat from its own median calculation
      win_idx <- setdiff(win_start:win_end, i)
      if (length(win_idx) == 0) next
      local_med <- stats::median(rr_vals[win_idx])
      if (abs(rr_vals[i] - local_med) > threshold_ms) {
        rr$is_ectopic[idx[i]] <- TRUE
      }
    }
  }

  rr
}

#' Correct Ectopic Beats in RR Interval Data
#'
#' Replaces or removes ectopic beats identified by
#' \code{\link{ecgQualityCheck}}. The \code{"interpolate"} method uses linear
#' interpolation from surrounding non-ectopic intervals, while
#' \code{"remove"} simply drops ectopic rows.
#'
#' @param rr A data.frame with an \code{is_ectopic} logical column, as
#'   returned by \code{\link{ecgQualityCheck}}.
#' @param method Correction method: \code{"interpolate"} (default) replaces
#'   ectopic values with linear interpolation; \code{"remove"} drops ectopic
#'   rows.
#' @return A data.frame with corrected RR intervals. The \code{is_ectopic}
#'   column is removed from the output.
#'
#' @references Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006).
#'   \emph{Advanced Methods and Tools for ECG Data Analysis}. Artech House.
#'
#' @seealso \code{\link{ecgQualityCheck}} for detecting ectopic beats,
#'   \code{\link{ecgRRintervals}} for computing RR intervals,
#'   \code{\link{ecgHRVtime}} for time-domain HRV analysis.
#'
#' @export
ecgRRcorrect <- function(rr, method = c("interpolate", "remove")) {
  stopifnot(is.data.frame(rr))
  if (!"is_ectopic" %in% names(rr)) {
    stop("Missing required column: is_ectopic. Run ecgQualityCheck first.")
  }
  method <- match.arg(method)

  if (method == "remove") {
    out <- rr[!rr$is_ectopic, , drop = FALSE]
    out$is_ectopic <- NULL
    rownames(out) <- NULL
    return(out)
  }

  # Interpolation method
  channels <- unique(rr$channel)

  for (ch in channels) {
    idx <- which(rr$channel == ch)
    rr_vals <- rr$rr_ms[idx]
    ectopic <- rr$is_ectopic[idx]
    n <- length(rr_vals)

    valid <- which(!ectopic)
    invalid <- which(ectopic)

    if (length(invalid) == 0 || length(valid) == 0) next

    if (length(valid) >= 2) {
      interpolated <- stats::approx(
        x = valid, y = rr_vals[valid],
        xout = invalid, rule = 2
      )$y
    } else {
      interpolated <- rep(rr_vals[valid[1]], length(invalid))
    }

    rr$rr_ms[idx[invalid]] <- interpolated
  }

  rr$is_ectopic <- NULL
  rr
}


#' Assess ECG Signal Quality Per Channel
#'
#' Computes per-channel signal quality metrics for ECG data stored in a
#' PhysioExperiment object.  When detected R-peak locations are provided the
#' signal-to-noise ratio is estimated from QRS vs.\ baseline power; otherwise a
#' variance-based estimate is used.
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param peaks Optional data.frame of detected R-peaks as returned by
#'   \code{\link{ecgDetectRpeaks}}, with at least columns \code{channel} and
#'   \code{sample}.
#' @param assay_name Name of the assay to use.
#'   If \code{NULL} the default assay is used.
#' @return A data.frame with one row per channel and columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{snr_db}{Signal-to-noise ratio in decibels.}
#'     \item{baseline_wander}{RMS amplitude of the low-frequency drift.}
#'     \item{saturation_ratio}{Fraction of samples within 1\% of signal
#'       min or max.}
#'     \item{quality_score}{Composite quality score in the range \[0, 1\],
#'       where 1 indicates excellent quality.}
#'   }
#'
#' @references Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006).
#'   \emph{Advanced Methods and Tools for ECG Data Analysis}. Artech House.
#'
#'   Shaffer, F. & Ginsberg, J.P. (2017). "An overview of heart rate
#'   variability metrics and norms." \emph{Frontiers in Public Health}, 5, 258.
#'   \doi{10.3389/fpubh.2017.00258}
#'
#' @seealso \code{\link{ecgDetectRpeaks}} for R-peak detection,
#'   \code{\link{ecgBaselineCorrect}} for baseline wander correction,
#'   \code{\link{ecgQualityCheck}} for ectopic beat detection.
#'
#' @export
ecgSignalQuality <- function(x, peaks = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (!is.null(peaks)) {
    stopifnot(is.data.frame(peaks))
    stopifnot(all(c("channel", "sample") %in% names(peaks)))
  }

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  results <- vector("list", n_channels)

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    # --- SNR estimation ---
    if (!is.null(peaks)) {
      ch_peaks <- peaks$sample[peaks$channel == ch]
      if (length(ch_peaks) > 0) {
        # QRS half-window: ~50 ms each side
        qrs_half <- max(1L, as.integer(round(0.05 * sr)))
        qrs_idx <- unique(unlist(lapply(ch_peaks, function(pk) {
          seq(max(1L, pk - qrs_half), min(n_time, pk + qrs_half))
        })))
        baseline_idx <- setdiff(seq_len(n_time), qrs_idx)
        if (length(baseline_idx) > 0 && length(qrs_idx) > 0) {
          p_signal <- mean(sig[qrs_idx]^2)
          p_noise <- mean(sig[baseline_idx]^2)
          snr_db <- if (p_noise > 0) 10 * log10(p_signal / p_noise) else Inf
        } else {
          snr_db <- NA_real_
        }
      } else {
        snr_db <- NA_real_
      }
    } else {
      # Variance-based: overall variance vs residual after median filter
      win <- max(3L, as.integer(round(0.2 * sr)))
      if (win %% 2 == 0) win <- win + 1L
      med_filt <- stats::runmed(sig, k = win, endrule = "constant")
      residual <- sig - med_filt
      p_signal <- stats::var(sig)
      p_noise <- stats::var(residual)
      snr_db <- if (p_noise > 0) 10 * log10(p_signal / p_noise) else Inf
    }

    # --- Baseline wander: RMS of moving-average low-frequency component ---
    ma_win <- max(1L, as.integer(round(2 * sr)))
    kern <- rep(1 / ma_win, ma_win)
    pad_len <- (ma_win - 1L) %/% 2L
    padded_sig <- c(rev(sig[seq_len(pad_len)]), sig,
                    rev(sig[seq(n_time - pad_len + 1L, n_time)]))
    lf_padded <- as.numeric(stats::filter(padded_sig, kern, sides = 2))
    lf_component <- lf_padded[seq(pad_len + 1L, pad_len + n_time)]
    na_lf <- which(is.na(lf_component))
    if (length(na_lf) > 0) {
      valid_lf <- which(!is.na(lf_component))
      if (length(valid_lf) > 0) {
        lf_component[na_lf] <- lf_component[valid_lf[
          findInterval(na_lf, valid_lf, all.inside = TRUE)
        ]]
      } else {
        lf_component[na_lf] <- 0
      }
    }
    baseline_wander <- sqrt(mean(lf_component^2))

    # --- Saturation ratio ---
    sig_min <- min(sig)
    sig_max <- max(sig)
    sig_range <- sig_max - sig_min
    if (sig_range > 0) {
      tol <- 0.01 * sig_range
      n_saturated <- sum(sig <= sig_min + tol | sig >= sig_max - tol)
      saturation_ratio <- n_saturated / n_time
    } else {
      saturation_ratio <- 1.0
    }

    # --- Composite quality score (0-1) ---
    # SNR component: sigmoid mapping, 20 dB -> ~0.95, 0 dB -> ~0.5
    snr_score <- if (is.finite(snr_db)) {
      1 / (1 + exp(-0.2 * (snr_db - 5)))
    } else if (is.na(snr_db)) {
      0.5
    } else {
      1.0
    }

    # Baseline wander component: lower is better
    # Normalise by signal RMS
    sig_rms <- sqrt(mean(sig^2))
    wander_ratio <- if (sig_rms > 0) baseline_wander / sig_rms else 0
    wander_score <- max(0, 1 - wander_ratio)

    # Saturation component: lower is better
    sat_score <- max(0, 1 - 10 * saturation_ratio)

    quality_score <- max(0, min(1,
      0.4 * snr_score + 0.3 * wander_score + 0.3 * sat_score
    ))

    results[[ch]] <- data.frame(
      channel = as.integer(ch),
      snr_db = snr_db,
      baseline_wander = baseline_wander,
      saturation_ratio = saturation_ratio,
      quality_score = quality_score,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Correct Baseline Wander in ECG Signals
#'
#' Removes low-frequency baseline drift from ECG data using either a
#' high-pass moving-average subtraction or a running-median subtraction.
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param method Correction method: \code{"highpass"} (default) subtracts a
#'   moving average; \code{"median"} subtracts a running median.
#' @param cutoff Approximate cutoff frequency in Hz (default: 0.5).
#'   Controls the window size of the moving average or median filter.
#' @param assay_name Name of the input assay.
#'   If \code{NULL} the default assay is used.
#' @param output_assay Name of the assay in which to store the corrected
#'   signal (default: \code{"baseline_corrected"}).
#' @return A PhysioExperiment with the corrected signal stored in
#'   \code{output_assay}.
#'
#' @references Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006).
#'   \emph{Advanced Methods and Tools for ECG Data Analysis}. Artech House.
#'
#' @seealso \code{\link{ecgSignalQuality}} for signal quality assessment,
#'   \code{\link{ecgDetectRpeaks}} for R-peak detection,
#'   \code{\link{ecgQualityCheck}} for ectopic beat detection.
#'
#' @export
ecgBaselineCorrect <- function(x, method = c("highpass", "median"),
                               cutoff = 0.5,
                               assay_name = NULL,
                               output_assay = "baseline_corrected") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(is.numeric(cutoff) && length(cutoff) == 1 && cutoff > 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  corrected <- matrix(NA_real_, nrow = n_time, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    if (method == "highpass") {
      win <- max(1L, as.integer(round(sr / cutoff)))
      kern <- rep(1 / win, win)
      # Pad signal symmetrically to avoid edge NAs
      pad_len <- (win - 1L) %/% 2L
      padded <- c(rev(sig[seq_len(pad_len)]), sig,
                  rev(sig[seq(n_time - pad_len + 1L, n_time)]))
      baseline_padded <- as.numeric(stats::filter(padded, kern, sides = 2))
      baseline <- baseline_padded[seq(pad_len + 1L, pad_len + n_time)]
      # Handle any remaining NAs at boundaries
      na_idx <- which(is.na(baseline))
      if (length(na_idx) > 0) {
        valid_idx <- which(!is.na(baseline))
        if (length(valid_idx) > 0) {
          baseline[na_idx] <- baseline[valid_idx[
            findInterval(na_idx, valid_idx, all.inside = TRUE)
          ]]
        } else {
          baseline[na_idx] <- 0
        }
      }
      corrected[, ch] <- sig - baseline
    } else {
      # median method
      win <- max(3L, as.integer(round(2 * sr / cutoff)))
      if (win %% 2 == 0) win <- win + 1L
      baseline <- stats::runmed(sig, k = win, endrule = "constant")
      corrected[, ch] <- sig - baseline
    }
  }

  SummarizedExperiment::assay(x, output_assay) <- corrected
  x
}
