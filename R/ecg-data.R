#' Create Simulated ECG PhysioExperiment
#'
#' Generates a PhysioExperiment object containing synthetic ECG data with
#' Gaussian-shaped R-peaks at regular intervals. The resulting object is
#' suitable for testing R-peak detection, RR interval computation, and HRV
#' analysis pipelines.
#'
#' @param n_time Number of time points (default: 5000).
#' @param n_channels Number of ECG channels (default: 1).
#' @param sr Sampling rate in Hz (default: 500).
#' @param heart_rate Heart rate in beats per minute (default: 72).
#' @return A PhysioExperiment object with a single \code{"raw"} assay
#'   containing the simulated ECG signal.
#'
#' @references Pan, J. & Tompkins, W.J. (1985). "A real-time QRS detection
#'   algorithm." \emph{IEEE Transactions on Biomedical Engineering}, 32(3),
#'   230--236. \doi{10.1109/TBME.1985.325532}
#'
#'   Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006).
#'   \emph{Advanced Methods and Tools for ECG Data Analysis}. Artech House.
#'
#' @seealso \code{\link{make_ecg_irregular}} for ECG with ectopic beats,
#'   \code{\link{make_ecg_pqrst}} for ECG with full PQRST morphology,
#'   \code{\link{make_ecg_noisy}} for ECG with noise artifacts,
#'   \code{\link{ecgDetectRpeaks}} for R-peak detection.
#'
#' @export
#' @examples
#' pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
#' pe
make_ecg <- function(n_time = 5000, n_channels = 1, sr = 500, heart_rate = 72) {
  rr_sec <- 60 / heart_rate
  rr_samples <- as.integer(round(rr_sec * sr))

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    signal <- stats::rnorm(n_time, sd = 0.05)  # baseline noise

    # Add QRS-like complexes at regular intervals
    qrs_width <- as.integer(round(0.04 * sr))  # ~40ms QRS
    peak_locs <- seq(rr_samples, n_time - qrs_width, by = rr_samples)

    for (pk in peak_locs) {
      idx <- seq(max(1L, pk - qrs_width), min(n_time, pk + qrs_width))
      # Gaussian-shaped R-peak
      signal[idx] <- signal[idx] +
        1.5 * exp(-((idx - pk)^2) / (2 * (qrs_width / 3)^2))
    }

    data[, ch] <- signal
  }

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0("ECG", seq_len(n_channels)),
      type = rep("ECG", n_channels)
    ),
    samplingRate = sr
  )
}


#' Create Simulated ECG with Irregular R-R Intervals
#'
#' Generates a PhysioExperiment object containing synthetic ECG data with
#' irregular beat timing. Every 5th beat is premature (60\% of normal RR
#' interval), followed by a compensatory pause (140\% of normal RR interval).
#' Useful for testing ectopic beat detection and RR interval correction.
#'
#' @param n_time Number of time points (default: 5000).
#' @param sr Sampling rate in Hz (default: 500).
#' @param heart_rate Base heart rate in beats per minute (default: 72).
#' @return A PhysioExperiment object with a single \code{"raw"} assay
#'   containing the simulated ECG signal with irregular beats.
#'
#' @references Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006).
#'   \emph{Advanced Methods and Tools for ECG Data Analysis}. Artech House.
#'
#'   Task Force of the European Society of Cardiology and the North American
#'   Society of Pacing and Electrophysiology (1996). "Heart rate variability:
#'   Standards of measurement, physiological interpretation and clinical use."
#'   \emph{Circulation}, 93(5), 1043--1065.
#'
#' @seealso \code{\link{make_ecg}} for regular ECG data,
#'   \code{\link{ecgQualityCheck}} for ectopic beat detection,
#'   \code{\link{ecgRRcorrect}} for ectopic beat correction,
#'   \code{\link{ecgRRintervals}} for RR interval computation.
#'
#' @export
#' @examples
#' pe <- make_ecg_irregular(n_time = 5000, sr = 500, heart_rate = 72)
#' pe
make_ecg_irregular <- function(n_time = 5000, sr = 500, heart_rate = 72) {
  rr_sec <- 60 / heart_rate
  rr_samples <- as.integer(round(rr_sec * sr))

  signal <- stats::rnorm(n_time, sd = 0.05)
  qrs_width <- as.integer(round(0.04 * sr))

  # Build peak locations with some irregular intervals
  peak_locs <- integer(0)
  pos <- rr_samples
  beat_idx <- 0L
  while (pos < n_time - qrs_width) {
    beat_idx <- beat_idx + 1L
    peak_locs <- c(peak_locs, pos)

    # Every 5th beat: premature (short RR) followed by compensatory (long RR)
    if (beat_idx %% 5 == 0) {
      pos <- pos + as.integer(round(rr_samples * 0.6))  # premature
    } else if (beat_idx %% 5 == 1 && beat_idx > 1) {
      pos <- pos + as.integer(round(rr_samples * 1.4))  # compensatory
    } else {
      pos <- pos + rr_samples
    }
  }

  for (pk in peak_locs) {
    idx <- seq(max(1L, pk - qrs_width), min(n_time, pk + qrs_width))
    signal[idx] <- signal[idx] +
      1.5 * exp(-((idx - pk)^2) / (2 * (qrs_width / 3)^2))
  }

  data <- matrix(signal, ncol = 1)

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )
}


#' Create Simulated ECG with PQRST Morphology
#'
#' Generates a PhysioExperiment object containing synthetic ECG data with
#' physiologically realistic P, Q, R, S, and T wave components. Returns both
#' the PhysioExperiment and a data.frame of known fiducial points for
#' validation testing of waveform delineation algorithms.
#'
#' @param n_time Number of time points (default: 10000, i.e., 20 seconds at
#'   500 Hz).
#' @param n_channels Number of ECG channels (default: 1).
#' @param sr Sampling rate in Hz (default: 500).
#' @param heart_rate Heart rate in beats per minute (default: 72).
#' @param noise_sd Standard deviation of baseline noise (default: 0.02).
#' @return A list with two components:
#'   \describe{
#'     \item{pe}{PhysioExperiment object with the simulated ECG signal.}
#'     \item{fiducials}{A data.frame with columns: \code{beat}, \code{r_peak},
#'       \code{p_peak}, \code{q_point}, \code{s_point}, \code{t_peak},
#'       \code{qrs_onset}, \code{qrs_offset} (all in sample indices).}
#'   }
#'
#' @references Goldberger, A.L., et al. (2000). "PhysioBank, PhysioToolkit,
#'   and PhysioNet: Components of a new research resource for complex
#'   physiologic signals." \emph{Circulation}, 101(23), e215--e220.
#'   \doi{10.1161/01.CIR.101.23.e215}
#'
#'   Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006).
#'   \emph{Advanced Methods and Tools for ECG Data Analysis}. Artech House.
#'
#' @seealso \code{\link{make_ecg}} for basic ECG data,
#'   \code{\link{ecgDelineate}} for PQRST waveform delineation,
#'   \code{\link{ecgIntervals}} for computing clinical ECG intervals,
#'   \code{\link{ecgDetectRpeaks}} for R-peak detection.
#'
#' @export
#' @examples
#' result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72)
#' pe <- result$pe
#' fiducials <- result$fiducials
#' head(fiducials)
make_ecg_pqrst <- function(n_time = 10000, n_channels = 1, sr = 500,
                            heart_rate = 72, noise_sd = 0.02) {
  rr_sec <- 60 / heart_rate
  rr_samples <- as.integer(round(rr_sec * sr))

  # PQRST wave parameters (offsets in seconds relative to R-peak)
  waves <- list(
    P = list(amp = 0.15, sigma_ms = 20, offset_ms = -160),
    Q = list(amp = -0.10, sigma_ms = 8,  offset_ms = -30),
    R = list(amp = 1.50,  sigma_ms = 7,  offset_ms = 0),
    S = list(amp = -0.30, sigma_ms = 8,  offset_ms = 30),
    T = list(amp = 0.30,  sigma_ms = 40, offset_ms = 300)
  )

  # Compute peak locations (leaving margin for P wave before first beat)
  margin <- as.integer(round(0.3 * sr))  # 300ms margin
  peak_locs <- seq(margin + rr_samples, n_time - margin, by = rr_samples)

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  fiducials_list <- list()

  for (ch in seq_len(n_channels)) {
    signal <- stats::rnorm(n_time, sd = noise_sd)

    for (beat_idx in seq_along(peak_locs)) {
      r_pos <- peak_locs[beat_idx]
      fid <- list(beat = beat_idx)

      for (wname in names(waves)) {
        w <- waves[[wname]]
        center <- r_pos + as.integer(round(w$offset_ms / 1000 * sr))
        sigma_samples <- w$sigma_ms / 1000 * sr

        # Compute contribution over +/- 4 sigma
        half_win <- as.integer(ceiling(4 * sigma_samples))
        idx_range <- seq(max(1L, center - half_win),
                         min(n_time, center + half_win))

        gauss <- w$amp * exp(-((idx_range - center)^2) /
                               (2 * sigma_samples^2))
        signal[idx_range] <- signal[idx_range] + gauss

        # Record fiducial point (peak/trough location)
        if (center >= 1 && center <= n_time) {
          fid[[tolower(paste0(wname, "_peak"))]] <- center
        }
      }

      # Derive QRS onset/offset from Q and S positions
      q_center <- r_pos + as.integer(round(waves$Q$offset_ms / 1000 * sr))
      s_center <- r_pos + as.integer(round(waves$S$offset_ms / 1000 * sr))
      q_sigma <- waves$Q$sigma_ms / 1000 * sr
      s_sigma <- waves$S$sigma_ms / 1000 * sr

      fid$r_peak <- r_pos
      fid$q_point <- q_center
      fid$s_point <- s_center
      fid$qrs_onset <- as.integer(round(q_center - 2 * q_sigma))
      fid$qrs_offset <- as.integer(round(s_center + 2 * s_sigma))

      if (ch == 1) {
        fiducials_list[[beat_idx]] <- fid
      }
    }

    data[, ch] <- signal
  }

  fiducials <- do.call(rbind, lapply(fiducials_list, as.data.frame))

  pe <- PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0("ECG", seq_len(n_channels)),
      type = rep("ECG", n_channels)
    ),
    samplingRate = sr
  )

  list(pe = pe, fiducials = fiducials)
}


#' Create Simulated ECG with Noise Artifacts
#'
#' Generates a PhysioExperiment object containing synthetic ECG data
#' contaminated with multiple noise sources: baseline wander (0.3 Hz
#' sinusoidal drift), powerline interference (50 Hz), and broadband
#' Gaussian noise. Useful for testing signal quality assessment, baseline
#' correction, and filtering pipelines.
#'
#' @param n_time Number of time points (default: 5000).
#' @param n_channels Number of ECG channels (default: 1).
#' @param sr Sampling rate in Hz (default: 500).
#' @param heart_rate Heart rate in beats per minute (default: 72).
#' @param baseline_amp Amplitude of baseline wander in arbitrary units
#'   (default: 0.3).
#' @param powerline_amp Amplitude of 50 Hz powerline noise (default: 0.1).
#' @param noise_sd Standard deviation of broadband Gaussian noise
#'   (default: 0.15).
#' @return A PhysioExperiment object with a single \code{"raw"} assay
#'   containing the noisy ECG signal.
#'
#' @references Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006).
#'   \emph{Advanced Methods and Tools for ECG Data Analysis}. Artech House.
#'
#'   Shaffer, F. & Ginsberg, J.P. (2017). "An overview of heart rate
#'   variability metrics and norms." \emph{Frontiers in Public Health}, 5, 258.
#'   \doi{10.3389/fpubh.2017.00258}
#'
#' @seealso \code{\link{make_ecg}} for clean ECG data,
#'   \code{\link{ecgSignalQuality}} for signal quality assessment,
#'   \code{\link{ecgBaselineCorrect}} for baseline wander removal,
#'   \code{\link{ecgDetectRpeaks}} for R-peak detection.
#'
#' @export
#' @examples
#' pe <- make_ecg_noisy(n_time = 5000, sr = 500)
#' pe
make_ecg_noisy <- function(n_time = 5000, n_channels = 1, sr = 500,
                            heart_rate = 72, baseline_amp = 0.3,
                            powerline_amp = 0.1, noise_sd = 0.15) {
  rr_sec <- 60 / heart_rate
  rr_samples <- as.integer(round(rr_sec * sr))
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    # Start with broadband noise
    signal <- stats::rnorm(n_time, sd = noise_sd)

    # Add QRS-like complexes at regular intervals
    qrs_width <- as.integer(round(0.04 * sr))  # ~40ms QRS
    peak_locs <- seq(rr_samples, n_time - qrs_width, by = rr_samples)

    for (pk in peak_locs) {
      idx <- seq(max(1L, pk - qrs_width), min(n_time, pk + qrs_width))
      signal[idx] <- signal[idx] +
        1.5 * exp(-((idx - pk)^2) / (2 * (qrs_width / 3)^2))
    }

    # Add baseline wander (low-frequency sinusoidal drift)
    signal <- signal + baseline_amp * sin(2 * pi * 0.3 * t)

    # Add powerline interference (50 Hz)
    signal <- signal + powerline_amp * sin(2 * pi * 50 * t)

    data[, ch] <- signal
  }

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0("ECG", seq_len(n_channels)),
      type = rep("ECG", n_channels)
    ),
    samplingRate = sr
  )
}
