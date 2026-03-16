#' Create test ECG PhysioExperiment with known R-peak locations
#' @param n_time Number of time points
#' @param n_channels Number of ECG channels
#' @param sr Sampling rate in Hz
#' @param heart_rate Heart rate in bpm (default: 72)
#' @return PhysioExperiment with simulated ECG data
make_ecg <- function(n_time = 5000, n_channels = 1, sr = 500, heart_rate = 72) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)
  rr_sec <- 60 / heart_rate
  rr_samples <- as.integer(round(rr_sec * sr))

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)

  for (ch in seq_len(n_channels)) {
    signal <- rnorm(n_time, sd = 0.05)  # baseline noise

    # Add QRS-like complexes at regular intervals
    qrs_width <- as.integer(round(0.04 * sr))  # ~40ms QRS
    peak_locs <- seq(rr_samples, n_time - qrs_width, by = rr_samples)

    for (pk in peak_locs) {
      idx <- seq(max(1L, pk - qrs_width), min(n_time, pk + qrs_width))
      # Gaussian-shaped R-peak
      signal[idx] <- signal[idx] + 1.5 * exp(-((idx - pk)^2) / (2 * (qrs_width / 3)^2))
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

#' Create test ECG data with irregular beats (ectopic)
#' @param n_time Number of time points
#' @param sr Sampling rate in Hz
#' @param heart_rate Base heart rate in bpm
#' @return PhysioExperiment with simulated irregular ECG
make_ecg_irregular <- function(n_time = 5000, sr = 500, heart_rate = 72) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)
  rr_sec <- 60 / heart_rate
  rr_samples <- as.integer(round(rr_sec * sr))

  signal <- rnorm(n_time, sd = 0.05)

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
    signal[idx] <- signal[idx] + 1.5 * exp(-((idx - pk)^2) / (2 * (qrs_width / 3)^2))
  }

  data <- matrix(signal, ncol = 1)

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )
}


#' Create realistic PQRST-morphology ECG with known fiducial points
#'
#' Generates synthetic ECG with physiologically realistic P, Q, R, S, T waves.
#' Returns both the PhysioExperiment and a data.frame of known fiducial points
#' for validation testing.
#'
#' @param n_time Number of time points (default: 10000 = 20 sec at 500 Hz)
#' @param n_channels Number of ECG channels
#' @param sr Sampling rate in Hz
#' @param heart_rate Heart rate in bpm (default: 72)
#' @param noise_sd Baseline noise standard deviation (default: 0.02)
#' @return A list with components:
#'   \describe{
#'     \item{pe}{PhysioExperiment object with the ECG signal}
#'     \item{fiducials}{data.frame with columns: beat, r_peak, p_peak, q_point,
#'       s_point, t_peak, qrs_onset, qrs_offset (all in sample indices)}
#'   }
make_ecg_pqrst <- function(n_time = 10000, n_channels = 1, sr = 500,
                            heart_rate = 72, noise_sd = 0.02) {
  rr_sec <- 60 / heart_rate
  rr_samples <- as.integer(round(rr_sec * sr))

  # PQRST wave parameters (offsets in seconds relative to R-peak)
  # P wave: amplitude 0.15 mV, width ~80ms, center -160ms before R
  # Q wave: amplitude -0.10 mV, width ~30ms, center -30ms before R
  # R wave: amplitude 1.50 mV, width ~40ms, center at 0

  # S wave: amplitude -0.30 mV, width ~30ms, center +30ms after R
  # T wave: amplitude 0.30 mV, width ~160ms, center +300ms after R
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
    signal <- rnorm(n_time, sd = noise_sd)

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

        gauss <- w$amp * exp(-((idx_range - center)^2) / (2 * sigma_samples^2))
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


#' Create ECG with amplitude variation for adaptive threshold testing
#'
#' Signal amplitude gradually decreases then increases, testing whether
#' the detector adapts its threshold.
#'
#' @param n_time Number of time points
#' @param sr Sampling rate in Hz
#' @param heart_rate Heart rate in bpm
#' @return PhysioExperiment with amplitude-varying ECG
make_ecg_amplitude_varying <- function(n_time = 10000, sr = 500, heart_rate = 72) {
  rr_sec <- 60 / heart_rate
  rr_samples <- as.integer(round(rr_sec * sr))

  signal <- rnorm(n_time, sd = 0.03)
  qrs_width <- as.integer(round(0.04 * sr))

  peak_locs <- seq(rr_samples, n_time - qrs_width, by = rr_samples)
  n_peaks <- length(peak_locs)

  # Amplitude envelope: 1.5 -> 0.4 -> 1.5 (V-shape)
  amp_envelope <- 1.5 - 1.1 * (1 - abs(2 * seq_len(n_peaks) / n_peaks - 1))

  for (i in seq_along(peak_locs)) {
    pk <- peak_locs[i]
    idx <- seq(max(1L, pk - qrs_width), min(n_time, pk + qrs_width))
    signal[idx] <- signal[idx] +
      amp_envelope[i] * exp(-((idx - pk)^2) / (2 * (qrs_width / 3)^2))
  }

  PhysioExperiment(
    assays = list(raw = matrix(signal, ncol = 1)),
    colData = S4Vectors::DataFrame(label = "ECG1", type = "ECG"),
    samplingRate = sr
  )
}
