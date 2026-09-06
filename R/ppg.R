# Photoplethysmography analysis

#' Validate a PPG scalar
#' @keywords internal
#' @noRd
.ppg_scalar <- function(x, name, lower = 0, inclusive = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (valid) {
    valid <- if (inclusive) x >= lower else x > lower
  }
  if (!valid) {
    relation <- if (inclusive) "at least" else "greater than"
    stop(sprintf("%s must be a finite number %s %s.",
                 name, relation, format(lower)), call. = FALSE)
  }
  as.numeric(x)
}


#' Resolve PPG assay data
#' @keywords internal
#' @noRd
.ppg_data <- function(x, assay_name = NULL) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("x must be a PhysioExperiment object.", call. = FALSE)
  }
  if (is.null(assay_name)) {
    assay_name <- defaultAssay(x)
  }
  data <- SummarizedExperiment::assay(x, assay_name)
  if (!is.matrix(data) || !is.numeric(data) ||
      nrow(data) < 1L || ncol(data) < 1L ||
      any(!is.finite(data))) {
    stop("The PPG assay must be a non-empty finite numeric matrix.",
         call. = FALSE)
  }
  sr <- .ppg_scalar(samplingRate(x), "sampling rate")
  list(data = data, sr = sr, assay_name = assay_name)
}


#' Validate a shared peak table
#' @keywords internal
#' @noRd
.ppg_peaks <- function(peaks,
                       n_channels,
                       n_time = NULL,
                       require_time = FALSE) {
  required <- c("channel", "sample")
  if (require_time) {
    required <- c(required, "time_sec")
  }
  if (!is.data.frame(peaks) || !all(required %in% names(peaks))) {
    stop(sprintf(
      "peaks must be a data frame with columns %s.",
      paste(required, collapse = ", ")
    ), call. = FALSE)
  }
  if (nrow(peaks) == 0L) {
    return(peaks)
  }
  if (!is.numeric(peaks$channel) || !is.numeric(peaks$sample) ||
      any(!is.finite(peaks$channel)) ||
      any(!is.finite(peaks$sample)) ||
      any(peaks$channel != as.integer(peaks$channel)) ||
      any(peaks$sample != as.integer(peaks$sample)) ||
      any(peaks$channel < 1L | peaks$channel > n_channels) ||
      any(peaks$sample < 1L)) {
    stop("peaks contain invalid channel or sample indices.",
         call. = FALSE)
  }
  if (!is.null(n_time) && any(peaks$sample > n_time)) {
    stop("peak samples lie outside the PPG assay.", call. = FALSE)
  }
  if (require_time &&
      (!is.numeric(peaks$time_sec) || any(!is.finite(peaks$time_sec)))) {
    stop("peaks$time_sec must contain finite numeric values.",
         call. = FALSE)
  }
  split_samples <- split(peaks$sample, peaks$channel)
  if (any(vapply(split_samples, anyDuplicated, integer(1)) > 0L)) {
    stop("peaks must not contain duplicate samples within a channel.",
         call. = FALSE)
  }
  peaks$channel <- as.integer(peaks$channel)
  peaks$sample <- as.integer(peaks$sample)
  peaks
}


#' Empty PPG peak table
#' @keywords internal
#' @noRd
.ppg_empty_peaks <- function() {
  data.frame(
    channel = integer(),
    sample = integer(),
    time_sec = numeric(),
    amplitude = numeric()
  )
}


#' Detect systolic peaks in PPG
#'
#' Implements the Elgendi two event-related moving-average (TERMA) detector.
#' The shipped windowed-sinc FIR is used in place of the original Butterworth
#' filter. Its linear-phase group delay is explicitly undone by snapping each
#' event-block candidate to the corresponding raw-signal maximum.
#'
#' @details
#' This implementation uses the package's pure-R windowed-sinc FIR rather than
#' the Butterworth filter in the original method. A strong late dicrotic wave
#' may require increasing `refractory_ms` so that it is not counted as a second
#' systolic pulse.
#'
#' @param x A PPG `PhysioExperiment`.
#' @param w1_ms Systolic event window in milliseconds.
#' @param w2_ms Beat window in milliseconds.
#' @param beta Moving-average threshold offset.
#' @param low,high PPG band-pass cutoffs in Hz.
#' @param order FIR order; even values are promoted to the next odd value.
#' @param refractory_ms Minimum systolic-peak separation in milliseconds.
#' @param assay_name Optional raw PPG assay name.
#'
#' @return A data frame with `channel`, `sample`, `time_sec`, and raw
#'   `amplitude`.
#'
#' @references
#' Elgendi M, Norton I, Brearley M, Abbott D, Schuurmans D (2013).
#' Systolic peak detection in acceleration photoplethysmograms measured from
#' emergency responders in tropical conditions. *PLoS ONE*, 8:e76585.
#' \doi{10.1371/journal.pone.0076585}
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' simulation <- make_ppg(n_time = 3000)
#' head(ppgDetectPulses(simulation$pe))
ppgDetectPulses <- function(x,
                            w1_ms = 111,
                            w2_ms = 667,
                            beta = 0.02,
                            low = 0.5,
                            high = 8,
                            order = 65,
                            refractory_ms = 300,
                            assay_name = NULL) {
  source <- .ppg_data(x, assay_name)
  sr <- source$sr
  w1_ms <- .ppg_scalar(w1_ms, "w1_ms")
  w2_ms <- .ppg_scalar(w2_ms, "w2_ms")
  beta <- .ppg_scalar(beta, "beta", 0, inclusive = TRUE)
  low <- .ppg_scalar(low, "low")
  high <- .ppg_scalar(high, "high")
  refractory_ms <- .ppg_scalar(refractory_ms, "refractory_ms")
  if (high <= low || high >= sr / 2) {
    stop("low and high must be increasing frequencies below Nyquist.",
         call. = FALSE)
  }
  if (!is.numeric(order) || length(order) != 1L ||
      !is.finite(order) || order < 3L ||
      order != as.integer(order)) {
    stop("order must be an integer of at least 3.", call. = FALSE)
  }
  order <- as.integer(order)
  order_odd <- if (order %% 2L == 0L) order + 1L else order
  group_delay <- (order_odd - 1L) %/% 2L
  w1 <- max(1L, as.integer(round(w1_ms / 1000 * sr)))
  w2 <- max(1L, as.integer(round(w2_ms / 1000 * sr)))
  refractory <- max(
    1L, as.integer(round(refractory_ms / 1000 * sr))
  )
  if (nrow(source$data) < max(order_odd, w1, w2)) {
    return(.ppg_empty_peaks())
  }

  output <- list()
  for (channel in seq_len(ncol(source$data))) {
    signal <- source$data[, channel]
    if (max(signal) == min(signal)) {
      next
    }
    filtered <- .fir_bandpass(
      signal - mean(signal),
      sr,
      low = low,
      high = high,
      order = order
    )
    rectified <- pmax(filtered, 0)^2
    moving_peak <- as.numeric(stats::filter(
      rectified, rep(1 / w1, w1), sides = 2
    ))
    moving_beat <- as.numeric(stats::filter(
      rectified, rep(1 / w2, w2), sides = 2
    ))
    moving_peak[is.na(moving_peak)] <- 0
    moving_beat[is.na(moving_beat)] <- 0
    threshold <- moving_beat + beta * mean(rectified)
    blocks <- moving_peak > threshold
    runs <- rle(blocks)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L

    candidates <- integer()
    for (i in seq_along(runs$values)) {
      if (isTRUE(runs$values[i]) && runs$lengths[i] >= w1) {
        lo <- starts[i]
        hi <- ends[i]
        candidates <- c(
          candidates,
          lo + which.max(filtered[lo:hi]) - 1L
        )
      }
    }
    if (length(candidates) == 0L) {
      next
    }

    peaks <- vapply(candidates, function(candidate) {
      lo <- max(1L, candidate - group_delay - w1)
      hi <- min(length(signal), candidate - group_delay + w1)
      if (hi < lo) {
        return(NA_integer_)
      }
      as.integer(lo + which.max(signal[lo:hi]) - 1L)
    }, integer(1))
    peaks <- peaks[!is.na(peaks)]
    peaks <- .ecg_apply_refractory(
      sort(unique(peaks)), refractory
    )
    if (length(peaks) == 0L) {
      next
    }
    output[[length(output) + 1L]] <- data.frame(
      channel = rep(as.integer(channel), length(peaks)),
      sample = as.integer(peaks),
      time_sec = (peaks - 1) / sr,
      amplitude = signal[peaks]
    )
  }

  if (length(output) == 0L) {
    .ppg_empty_peaks()
  } else {
    do.call(rbind, output)
  }
}


#' Assess PPG signal quality
#'
#' Computes population skewness and raw kurtosis, either per channel or in
#' non-overlapping windows.
#'
#' @param x A PPG `PhysioExperiment`.
#' @param window_sec Optional non-overlapping window duration in seconds.
#' @param skew_threshold Minimum skewness SQI accepted as usable.
#' @param assay_name Optional raw PPG assay name.
#'
#' @return A data frame with channel, segment, start time, skewness SQI,
#'   kurtosis SQI, and acceptance flag.
#'
#' @references
#' Elgendi M (2016). Optimal signal quality index for photoplethysmogram
#' signals. *Bioengineering*, 3:21. \doi{10.3390/bioengineering3040021}
#'
#' @export
#'
#' @examples
#' simulation <- make_ppg(n_time = 2000)
#' ppgQuality(simulation$pe)
ppgQuality <- function(x,
                       window_sec = NULL,
                       skew_threshold = -0.3,
                       assay_name = NULL) {
  source <- .ppg_data(x, assay_name)
  skew_threshold <- .ppg_scalar(
    skew_threshold, "skew_threshold", -Inf, inclusive = TRUE
  )
  if (is.null(window_sec)) {
    window_samples <- nrow(source$data)
  } else {
    window_sec <- .ppg_scalar(window_sec, "window_sec")
    window_samples <- max(
      1L, as.integer(round(window_sec * source$sr))
    )
  }

  rows <- list()
  starts <- seq.int(1L, nrow(source$data), by = window_samples)
  for (channel in seq_len(ncol(source$data))) {
    for (segment in seq_along(starts)) {
      start <- starts[segment]
      end <- min(nrow(source$data), start + window_samples - 1L)
      values <- source$data[start:end, channel]
      centered <- values - mean(values)
      sigma <- sqrt(mean(centered^2))
      if (!is.finite(sigma) || sigma == 0) {
        skewness <- NA_real_
        kurtosis <- NA_real_
        accept <- FALSE
      } else {
        standardized <- centered / sigma
        skewness <- mean(standardized^3)
        kurtosis <- mean(standardized^4)
        accept <- is.finite(skewness) && skewness >= skew_threshold
      }
      rows[[length(rows) + 1L]] <- data.frame(
        channel = as.integer(channel),
        segment = as.integer(segment),
        start_sec = (start - 1) / source$sr,
        s_sqi = skewness,
        k_sqi = kurtosis,
        accept = accept
      )
    }
  }
  do.call(rbind, rows)
}


#' Calculate pulse-rate variability
#'
#' Converts detected PPG pulses to pulse-to-pulse intervals and delegates
#' time-, frequency-, and optional nonlinear-domain metrics to the existing HRV
#' engine with ECG rhythm gating disabled.
#'
#' @param x A PPG `PhysioExperiment`.
#' @param peaks Optional peak table from [ppgDetectPulses()].
#' @param freq Logical; calculate frequency-domain PRV.
#' @param nonlinear Logical; calculate nonlinear PRV.
#' @param assay_name Optional raw PPG assay name.
#' @param ... Detection arguments passed to [ppgDetectPulses()] when `peaks`
#'   is `NULL`.
#'
#' @return A named list containing `ppi`, `time`, `freq`, and `nonlinear`.
#'
#' @references
#' Schafer A, Vagedes J (2013). How accurate is pulse rate variability as an
#' estimate of heart rate variability? *International Journal of Cardiology*,
#' 166:15-29. \doi{10.1016/j.ijcard.2012.03.119}
#'
#' @export
#'
#' @examples
#' simulation <- make_ppg(n_time = 15000)
#' pulseRateVariability(simulation$pe, freq = FALSE)
pulseRateVariability <- function(x,
                                 peaks = NULL,
                                 freq = TRUE,
                                 nonlinear = FALSE,
                                 assay_name = NULL,
                                 ...) {
  source <- .ppg_data(x, assay_name)
  if (!is.logical(freq) || length(freq) != 1L || is.na(freq) ||
      !is.logical(nonlinear) || length(nonlinear) != 1L ||
      is.na(nonlinear)) {
    stop("freq and nonlinear must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(peaks)) {
    peaks <- ppgDetectPulses(
      x, assay_name = assay_name, ...
    )
  }
  peaks <- .ppg_peaks(
    peaks,
    n_channels = ncol(source$data),
    require_time = TRUE
  )
  intervals <- ecgRRintervals(x, peaks)
  names(intervals)[names(intervals) == "rr_ms"] <- "ppi_ms"

  empty_time <- data.frame(
    channel = integer(),
    mean_pp = numeric(),
    sdpp = numeric(),
    rmssd = numeric(),
    pnn50 = numeric(),
    mean_pr = numeric()
  )
  if (nrow(intervals) == 0L) {
    return(list(
      ppi = intervals,
      time = empty_time,
      freq = NULL,
      nonlinear = NULL
    ))
  }

  rr <- data.frame(
    channel = intervals$channel,
    rr_ms = intervals$ppi_ms,
    time_sec = intervals$time_sec
  )
  time_metrics <- ecgHRVtime(rr, rhythm_check = FALSE)
  names(time_metrics)[names(time_metrics) == "mean_rr"] <- "mean_pp"
  names(time_metrics)[names(time_metrics) == "sdnn"] <- "sdpp"
  names(time_metrics)[names(time_metrics) == "mean_hr"] <- "mean_pr"
  frequency_metrics <- if (freq) {
    ecgHRVfreq(rr, rhythm_check = FALSE)
  } else {
    NULL
  }
  nonlinear_metrics <- if (nonlinear) {
    ecgHRVnonlinear(rr, rhythm_check = FALSE)
  } else {
    NULL
  }

  list(
    ppi = intervals,
    time = time_metrics,
    freq = frequency_metrics,
    nonlinear = nonlinear_metrics
  )
}


#' Locate pulse feet
#' @keywords internal
#' @noRd
.ppg_feet <- function(signal, systolic, sr) {
  systolic <- sort(as.integer(systolic))
  feet <- integer(length(systolic))
  back <- max(1L, as.integer(round(sr)))
  for (i in seq_along(systolic)) {
    lo <- if (i == 1L) {
      max(1L, systolic[i] - back)
    } else {
      systolic[i - 1L] + 1L
    }
    hi <- systolic[i]
    feet[i] <- lo + which.min(signal[lo:hi]) - 1L
  }
  feet
}


#' Calculate the PPG perfusion index
#'
#' Perfusion index is 100 times the pulsatile systolic-to-foot amplitude
#' divided by the raw-signal DC level.
#'
#' @param x A raw, DC-preserving PPG `PhysioExperiment`.
#' @param peaks Optional systolic peak table from [ppgDetectPulses()].
#' @param assay_name Optional raw PPG assay name.
#' @param agg Aggregation of per-pulse AC amplitudes.
#'
#' @return A per-channel data frame containing AC, DC, perfusion index, and
#'   pulse count.
#'
#' @references
#' Charlton PH, Kyriacou PA, Mant J, et al. (2022). Wearable
#' photoplethysmography for cardiovascular monitoring. *Proceedings of the
#' IEEE*, 110:355-381. \doi{10.1109/JPROC.2022.3149785}
#'
#' @export
#'
#' @examples
#' simulation <- make_ppg(n_time = 5000)
#' perfusionIndex(simulation$pe)
perfusionIndex <- function(x,
                           peaks = NULL,
                           assay_name = NULL,
                           agg = c("median", "mean")) {
  source <- .ppg_data(x, assay_name)
  agg <- match.arg(agg)
  if (is.null(peaks)) {
    peaks <- ppgDetectPulses(x, assay_name = assay_name)
  }
  peaks <- .ppg_peaks(
    peaks, n_channels = ncol(source$data), n_time = nrow(source$data)
  )

  rows <- lapply(seq_len(ncol(source$data)), function(channel) {
    signal <- source$data[, channel]
    selected <- sort(peaks$sample[peaks$channel == channel])
    dc <- mean(signal)
    if (length(selected) == 0L) {
      ac <- NA_real_
    } else {
      feet <- .ppg_feet(signal, selected, source$sr)
      amplitudes <- signal[selected] - signal[feet]
      ac <- if (agg == "median") {
        stats::median(amplitudes)
      } else {
        mean(amplitudes)
      }
    }
    dc_tolerance <- sqrt(.Machine$double.eps) *
      max(1, max(abs(signal)))
    pi_pct <- if (length(selected) < 2L ||
                  !is.finite(ac) || abs(dc) <= dc_tolerance) {
      NA_real_
    } else {
      100 * ac / dc
    }
    data.frame(
      channel = as.integer(channel),
      ac = ac,
      dc = dc,
      pi_pct = pi_pct,
      n_pulses = as.integer(length(selected))
    )
  })
  do.call(rbind, rows)
}


#' Find a dicrotic notch and diastolic peak
#' @keywords internal
#' @noRd
.ppg_diastolic_fiducials <- function(signal, systolic, next_foot) {
  if (is.na(next_foot) || next_foot - systolic < 4L) {
    return(c(notch = NA_integer_, diastolic = NA_integer_))
  }
  candidates <- seq.int(systolic + 1L, next_foot - 1L)
  candidates <- candidates[
    candidates > 1L & candidates < length(signal)
  ]
  if (length(candidates) < 2L) {
    return(c(notch = NA_integer_, diastolic = NA_integer_))
  }
  local_minimum <- candidates[
    signal[candidates] <= signal[candidates - 1L] &
      signal[candidates] < signal[candidates + 1L]
  ]
  if (length(local_minimum) == 0L) {
    return(c(notch = NA_integer_, diastolic = NA_integer_))
  }
  notch <- local_minimum[1L]
  diastolic_range <- seq.int(notch + 1L, next_foot - 1L)
  if (length(diastolic_range) == 0L) {
    return(c(notch = notch, diastolic = NA_integer_))
  }
  diastolic <- diastolic_range[
    which.max(signal[diastolic_range])
  ]
  c(notch = as.integer(notch), diastolic = as.integer(diastolic))
}


#' Extract PPG pulse-wave features
#'
#' Measures foot-to-systolic rise time, pulse amplitude, dicrotic and
#' diastolic fiducials, reflection and augmentation indices, pulse width, and
#' optional ECG-to-PPG pulse transit time and pulse-wave velocity.
#'
#' @param x A raw PPG `PhysioExperiment`.
#' @param peaks Optional systolic peak table from [ppgDetectPulses()].
#' @param ecg_peaks Optional ECG R-peak table containing `sample` and
#'   optionally `channel`.
#' @param distance_m Optional ECG-to-PPG path length in metres for PWV.
#' @param assay_name Optional raw PPG assay name.
#'
#' @return A data frame with one row per pulse and waveform/PTT/PWV features.
#'
#' @references
#' Elgendi M (2012). On the analysis of fingertip photoplethysmogram signals.
#' *Current Cardiology Reviews*, 8:14-25.
#' \doi{10.2174/157340312801215782}
#'
#' @export
#'
#' @examples
#' simulation <- make_ppg(n_time = 5000)
#' head(pulseWaveFeatures(simulation$pe))
pulseWaveFeatures <- function(x,
                              peaks = NULL,
                              ecg_peaks = NULL,
                              distance_m = NULL,
                              assay_name = NULL) {
  source <- .ppg_data(x, assay_name)
  if (is.null(peaks)) {
    peaks <- ppgDetectPulses(x, assay_name = assay_name)
  }
  peaks <- .ppg_peaks(
    peaks, n_channels = ncol(source$data), n_time = nrow(source$data)
  )
  if (!is.null(distance_m)) {
    distance_m <- .ppg_scalar(distance_m, "distance_m")
  }
  if (!is.null(ecg_peaks)) {
    if (!is.data.frame(ecg_peaks) || !"sample" %in% names(ecg_peaks) ||
        !is.numeric(ecg_peaks$sample) ||
        any(!is.finite(ecg_peaks$sample)) ||
        any(ecg_peaks$sample != as.integer(ecg_peaks$sample)) ||
        any(ecg_peaks$sample < 1L)) {
      stop("ecg_peaks must contain positive finite integer sample indices.",
           call. = FALSE)
    }
    if ("channel" %in% names(ecg_peaks) &&
        (!is.numeric(ecg_peaks$channel) ||
         any(!is.finite(ecg_peaks$channel)) ||
         any(ecg_peaks$channel != as.integer(ecg_peaks$channel)) ||
         any(ecg_peaks$channel < 1L))) {
      stop("ecg_peaks$channel must contain positive integer indices.",
           call. = FALSE)
    }
  }

  rows <- list()
  ptt_window <- max(1L, as.integer(round(0.4 * source$sr)))
  for (channel in seq_len(ncol(source$data))) {
    signal <- source$data[, channel]
    systolic <- sort(peaks$sample[peaks$channel == channel])
    if (length(systolic) == 0L) {
      next
    }
    feet <- .ppg_feet(signal, systolic, source$sr)
    next_feet <- c(feet[-1L], NA_integer_)
    for (beat in seq_along(systolic)) {
      amplitude <- signal[systolic[beat]] - signal[feet[beat]]
      fiducials <- .ppg_diastolic_fiducials(
        signal, systolic[beat], next_feet[beat]
      )
      if (is.finite(amplitude) && amplitude > 0 &&
          !is.na(fiducials["diastolic"])) {
        diastolic_amplitude <- signal[fiducials["diastolic"]] -
          signal[feet[beat]]
        reflection_index <- 100 * diastolic_amplitude / amplitude
        augmentation_index <- 100 *
          (signal[systolic[beat]] -
             signal[fiducials["diastolic"]]) / amplitude
      } else {
        reflection_index <- NA_real_
        augmentation_index <- NA_real_
      }

      ptt_ms <- NA_real_
      if (!is.null(ecg_peaks)) {
        candidates <- ecg_peaks
        if ("channel" %in% names(candidates)) {
          candidates <- candidates[candidates$channel == channel, ,
                                   drop = FALSE]
        }
        preceding <- candidates$sample[
          candidates$sample < feet[beat] &
            candidates$sample >= feet[beat] - ptt_window
        ]
        if (length(preceding) > 0L) {
          r_peak <- max(preceding)
          ptt_ms <- (feet[beat] - r_peak) / source$sr * 1000
        }
      }
      pwv <- if (!is.null(distance_m) &&
                 is.finite(ptt_ms) && ptt_ms > 0) {
        distance_m / (ptt_ms / 1000)
      } else {
        NA_real_
      }

      rows[[length(rows) + 1L]] <- data.frame(
        channel = as.integer(channel),
        beat = as.integer(beat),
        foot_sample = as.integer(feet[beat]),
        systolic_sample = as.integer(systolic[beat]),
        systolic_time_sec = (systolic[beat] - 1) / source$sr,
        rise_time_ms = (systolic[beat] - feet[beat]) /
          source$sr * 1000,
        pulse_amplitude = amplitude,
        dicrotic_notch_sample = as.integer(fiducials["notch"]),
        diastolic_sample = as.integer(fiducials["diastolic"]),
        reflection_index = reflection_index,
        augmentation_index = augmentation_index,
        pulse_width_ms = if (is.na(next_feet[beat])) {
          NA_real_
        } else {
          (next_feet[beat] - feet[beat]) / source$sr * 1000
        },
        ptt_ms = ptt_ms,
        pwv_m_s = pwv
      )
    }
  }

  if (length(rows) == 0L) {
    data.frame(
      channel = integer(),
      beat = integer(),
      foot_sample = integer(),
      systolic_sample = integer(),
      systolic_time_sec = numeric(),
      rise_time_ms = numeric(),
      pulse_amplitude = numeric(),
      dicrotic_notch_sample = integer(),
      diastolic_sample = integer(),
      reflection_index = numeric(),
      augmentation_index = numeric(),
      pulse_width_ms = numeric(),
      ptt_ms = numeric(),
      pwv_m_s = numeric()
    )
  } else {
    do.call(rbind, rows)
  }
}


#' Estimate respiration from PPG modulation
#'
#' Extracts amplitude modulation, baseline wander, or frequency modulation
#' from successive PPG pulses, resamples the feature to a uniform grid, and
#' estimates its respiratory-band spectral peak.
#'
#' @param x A raw PPG `PhysioExperiment`.
#' @param peaks Optional systolic peak table from [ppgDetectPulses()].
#' @param method `"am"` for pulse amplitude, `"bw"` for raw systolic level, or
#'   `"fm"` for instantaneous pulse rate.
#' @param fs Uniform feature resampling rate in Hz.
#' @param resp_band Two-element respiratory band in Hz.
#' @param assay_name Optional raw PPG assay name.
#'
#' @return A list containing method, per-pulse features, uniformly resampled
#'   EDR-like traces, and respiratory rates.
#'
#' @references
#' Charlton PH, Bonnici T, Mainardi L, et al. (2016). An assessment of
#' algorithms to estimate respiratory rate from the electrocardiogram and
#' photoplethysmogram. *Physiological Measurement*, 37:610-626.
#' \doi{10.1088/0967-3334/37/4/610}
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' simulation <- make_ppg(n_time = 15000, resp_freq = 0.25)
#' ppgRespiration(simulation$pe)
ppgRespiration <- function(x,
                           peaks = NULL,
                           method = c("am", "bw", "fm"),
                           fs = 4,
                           resp_band = c(0.1, 0.5),
                           assay_name = NULL) {
  source <- .ppg_data(x, assay_name)
  method <- match.arg(method)
  fs <- .ppg_scalar(fs, "fs")
  if (!is.numeric(resp_band) || length(resp_band) != 2L ||
      any(!is.finite(resp_band)) || resp_band[1L] <= 0 ||
      resp_band[2L] <= resp_band[1L] ||
      resp_band[2L] >= fs / 2) {
    stop("resp_band must contain increasing positive frequencies below Nyquist.",
         call. = FALSE)
  }
  if (is.null(peaks)) {
    peaks <- ppgDetectPulses(x, assay_name = assay_name)
  }
  peaks <- .ppg_peaks(
    peaks,
    n_channels = ncol(source$data),
    n_time = nrow(source$data),
    require_time = TRUE
  )

  features <- list()
  if (method %in% c("am", "bw")) {
    wave <- pulseWaveFeatures(
      x, peaks = peaks, assay_name = assay_name
    )
    for (channel in seq_len(ncol(source$data))) {
      selected <- wave[wave$channel == channel, , drop = FALSE]
      if (nrow(selected) == 0L) {
        next
      }
      value <- if (method == "am") {
        selected$pulse_amplitude
      } else {
        source$data[
          cbind(selected$systolic_sample, rep(channel, nrow(selected)))
        ]
      }
      features[[length(features) + 1L]] <- data.frame(
        channel = as.integer(channel),
        time_sec = selected$systolic_time_sec,
        feature = value
      )
    }
  } else {
    intervals <- ecgRRintervals(x, peaks)
    for (channel in seq_len(ncol(source$data))) {
      selected <- intervals[intervals$channel == channel, , drop = FALSE]
      if (nrow(selected) == 0L) {
        next
      }
      features[[length(features) + 1L]] <- data.frame(
        channel = as.integer(channel),
        time_sec = selected$time_sec,
        feature = 60000 / selected$rr_ms
      )
    }
  }
  beats <- if (length(features)) {
    do.call(rbind, features)
  } else {
    data.frame(
      channel = integer(),
      time_sec = numeric(),
      feature = numeric()
    )
  }

  resampled <- rates <- list()
  for (channel in seq_len(ncol(source$data))) {
    selected <- beats[beats$channel == channel, , drop = FALSE]
    uniform <- .edr_resample(
      selected$time_sec, selected$feature, fs
    )
    if (is.null(uniform)) {
      frequency <- NA_real_
    } else {
      tolerance <- sqrt(.Machine$double.eps) *
        max(1, max(abs(uniform$value)))
      if (max(abs(uniform$value)) <= tolerance) {
        frequency <- NA_real_
      } else {
        frequency <- .edr_peak_freq(
          uniform$value, fs, resp_band
        )
      }
      resampled[[length(resampled) + 1L]] <- data.frame(
        channel = as.integer(channel),
        time_sec = uniform$time_sec,
        edr = uniform$value
      )
    }
    rates[[length(rates) + 1L]] <- data.frame(
      channel = as.integer(channel),
      resp_rate_hz = frequency,
      resp_rate_bpm = frequency * 60
    )
  }

  list(
    method = method,
    beats = beats,
    edr = if (length(resampled)) do.call(rbind, resampled) else NULL,
    resp_rate = do.call(rbind, rates)
  )
}


#' Simulate PPG with known fiducials
#'
#' Generates DC-offset PPG pulses from systolic and diastolic Gaussian waves
#' plus a sharp foot dip, with optional respiratory amplitude modulation.
#'
#' @param n_time Number of samples.
#' @param n_channels Number of PPG channels.
#' @param sr Sampling rate in Hz.
#' @param heart_rate Pulse rate in beats per minute.
#' @param dc Non-pulsatile baseline level.
#' @param sys_amp,dia_amp Systolic and diastolic wave amplitudes.
#' @param crest_ms Foot-to-systolic interval in milliseconds.
#' @param dia_offset_ms Foot-to-diastolic interval in milliseconds.
#' @param noise_sd Gaussian noise standard deviation.
#' @param foot_dip Depth of the local foot minimum.
#' @param resp_freq Optional respiratory amplitude-modulation frequency in Hz.
#' @param resp_depth Fractional respiratory modulation depth.
#'
#' @return A list with `pe`, a `PhysioExperiment`, and `truth`, a fiducial
#'   data frame (`beat`, `foot`, `systolic`, `diastolic`).
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' simulation <- make_ppg(n_time = 3000)
#' head(simulation$truth)
make_ppg <- function(n_time = 15000,
                     n_channels = 1,
                     sr = 125,
                     heart_rate = 72,
                     dc = 2,
                     sys_amp = 1,
                     dia_amp = 0.4,
                     crest_ms = 150,
                     dia_offset_ms = 350,
                     noise_sd = 0.02,
                     foot_dip = 0.1,
                     resp_freq = NULL,
                     resp_depth = 0.3) {
  integer_scalar <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 1L ||
        value != as.integer(value)) {
      stop(sprintf("%s must be a positive integer.", name),
           call. = FALSE)
    }
    as.integer(value)
  }
  n_time <- integer_scalar(n_time, "n_time")
  n_channels <- integer_scalar(n_channels, "n_channels")
  sr <- .ppg_scalar(sr, "sr")
  heart_rate <- .ppg_scalar(heart_rate, "heart_rate")
  if (!is.numeric(dc) || length(dc) != 1L || !is.finite(dc)) {
    stop("dc must be a finite numeric scalar.", call. = FALSE)
  }
  sys_amp <- .ppg_scalar(sys_amp, "sys_amp")
  dia_amp <- .ppg_scalar(dia_amp, "dia_amp", 0, inclusive = TRUE)
  crest_ms <- .ppg_scalar(crest_ms, "crest_ms")
  dia_offset_ms <- .ppg_scalar(dia_offset_ms, "dia_offset_ms")
  noise_sd <- .ppg_scalar(noise_sd, "noise_sd", 0, inclusive = TRUE)
  foot_dip <- .ppg_scalar(foot_dip, "foot_dip", 0, inclusive = TRUE)
  resp_depth <- .ppg_scalar(
    resp_depth, "resp_depth", 0, inclusive = TRUE
  )
  if (resp_depth >= 1) {
    stop("resp_depth must be less than 1.", call. = FALSE)
  }
  if (!is.null(resp_freq)) {
    resp_freq <- .ppg_scalar(resp_freq, "resp_freq")
    if (resp_freq >= sr / 2) {
      stop("resp_freq must be below Nyquist.", call. = FALSE)
    }
  }

  rr <- max(1L, as.integer(round(60 / heart_rate * sr)))
  crest_offset <- as.integer(round(crest_ms / 1000 * sr))
  diastolic_offset <- as.integer(round(dia_offset_ms / 1000 * sr))
  if (diastolic_offset <= crest_offset) {
    stop("dia_offset_ms must be greater than crest_ms.",
         call. = FALSE)
  }
  systolic_sigma <- max(1L, as.integer(round(0.030 * sr)))
  diastolic_sigma <- max(1L, as.integer(round(0.045 * sr)))
  foot_sigma <- max(1L, as.integer(round(0.010 * sr)))
  final_foot <- floor(
    n_time - diastolic_offset - 4 * diastolic_sigma
  )
  if (final_foot < rr) {
    stop("n_time is too short for one complete simulated pulse.",
         call. = FALSE)
  }
  feet <- seq.int(rr, final_foot, by = rr)
  systolic <- feet + crest_offset
  diastolic <- feet + diastolic_offset
  truth <- data.frame(
    beat = seq_along(feet),
    foot = as.integer(feet),
    systolic = as.integer(systolic),
    diastolic = as.integer(diastolic)
  )

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  for (channel in seq_len(n_channels)) {
    signal <- rep(dc, n_time) + stats::rnorm(n_time, sd = noise_sd)
    for (beat in seq_along(feet)) {
      modulation <- if (is.null(resp_freq)) {
        1
      } else {
        1 + resp_depth * sin(
          2 * pi * resp_freq * (feet[beat] - 1) / sr
        )
      }
      components <- list(
        c(center = systolic[beat],
          amplitude = sys_amp * modulation,
          sigma = systolic_sigma),
        c(center = diastolic[beat],
          amplitude = dia_amp,
          sigma = diastolic_sigma),
        c(center = feet[beat],
          amplitude = -foot_dip,
          sigma = foot_sigma)
      )
      for (component in components) {
        radius <- as.integer(ceiling(4 * component["sigma"]))
        index <- seq.int(
          max(1L, as.integer(component["center"]) - radius),
          min(n_time, as.integer(component["center"]) + radius)
        )
        signal[index] <- signal[index] +
          component["amplitude"] * exp(
            -((index - component["center"])^2) /
              (2 * component["sigma"]^2)
          )
      }
    }
    data[, channel] <- signal
  }
  colnames(data) <- paste0("PPG", seq_len(n_channels))

  experiment <- PhysioCore::PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = colnames(data),
      type = rep("PPG", n_channels)
    ),
    samplingRate = sr
  )
  list(pe = experiment, truth = truth)
}
