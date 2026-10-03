# Respiratory inductance plethysmography and ventilatory analysis

#' Validate a respiratory-analysis scalar
#' @keywords internal
#' @noRd
.resp_scalar <- function(x, name, lower = 0, inclusive = FALSE) {
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


#' Validate respiratory sampling and frequency parameters
#' @keywords internal
#' @noRd
.resp_parameters <- function(sr, band, min_period_sec = NULL) {
  sr <- .resp_scalar(sr, "sr")
  if (!is.numeric(band) || length(band) != 2L ||
      any(!is.finite(band)) || band[1L] <= 0 ||
      band[2L] <= band[1L] || band[2L] >= sr / 2) {
    stop("band must contain increasing positive frequencies below Nyquist.",
         call. = FALSE)
  }
  out <- list(sr = sr, band = as.numeric(band))
  if (!is.null(min_period_sec)) {
    out$min_period_sec <- .resp_scalar(
      min_period_sec, "min_period_sec"
    )
  }
  out
}


#' Resolve a respiration signal source
#' @keywords internal
#' @noRd
.resp_signal <- function(x, channel = NULL, sr = NULL, assay_name = NULL) {
  if (inherits(x, "rip_calibration")) {
    if (!is.numeric(x$volume) || length(x$volume) < 1L ||
        any(!is.finite(x$volume))) {
      stop("x contains an invalid calibrated volume signal.",
           call. = FALSE)
    }
    data <- matrix(as.numeric(x$volume), ncol = 1L)
    source_sr <- x$sr
    channel_ids <- 1L
    labels <- "volume"
  } else if (inherits(x, "PhysioExperiment")) {
    if (is.null(assay_name)) {
      assay_name <- defaultAssay(x)
    }
    data <- SummarizedExperiment::assay(x, assay_name)
    if (!is.matrix(data) || !is.numeric(data) || nrow(data) < 1L ||
        ncol(data) < 1L) {
      stop("The selected assay must be a non-empty numeric matrix.",
           call. = FALSE)
    }
    source_sr <- samplingRate(x)
    channel_ids <- seq_len(ncol(data))
    labels <- colnames(data)
    col_data <- SummarizedExperiment::colData(x)
    if ("label" %in% colnames(col_data)) {
      labels <- as.character(col_data$label)
    }
    if (is.null(labels) || length(labels) != ncol(data)) {
      labels <- as.character(channel_ids)
    }
  } else if (is.numeric(x) && is.null(dim(x))) {
    if (length(x) < 1L) {
      stop("x must be a non-empty numeric vector.", call. = FALSE)
    }
    data <- matrix(as.numeric(x), ncol = 1L)
    source_sr <- sr
    channel_ids <- 1L
    labels <- "signal"
  } else {
    stop("x must be a PhysioExperiment, numeric vector, or rip_calibration.",
         call. = FALSE)
  }

  if (any(!is.finite(data))) {
    stop("Respiration signals must contain only finite values.",
         call. = FALSE)
  }
  source_sr <- .resp_scalar(source_sr, "sr")

  if (!is.null(channel)) {
    if (is.character(channel)) {
      if (anyNA(channel) || any(!nzchar(channel))) {
        stop("channel labels must be non-empty.", call. = FALSE)
      }
      selected <- match(channel, labels)
      if (anyNA(selected)) {
        stop(sprintf(
          "Unknown respiration channel label: %s.",
          paste(channel[is.na(selected)], collapse = ", ")
        ), call. = FALSE)
      }
    } else if (is.numeric(channel)) {
      if (any(!is.finite(channel)) ||
          any(channel != as.integer(channel)) ||
          any(channel < 1L | channel > ncol(data))) {
        stop("channel indices are outside the selected assay.",
             call. = FALSE)
      }
      selected <- as.integer(channel)
    } else {
      stop("channel must contain integer indices or channel labels.",
           call. = FALSE)
    }
    data <- data[, selected, drop = FALSE]
    channel_ids <- channel_ids[selected]
    labels <- labels[selected]
  }

  list(
    data = data,
    sr = source_sr,
    channel = as.integer(channel_ids),
    label = labels
  )
}


#' Empty respiratory onset table
#' @keywords internal
#' @noRd
.resp_empty_onsets <- function() {
  data.frame(
    channel = integer(),
    breath = integer(),
    onset_sample = integer(),
    peak_sample = integer(),
    offset_sample = integer(),
    onset_time = numeric(),
    peak_time = numeric(),
    offset_time = numeric(),
    insp_amplitude = numeric()
  )
}


#' Empty per-breath respiratory metrics
#' @keywords internal
#' @noRd
.resp_empty_breaths <- function() {
  data.frame(
    channel = integer(),
    breath = integer(),
    Ti = numeric(),
    Te = numeric(),
    Ttot = numeric(),
    duty_cycle = numeric(),
    tidal_volume = numeric(),
    mean_insp_flow = numeric(),
    rate_bpm = numeric()
  )
}


#' Mean inspiratory amplitude of a respiratory signal
#' @keywords internal
#' @noRd
.mean_insp_amp <- function(signal, sr, band, min_period_sec = 1) {
  onsets <- respiratoryOnsets(
    signal,
    band = band,
    min_period_sec = min_period_sec,
    sr = sr
  )
  amplitudes <- onsets$insp_amplitude
  amplitudes <- amplitudes[is.finite(amplitudes) & amplitudes > 0]
  if (length(amplitudes) == 0L) {
    return(NA_real_)
  }
  mean(amplitudes)
}


#' Calibrate two-belt respiratory inductance plethysmography
#'
#' Fits the Konno-Mead two-degree-of-freedom volume model from ribcage and
#' abdominal belt signals. Least-squares calibration uses a reference volume
#' trace; qualitative diagnostic calibration (QDC) equalizes the variability
#' of the two belt contributions.
#'
#' @param rc Ribcage signal as a numeric vector, or a `PhysioExperiment`.
#' @param ab Abdominal numeric vector, or a channel index/label when `rc` is a
#'   `PhysioExperiment`.
#' @param reference Optional reference volume signal, required for `"lsq"`.
#' @param sr Sampling rate in Hz for numeric signals. A `PhysioExperiment`
#'   supplies its own sampling rate.
#' @param method Calibration method, `"qdc"` or `"lsq"`.
#' @param band Two-element respiratory band in Hz.
#' @param window Optional inclusive calibration sample range `c(lo, hi)`.
#' @param fit_intercept Logical; include an intercept in least squares.
#' @param known_volume Optional known inspiratory volume used for absolute
#'   scaling.
#'
#' @return A `rip_calibration` object containing coefficients, calibration
#'   diagnostics, and the calibrated volume-proportional signal.
#'
#' @references
#' Konno K, Mead J (1967). Measurement of the separate volume changes of rib
#' cage and abdomen during breathing. *Journal of Applied Physiology*,
#' 22:407-422.
#'
#' Sackner MA, Watson H, Belsito AS, et al. (1989). Calibration of respiratory
#' inductive plethysmograph during natural breathing. *Journal of Applied
#' Physiology*, 66:410-420. \doi{10.1152/jappl.1989.66.1.410}
#'
#' @export
#'
#' @examples
#' sr <- 10
#' time <- (0:599) / sr
#' rc <- 2 * sin(2 * pi * 0.25 * time)
#' ab <- sin(2 * pi * 0.25 * time)
#' ripCalibrate(rc, ab, sr = sr)
ripCalibrate <- function(rc,
                         ab,
                         reference = NULL,
                         sr = NULL,
                         method = c("qdc", "lsq"),
                         band = c(0.1, 0.6),
                         window = NULL,
                         fit_intercept = FALSE,
                         known_volume = NULL) {
  method <- match.arg(method)
  if (!is.logical(fit_intercept) || length(fit_intercept) != 1L ||
      is.na(fit_intercept)) {
    stop("fit_intercept must be TRUE or FALSE.", call. = FALSE)
  }

  if (inherits(rc, "PhysioExperiment")) {
    if (missing(ab)) {
      stop("ab must identify the abdominal channel.", call. = FALSE)
    }
    source <- .resp_signal(rc, channel = ab)
    if (ncol(source$data) != 1L) {
      stop("ab must identify one abdominal channel.", call. = FALSE)
    }
    abdominal <- source$data[, 1L]

    rc_labels <- c("RC", "ribcage", "rib_cage")
    assay_name <- defaultAssay(rc)
    all_data <- SummarizedExperiment::assay(rc, assay_name)
    col_data <- SummarizedExperiment::colData(rc)
    labels <- colnames(all_data)
    if ("label" %in% colnames(col_data)) {
      labels <- as.character(col_data$label)
    }
    rc_index <- if (!is.null(labels)) {
      match(tolower(rc_labels), tolower(labels), nomatch = 0L)
    } else {
      integer()
    }
    rc_index <- rc_index[rc_index > 0L][1L]
    if (length(rc_index) == 0L || is.na(rc_index)) {
      remaining <- setdiff(seq_len(ncol(all_data)), source$channel)
      if (length(remaining) != 1L) {
        stop("Could not identify a unique ribcage channel.",
             call. = FALSE)
      }
      rc_index <- remaining
    }
    if (rc_index %in% source$channel) {
      stop("Ribcage and abdominal channels must be different.",
           call. = FALSE)
    }
    ribcage <- all_data[, rc_index]
    sr <- source$sr
  } else {
    if (missing(ab) || !is.numeric(rc) || !is.null(dim(rc)) ||
        !is.numeric(ab) || !is.null(dim(ab))) {
      stop("rc and ab must be numeric vectors of equal length.",
           call. = FALSE)
    }
    ribcage <- as.numeric(rc)
    abdominal <- as.numeric(ab)
  }

  if (length(ribcage) != length(abdominal) || length(ribcage) < 4L ||
      any(!is.finite(ribcage)) || any(!is.finite(abdominal))) {
    stop("rc and ab must be finite, equal-length vectors with at least 4 samples.",
         call. = FALSE)
  }
  parameters <- .resp_parameters(sr, band)
  sr <- parameters$sr
  band <- parameters$band
  n <- length(ribcage)

  if (is.null(window)) {
    window_index <- seq_len(n)
    window_echo <- c(1L, n)
  } else {
    if (!is.numeric(window) || length(window) != 2L ||
        any(!is.finite(window)) ||
        any(window != as.integer(window)) ||
        window[1L] < 1L || window[2L] > n ||
        window[1L] > window[2L]) {
      stop("window must be an inclusive sample range within the signals.",
           call. = FALSE)
    }
    window_echo <- as.integer(window)
    window_index <- seq.int(window_echo[1L], window_echo[2L])
  }
  minimum_window <- if (method == "lsq" && fit_intercept) 4L else 3L
  if (length(window_index) < minimum_window) {
    stop("The calibration window is too short for the selected model.",
         call. = FALSE)
  }

  rc_filtered <- .edr_bandpass(ribcage, sr, band[1L], band[2L])
  ab_filtered <- .edr_bandpass(abdominal, sr, band[1L], band[2L])
  intercept <- 0
  r_squared <- NA_real_

  if (method == "qdc") {
    sd_rc <- stats::sd(rc_filtered[window_index])
    sd_ab <- stats::sd(ab_filtered[window_index])
    tolerance <- sqrt(.Machine$double.eps) *
      max(1, max(abs(ab_filtered[window_index])))
    if (!is.finite(sd_ab) || sd_ab <= tolerance) {
      stop("QDC undefined: abdomen belt has zero variance.",
           call. = FALSE)
    }
    tolerance_rc <- sqrt(.Machine$double.eps) *
      max(1, max(abs(rc_filtered[window_index])))
    if (!is.finite(sd_rc) || sd_rc <= tolerance_rc) {
      stop("QDC undefined: ribcage belt has zero variance.",
           call. = FALSE)
    }
    coefficients <- c(rc = 1, ab = sd_rc / sd_ab)
  } else {
    if (is.null(reference) || !is.numeric(reference) ||
        !is.null(dim(reference)) || length(reference) != n ||
        any(!is.finite(reference))) {
      stop("reference must be a finite numeric vector matching rc for LSQ.",
           call. = FALSE)
    }
    reference_filtered <- .edr_bandpass(
      as.numeric(reference), sr, band[1L], band[2L]
    )
    design <- cbind(
      rc = rc_filtered[window_index],
      ab = ab_filtered[window_index]
    )
    if (fit_intercept) {
      design <- cbind(design, `(intercept)` = 1)
    }
    decomposition <- qr(design, tol = sqrt(.Machine$double.eps))
    if (decomposition$rank < ncol(design)) {
      stop("RIP belts are collinear; cannot solve LSQ.",
           call. = FALSE)
    }
    beta <- qr.coef(decomposition, reference_filtered[window_index])
    coefficients <- beta[c("rc", "ab")]
    if (fit_intercept) {
      intercept <- unname(beta["(intercept)"])
    }
    fitted <- as.numeric(design %*% beta)
    residual <- reference_filtered[window_index] - fitted
    total <- sum(
      (reference_filtered[window_index] -
         mean(reference_filtered[window_index]))^2
    )
    if (total <= .Machine$double.eps) {
      stop("reference has zero variance in the calibration window.",
           call. = FALSE)
    }
    r_squared <- 1 - sum(residual^2) / total
  }

  unscaled_volume <- coefficients["rc"] * rc_filtered +
    coefficients["ab"] * ab_filtered
  scale <- 1
  if (!is.null(known_volume)) {
    known_volume <- .resp_scalar(known_volume, "known_volume")
    amplitude <- .mean_insp_amp(unscaled_volume, sr, band)
    if (!is.finite(amplitude) || amplitude <= 0) {
      stop("known_volume scaling requires at least one detected breath.",
           call. = FALSE)
    }
    scale <- known_volume / amplitude
  }

  out <- list(
    method = method,
    coef = coefficients,
    k = unname(coefficients["ab"] / coefficients["rc"]),
    intercept = intercept,
    scale = scale,
    r_squared = r_squared,
    volume = as.numeric(scale * unscaled_volume + intercept),
    time_sec = (seq_len(n) - 1) / sr,
    sr = sr,
    band = band,
    window = window_echo
  )
  class(out) <- "rip_calibration"
  out
}


#' @export
print.rip_calibration <- function(x, ...) {
  cat(sprintf("<rip_calibration> method=%s, sr=%.4g Hz\n",
              x$method, x$sr))
  cat(sprintf(
    "  coefficients: RC=%.5g, AB=%.5g  K=%.5g  scale=%.5g\n",
    x$coef["rc"], x$coef["ab"], x$k, x$scale
  ))
  if (is.finite(x$r_squared)) {
    cat(sprintf("  R-squared: %.5f\n", x$r_squared))
  }
  invisible(x)
}


#' Detect respiratory cycle onsets
#'
#' Detects end-expiratory troughs in a zero-phase band-passed respiratory
#' signal. Consecutive troughs delimit complete breaths, and the intervening
#' maximum marks end inspiration.
#'
#' @param x A `PhysioExperiment`, numeric respiration vector, or
#'   `rip_calibration`.
#' @param channel Optional channel index or label.
#' @param band Two-element respiratory band in Hz.
#' @param min_period_sec Minimum breath period in seconds.
#' @param sr Sampling rate for a numeric vector.
#' @param assay_name Optional assay name for a `PhysioExperiment`.
#'
#' @return A data frame with one row per complete breath, containing sample and
#'   time coordinates for onset, peak, and offset plus inspiratory amplitude.
#'
#' @seealso [breathingRate()], [tidalMetrics()], [ripCalibrate()]
#'
#' @export
#'
#' @examples
#' sr <- 10
#' signal <- sin(2 * pi * 0.25 * (0:599) / sr)
#' head(respiratoryOnsets(signal, sr = sr))
respiratoryOnsets <- function(x,
                              channel = NULL,
                              band = c(0.1, 0.6),
                              min_period_sec = 1.0,
                              sr = NULL,
                              assay_name = NULL) {
  source <- .resp_signal(
    x, channel = channel, sr = sr, assay_name = assay_name
  )
  parameters <- .resp_parameters(
    source$sr, band, min_period_sec
  )
  sr <- parameters$sr
  band <- parameters$band
  min_period_sec <- parameters$min_period_sec
  if (nrow(source$data) < 4L) {
    return(.resp_empty_onsets())
  }

  rows <- list()
  min_distance <- max(1L, as.integer(round(min_period_sec * sr)))
  for (column in seq_len(ncol(source$data))) {
    filtered <- .edr_bandpass(
      source$data[, column], sr, band[1L], band[2L]
    )
    tolerance <- sqrt(.Machine$double.eps) *
      max(1, max(abs(filtered)))
    if (max(abs(filtered)) <= tolerance) {
      next
    }
    troughs <- .find_local_maxima(-filtered, min_distance)
    if (length(troughs) < 2L) {
      next
    }

    breath_number <- 0L
    for (i in seq_len(length(troughs) - 1L)) {
      onset <- troughs[i]
      offset <- troughs[i + 1L]
      if (offset - onset < 2L) {
        next
      }
      interior <- seq.int(onset + 1L, offset - 1L)
      peak <- interior[which.max(filtered[interior])]
      amplitude <- filtered[peak] - filtered[onset]
      if (!is.finite(amplitude) || amplitude <= tolerance) {
        next
      }
      breath_number <- breath_number + 1L
      rows[[length(rows) + 1L]] <- data.frame(
        channel = source$channel[column],
        breath = breath_number,
        onset_sample = as.integer(onset),
        peak_sample = as.integer(peak),
        offset_sample = as.integer(offset),
        onset_time = (onset - 1) / sr,
        peak_time = (peak - 1) / sr,
        offset_time = (offset - 1) / sr,
        insp_amplitude = amplitude
      )
    }
  }

  if (length(rows) == 0L) {
    .resp_empty_onsets()
  } else {
    do.call(rbind, rows)
  }
}


#' Estimate breathing rate
#'
#' Estimates rate from complete inter-trough breath periods or from the
#' dominant respiratory-band spectral peak.
#'
#' @param x,channel,band,min_period_sec,sr,assay_name As in
#'   [respiratoryOnsets()].
#' @param method `"peaks"` for breath intervals or `"spectral"` for the
#'   dominant respiratory frequency.
#'
#' @return A data frame with one row per selected channel and columns
#'   `channel`, `method`, `rate_bpm`, `n_breaths`, and `sd_bpm`.
#'
#' @export
#'
#' @examples
#' sr <- 10
#' signal <- sin(2 * pi * 0.25 * (0:599) / sr)
#' breathingRate(signal, sr = sr)
breathingRate <- function(x,
                          channel = NULL,
                          method = c("peaks", "spectral"),
                          band = c(0.1, 0.6),
                          min_period_sec = 1.0,
                          sr = NULL,
                          assay_name = NULL) {
  method <- match.arg(method)
  source <- .resp_signal(
    x, channel = channel, sr = sr, assay_name = assay_name
  )
  parameters <- .resp_parameters(
    source$sr, band, min_period_sec
  )
  sr <- parameters$sr
  band <- parameters$band

  if (method == "peaks") {
    onsets <- respiratoryOnsets(
      x,
      channel = channel,
      band = band,
      min_period_sec = parameters$min_period_sec,
      sr = sr,
      assay_name = assay_name
    )
    output <- lapply(source$channel, function(channel_id) {
      selected <- onsets[onsets$channel == channel_id, , drop = FALSE]
      rates <- 60 / (selected$offset_time - selected$onset_time)
      rates <- rates[is.finite(rates) & rates > 0]
      data.frame(
        channel = channel_id,
        method = method,
        rate_bpm = if (length(rates)) mean(rates) else NA_real_,
        n_breaths = as.integer(length(rates)),
        sd_bpm = if (length(rates) > 1L) stats::sd(rates) else NA_real_
      )
    })
  } else {
    output <- lapply(seq_len(ncol(source$data)), function(column) {
      filtered <- .edr_bandpass(
        source$data[, column], sr, band[1L], band[2L]
      )
      tolerance <- sqrt(.Machine$double.eps) *
        max(1, max(abs(filtered)))
      frequency <- if (length(filtered) < 4L ||
                       max(abs(filtered)) <= tolerance) {
        NA_real_
      } else {
        .edr_peak_freq(filtered, sr, band)
      }
      data.frame(
        channel = source$channel[column],
        method = method,
        rate_bpm = frequency * 60,
        n_breaths = NA_integer_,
        sd_bpm = NA_real_
      )
    })
  }
  do.call(rbind, output)
}


#' Calculate tidal breathing metrics
#'
#' Derives inspiratory and expiratory timing, duty cycle, tidal amplitude,
#' mean inspiratory flow, breathing rate, and minute ventilation.
#'
#' @param x A respiration signal source accepted by [respiratoryOnsets()].
#' @param onsets Optional precomputed table from [respiratoryOnsets()].
#' @param channel,band,min_period_sec,sr,assay_name As in
#'   [respiratoryOnsets()].
#' @param calibration Optional `rip_calibration` or positive numeric
#'   amplitude-to-liters scale.
#'
#' @return A named list containing per-breath `breaths` and per-channel
#'   `summary` data frames.
#'
#' @export
#'
#' @examples
#' sr <- 10
#' signal <- sin(2 * pi * 0.25 * (0:599) / sr)
#' tidalMetrics(signal, sr = sr)
tidalMetrics <- function(x,
                         onsets = NULL,
                         channel = NULL,
                         calibration = NULL,
                         band = c(0.1, 0.6),
                         min_period_sec = 1.0,
                         sr = NULL,
                         assay_name = NULL) {
  source <- .resp_signal(
    x, channel = channel, sr = sr, assay_name = assay_name
  )
  parameters <- .resp_parameters(
    source$sr, band, min_period_sec
  )
  if (is.null(onsets)) {
    onsets <- respiratoryOnsets(
      x,
      channel = channel,
      band = parameters$band,
      min_period_sec = parameters$min_period_sec,
      sr = parameters$sr,
      assay_name = assay_name
    )
  }
  required <- c(
    "channel", "breath", "onset_time", "peak_time", "offset_time",
    "insp_amplitude"
  )
  if (!is.data.frame(onsets) || !all(required %in% names(onsets))) {
    stop("onsets must be a respiratoryOnsets() data frame.",
         call. = FALSE)
  }
  if (nrow(onsets) > 0L &&
      (any(!is.finite(as.matrix(onsets[, required[-1L], drop = FALSE]))) ||
       any(!onsets$channel %in% source$channel))) {
    stop("onsets contain invalid values or channels.", call. = FALSE)
  }

  scale <- 1
  if (inherits(calibration, "rip_calibration")) {
    scale <- calibration$scale
  } else if (!is.null(calibration)) {
    scale <- .resp_scalar(calibration, "calibration")
  }
  if (!is.numeric(scale) || length(scale) != 1L ||
      !is.finite(scale) || scale <= 0) {
    stop("calibration must provide a positive finite scale.",
         call. = FALSE)
  }

  if (nrow(onsets) == 0L) {
    breaths <- .resp_empty_breaths()
  } else {
    Ti <- onsets$peak_time - onsets$onset_time
    Te <- onsets$offset_time - onsets$peak_time
    Ttot <- onsets$offset_time - onsets$onset_time
    valid <- Ti > 0 & Te > 0 & Ttot > 0
    if (!all(valid)) {
      warning("Invalid breath timing rows were omitted.", call. = FALSE)
      onsets <- onsets[valid, , drop = FALSE]
      Ti <- Ti[valid]
      Te <- Te[valid]
      Ttot <- Ttot[valid]
    }
    tidal_volume <- scale * onsets$insp_amplitude
    breaths <- data.frame(
      channel = as.integer(onsets$channel),
      breath = as.integer(onsets$breath),
      Ti = Ti,
      Te = Te,
      Ttot = Ttot,
      duty_cycle = Ti / Ttot,
      tidal_volume = tidal_volume,
      mean_insp_flow = tidal_volume / Ti,
      rate_bpm = 60 / Ttot
    )
  }

  summary_rows <- lapply(source$channel, function(channel_id) {
    selected <- breaths[breaths$channel == channel_id, , drop = FALSE]
    if (nrow(selected) == 0L) {
      data.frame(
        channel = channel_id,
        n_breaths = 0L,
        rate_bpm = NA_real_,
        tidal_volume_mean = NA_real_,
        minute_ventilation = NA_real_,
        duty_cycle_mean = NA_real_
      )
    } else {
      mean_rate <- mean(selected$rate_bpm)
      mean_volume <- mean(selected$tidal_volume)
      data.frame(
        channel = channel_id,
        n_breaths = as.integer(nrow(selected)),
        rate_bpm = mean_rate,
        tidal_volume_mean = mean_volume,
        minute_ventilation = mean_volume * mean_rate,
        duty_cycle_mean = mean(selected$duty_cycle)
      )
    }
  })

  list(breaths = breaths, summary = do.call(rbind, summary_rows))
}


#' Fit a two-segment linear breakpoint
#' @keywords internal
#' @noRd
.two_segment_break <- function(x, y, min_segment = 3L) {
  n <- length(x)
  if (n != length(y) || n < 2L * min_segment) {
    return(NULL)
  }
  best <- NULL
  best_sse <- Inf
  for (split in seq.int(min_segment, n - min_segment)) {
    left <- seq_len(split)
    right <- seq.int(split + 1L, n)
    design_left <- cbind(1, x[left])
    design_right <- cbind(1, x[right])
    qr_left <- qr(design_left, tol = sqrt(.Machine$double.eps))
    qr_right <- qr(design_right, tol = sqrt(.Machine$double.eps))
    if (qr_left$rank < 2L || qr_right$rank < 2L) {
      next
    }
    coef_left <- qr.coef(qr_left, y[left])
    coef_right <- qr.coef(qr_right, y[right])
    residual_left <- y[left] - as.numeric(design_left %*% coef_left)
    residual_right <- y[right] - as.numeric(design_right %*% coef_right)
    sse <- sum(residual_left^2) + sum(residual_right^2)
    if (is.finite(sse) && sse < best_sse) {
      slope_difference <- coef_right[2L] - coef_left[2L]
      x_break <- if (abs(slope_difference) >
                     sqrt(.Machine$double.eps)) {
        (coef_left[1L] - coef_right[1L]) / slope_difference
      } else {
        NA_real_
      }
      best_sse <- sse
      best <- list(
        index = as.integer(split),
        x_break = unname(x_break),
        slope1 = unname(coef_left[2L]),
        slope2 = unname(coef_right[2L]),
        sse = sse
      )
    }
  }
  best
}


#' Construct an unavailable ventilatory threshold
#' @keywords internal
#' @noRd
.empty_ventilatory_threshold <- function() {
  list(
    vo2 = NA_real_,
    vco2 = NA_real_,
    ve = NA_real_,
    index = NA_integer_,
    time = NA_real_,
    slope1 = NA_real_,
    slope2 = NA_real_,
    sse = NA_real_
  )
}


#' Test for a meaningful slope increase
#' @keywords internal
#' @noRd
.resp_slope_increase <- function(breakpoint) {
  if (is.null(breakpoint) ||
      any(!is.finite(c(
        breakpoint$x_break, breakpoint$slope1, breakpoint$slope2
      )))) {
    return(FALSE)
  }
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, abs(breakpoint$slope1), abs(breakpoint$slope2))
  breakpoint$slope2 > breakpoint$slope1 + tolerance
}


#' Detect ventilatory thresholds from incremental exercise data
#'
#' VT1 is estimated by V-slope breakpoint regression of VCO2 on VO2. VT2 is
#' estimated above VT1 from the breakpoint of VE on VCO2.
#'
#' @param vo2 Numeric VO2 vector, or a data frame containing `vo2`, `vco2`,
#'   and `ve`.
#' @param vco2,ve Numeric vectors matching `vo2` when it is not a data frame.
#' @param time Optional time vector. For data-frame input, an existing `time`
#'   column is used unless this argument is supplied.
#' @param average_n Positive integer number of ordered samples per block mean.
#' @param min_segment Minimum number of observations fitted on each side of a
#'   breakpoint.
#'
#' @return A `ventilatory_threshold` object with VT1, VT2, and the data actually
#'   fitted.
#'
#' @references
#' Beaver WL, Wasserman K, Whipp BJ (1986). A new method for detecting
#' anaerobic threshold by gas exchange. *Journal of Applied Physiology*,
#' 60:2020-2027. \doi{10.1152/jappl.1986.60.6.2020}
#'
#' @export
#'
#' @examples
#' vo2 <- seq(1, 3, by = 0.05)
#' vco2 <- ifelse(vo2 <= 2, 0.9 * vo2, 1.8 + 1.3 * (vo2 - 2))
#' ventilatoryThreshold(vo2, vco2, ve = 25 * vco2)
ventilatoryThreshold <- function(vo2,
                                 vco2 = NULL,
                                 ve = NULL,
                                 time = NULL,
                                 average_n = 1L,
                                 min_segment = 3L) {
  if (!is.numeric(average_n) || length(average_n) != 1L ||
      !is.finite(average_n) || average_n < 1L ||
      average_n != as.integer(average_n)) {
    stop("average_n must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(min_segment) || length(min_segment) != 1L ||
      !is.finite(min_segment) || min_segment < 2L ||
      min_segment != as.integer(min_segment)) {
    stop("min_segment must be an integer of at least 2.", call. = FALSE)
  }
  average_n <- as.integer(average_n)
  min_segment <- as.integer(min_segment)

  if (is.data.frame(vo2)) {
    data <- vo2
    required <- c("vo2", "vco2", "ve")
    if (!all(required %in% names(data))) {
      stop("vo2 data frame must contain vo2, vco2, and ve columns.",
           call. = FALSE)
    }
    if (is.null(time) && "time" %in% names(data)) {
      time <- data$time
    }
    vo2 <- data$vo2
    vco2 <- data$vco2
    ve <- data$ve
  }
  if (!is.numeric(vo2) || !is.null(dim(vo2)) ||
      !is.numeric(vco2) || !is.null(dim(vco2)) ||
      !is.numeric(ve) || !is.null(dim(ve)) ||
      length(vo2) != length(vco2) || length(vo2) != length(ve) ||
      length(vo2) < 1L ||
      any(!is.finite(vo2)) || any(!is.finite(vco2)) ||
      any(!is.finite(ve))) {
    stop("vo2, vco2, and ve must be finite equal-length numeric vectors.",
         call. = FALSE)
  }
  if (!is.null(time) &&
      (!is.numeric(time) || !is.null(dim(time)) ||
       length(time) != length(vo2) || any(!is.finite(time)))) {
    stop("time must be a finite numeric vector matching vo2.",
         call. = FALSE)
  }

  order_index <- order(vo2, seq_along(vo2))
  fit_data <- data.frame(
    vo2 = as.numeric(vo2[order_index]),
    vco2 = as.numeric(vco2[order_index]),
    ve = as.numeric(ve[order_index])
  )
  if (!is.null(time)) {
    fit_data$time <- as.numeric(time[order_index])
  }
  if (average_n > 1L) {
    block <- ceiling(seq_len(nrow(fit_data)) / average_n)
    fit_data <- as.data.frame(lapply(
      fit_data,
      function(values) as.numeric(tapply(values, block, mean))
    ))
  }

  empty <- .empty_ventilatory_threshold()
  vt1 <- empty
  vt2 <- empty
  first <- .two_segment_break(
    fit_data$vo2, fit_data$vco2, min_segment
  )
  first_valid <- .resp_slope_increase(first) &&
    first$x_break >= min(fit_data$vo2) &&
    first$x_break <= max(fit_data$vo2)

  if (first_valid) {
    vt1_index <- which.min(abs(fit_data$vo2 - first$x_break))
    vt1 <- list(
      vo2 = first$x_break,
      vco2 = stats::approx(
        fit_data$vo2, fit_data$vco2,
        xout = first$x_break, rule = 2, ties = mean
      )$y,
      ve = stats::approx(
        fit_data$vo2, fit_data$ve,
        xout = first$x_break, rule = 2, ties = mean
      )$y,
      index = as.integer(vt1_index),
      time = if ("time" %in% names(fit_data)) {
        fit_data$time[vt1_index]
      } else {
        NA_real_
      },
      slope1 = first$slope1,
      slope2 = first$slope2,
      sse = first$sse
    )

    above_index <- which(fit_data$vo2 >= first$x_break)
    above <- fit_data[above_index, , drop = FALSE]
    second <- .two_segment_break(above$vco2, above$ve, min_segment)
    monotone_vco2 <- all(diff(above$vco2) > 0)
    second_valid <- monotone_vco2 &&
      .resp_slope_increase(second) &&
      second$x_break >= min(above$vco2) &&
      second$x_break <= max(above$vco2)
    if (second_valid) {
      vt2_vo2 <- stats::approx(
        above$vco2, above$vo2,
        xout = second$x_break, rule = 2, ties = mean
      )$y
      vt2_index <- which.min(abs(fit_data$vo2 - vt2_vo2))
      vt2 <- list(
        vo2 = vt2_vo2,
        vco2 = second$x_break,
        ve = stats::approx(
          above$vco2, above$ve,
          xout = second$x_break, rule = 2, ties = mean
        )$y,
        index = as.integer(vt2_index),
        time = if ("time" %in% names(fit_data)) {
          fit_data$time[vt2_index]
        } else {
          NA_real_
        },
        slope1 = second$slope1,
        slope2 = second$slope2,
        sse = second$sse
      )
    }
  }

  out <- list(
    method = "v-slope",
    vt1 = vt1,
    vt2 = vt2,
    data = fit_data
  )
  class(out) <- "ventilatory_threshold"
  out
}


#' @export
print.ventilatory_threshold <- function(x, ...) {
  cat("<ventilatory_threshold> method=v-slope\n")
  if (is.finite(x$vt1$vo2)) {
    cat(sprintf(
      "  VT1: VO2 %.4g, VCO2 %.4g, VE %.4g (slopes %.4g -> %.4g)\n",
      x$vt1$vo2, x$vt1$vco2, x$vt1$ve,
      x$vt1$slope1, x$vt1$slope2
    ))
  } else {
    cat("  VT1: not detected\n")
  }
  if (is.finite(x$vt2$vo2)) {
    cat(sprintf(
      "  VT2: VO2 %.4g, VCO2 %.4g, VE %.4g (slopes %.4g -> %.4g)\n",
      x$vt2$vo2, x$vt2$vco2, x$vt2$ve,
      x$vt2$slope1, x$vt2$slope2
    ))
  } else {
    cat("  VT2: not detected\n")
  }
  invisible(x)
}
