# Symmetric (real-output) FFT band-pass via a folded-frequency mask.
.edr_bandpass <- function(sig, fs, flo, fhi) {
  sig <- sig - mean(sig)
  n <- length(sig)
  ft <- stats::fft(sig)
  freqs <- (0:(n - 1)) * fs / n
  f <- pmin(freqs, fs - freqs)
  Re(stats::fft(ft * as.numeric(f >= flo & f <= fhi), inverse = TRUE)) / n
}

# Frequency of the spectral peak within a band, with parabolic sub-bin
# interpolation for accuracy beyond the FFT resolution.
.edr_peak_freq <- function(sig, fs, band) {
  sig <- sig - mean(sig)
  n <- length(sig)
  if (n < 4L) return(NA_real_)
  psd <- Mod(stats::fft(sig))^2
  half <- seq_len(n %/% 2)
  f <- (half - 1) * fs / n
  p <- psd[half]
  inb <- which(f >= band[1] & f <= band[2])
  if (!length(inb)) return(NA_real_)
  pk <- inb[which.max(p[inb])]
  df <- if (length(f) > 1) f[2] - f[1] else fs / n
  if (pk > 1L && pk < length(p)) {
    y1 <- p[pk - 1L]; y2 <- p[pk]; y3 <- p[pk + 1L]
    denom <- y1 - 2 * y2 + y3
    delta <- if (denom != 0) 0.5 * (y1 - y3) / denom else 0
    f[pk] + delta * df
  } else {
    f[pk]
  }
}

# Uniform-grid resample of an irregular (time, value) series.
.edr_resample <- function(times, values, fs) {
  o <- order(times)
  times <- times[o]; values <- values[o]
  keep <- is.finite(times) & is.finite(values)
  times <- times[keep]; values <- values[keep]
  if (length(times) < 3L) return(NULL)
  tt <- seq(times[1], times[length(times)], by = 1 / fs)
  v <- stats::approx(times, values, xout = tt, rule = 2)$y
  v <- v - mean(v)                         # detrend: remove mean level
  list(time_sec = tt, value = v)
}

#' ECG-Derived Respiration (EDR)
#'
#' Estimates a respiration surrogate from beat-to-beat modulation of the QRS
#' complex (Moody et al. 1985). Respiration shifts the cardiac electrical axis,
#' modulating R-wave amplitude, QRS area or upslope; the per-beat feature series
#' is resampled to a uniform grid to form the EDR signal, and the respiratory
#' rate is read from its spectral peak.
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param peaks R-peak table from \code{\link{ecgDetectRpeaks}} (columns
#'   \code{channel}, \code{sample}, \code{time_sec}; \code{amplitude} is used by
#'   the "amplitude" method when present).
#' @param method Beat feature: "amplitude" (R-wave amplitude, default), "area"
#'   (QRS area) or "slope" (maximum QRS upslope).
#' @param fs Resampling rate (Hz) of the EDR signal (default 4).
#' @param resp_band Respiratory frequency band (Hz) for the rate estimate
#'   (default \code{c(0.1, 0.5)}, i.e. 6-30 breaths/min).
#' @param qrs_ms QRS half-window (ms) for the "area"/"slope" features
#'   (default 50).
#' @param assay_name Input assay name (default: first assay).
#' @return A list with:
#'   \describe{
#'     \item{method}{The feature used.}
#'     \item{beats}{A data.frame (\code{channel}, \code{time_sec},
#'       \code{feature}) of the per-beat feature.}
#'     \item{edr}{A data.frame (\code{channel}, \code{time_sec}, \code{edr}) of
#'       the uniform-resampled respiration surrogate.}
#'     \item{resp_rate}{A data.frame (\code{channel}, \code{resp_rate_hz},
#'       \code{resp_rate_bpm}) of the estimated respiratory rate.}
#'   }
#' @seealso \code{\link{ecgDetectRpeaks}}, \code{\link{ecgRSA}}
#' @references Moody, G.B. et al. (1985). Derivation of respiratory signals from
#'   multi-lead ECGs. \emph{Computers in Cardiology}, 12, 113-116.
#' @export
#' @examples
#' pe <- make_ecg(n_time = 10000, sr = 250)
#' pk <- ecgDetectRpeaks(pe)
#' edr <- ecgEDR(pe, pk, method = "area")
#' edr$resp_rate
ecgEDR <- function(x, peaks, method = c("amplitude", "area", "slope"),
                   fs = 4, resp_band = c(0.1, 0.5), qrs_ms = 50,
                   assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  half <- max(1L, as.integer(round(qrs_ms / 1000 * sr)))

  beats <- edr <- rate <- list()
  for (ch in sort(unique(peaks$channel))) {
    pk <- peaks[peaks$channel == ch, , drop = FALSE]
    pk <- pk[order(pk$sample), , drop = FALSE]
    sig <- data[, ch]

    feat <- vapply(seq_len(nrow(pk)), function(i) {
      p <- pk$sample[i]
      lo <- max(1L, p - half); hi <- min(n_time, p + half)
      seg <- sig[lo:hi]
      switch(method,
        amplitude = if (!is.null(pk$amplitude)) pk$amplitude[i] else sig[p],
        area = sum(abs(seg - stats::median(seg))),
        slope = max(abs(diff(seg))))
    }, numeric(1))

    beats[[length(beats) + 1]] <- data.frame(
      channel = as.integer(ch), time_sec = pk$time_sec, feature = feat,
      stringsAsFactors = FALSE)

    rs <- .edr_resample(pk$time_sec, feat, fs)
    if (is.null(rs)) {
      f_resp <- NA_real_
    } else {
      edr[[length(edr) + 1]] <- data.frame(
        channel = as.integer(ch), time_sec = rs$time_sec, edr = rs$value,
        stringsAsFactors = FALSE)
      f_resp <- .edr_peak_freq(rs$value, fs, resp_band)
    }
    rate[[length(rate) + 1]] <- data.frame(
      channel = as.integer(ch), resp_rate_hz = f_resp,
      resp_rate_bpm = f_resp * 60, stringsAsFactors = FALSE)
  }

  list(method = method,
       beats = do.call(rbind, beats),
       edr = if (length(edr)) do.call(rbind, edr) else NULL,
       resp_rate = do.call(rbind, rate))
}

#' Respiratory Sinus Arrhythmia (RSA)
#'
#' Quantifies respiratory sinus arrhythmia -- the respiration-linked oscillation
#' of the RR interval -- by the peak-valley method (Grossman) or the
#' Porges-Bohrer method (band-limited variance of the tachogram).
#'
#' @param rr A data.frame with columns \code{rr_ms} and \code{time_sec} (and
#'   optionally \code{channel}), or a numeric vector of RR intervals (ms), in
#'   which case beat times are taken as their cumulative sum.
#' @param resp Optional respiration surrogate (unused by the current methods,
#'   reserved for respiration-gated peak-valley).
#' @param method "peak_valley" (mean band-filtered peak-to-valley RR amplitude,
#'   default) or "porges_bohrer" (natural log of the band-filtered RR variance).
#' @param fs Resampling rate (Hz) of the tachogram (default 4).
#' @param resp_band Respiratory band (Hz) for filtering (default
#'   \code{c(0.12, 0.4)}).
#' @return A list with \code{method}, \code{rsa} (the RSA magnitude: mean
#'   peak-to-valley RR in ms, or \code{ln(variance)} for Porges-Bohrer), and,
#'   for "peak_valley", \code{n_cycles} (number of respiratory half-cycles used).
#' @seealso \code{\link{ecgEDR}}, \code{\link{ecgRRintervals}}
#' @references Grossman, P., van Beek, J. & Wientjes, C. (1990). A comparison of
#'   three quantification methods for estimation of respiratory sinus arrhythmia.
#'   \emph{Psychophysiology}, 27(6), 702-714. Lewis, G.F. et al. (2012).
#'   Statistical strategies to quantify RSA: are commonly used metrics
#'   equivalent? \emph{Biological Psychology}, 89(2), 349-364.
#' @export
#' @examples
#' set.seed(1)
#' t <- cumsum(rep(0.8, 300))
#' rr <- data.frame(rr_ms = 800 + 40 * sin(2 * pi * 0.25 * t) + rnorm(300, 5),
#'                  time_sec = t)
#' ecgRSA(rr)$rsa
ecgRSA <- function(rr, resp = NULL,
                   method = c("peak_valley", "porges_bohrer"),
                   fs = 4, resp_band = c(0.12, 0.4)) {
  method <- match.arg(method)
  if (is.data.frame(rr)) {
    stopifnot(all(c("rr_ms", "time_sec") %in% names(rr)))
    val <- rr$rr_ms; t <- rr$time_sec
  } else {
    val <- as.numeric(rr)
    t <- cumsum(val) / 1000
  }

  rs <- .edr_resample(t, val, fs)
  if (is.null(rs)) return(list(method = method, rsa = NA_real_))
  filt <- .edr_bandpass(rs$value, fs, resp_band[1], resp_band[2])

  if (method == "porges_bohrer") {
    v <- stats::var(filt)
    return(list(method = method, rsa = if (v > 0) log(v) else NA_real_))
  }

  # peak-valley: mean absolute difference between successive local extrema
  d <- diff(filt)
  s <- sign(d)
  s[s == 0] <- 1
  ext <- which(diff(s) != 0) + 1L            # local maxima/minima indices
  if (length(ext) < 2L) {
    return(list(method = method, rsa = NA_real_, n_cycles = 0L))
  }
  amps <- abs(diff(filt[ext]))
  list(method = method, rsa = mean(amps), n_cycles = length(amps))
}
