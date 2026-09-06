# Declare ggplot2 non-standard-evaluation column names.
utils::globalVariables(c(
  "time_sec", "rr_ms", "is_ectopic", "channel", "rr_n", "rr_n1", "x", "y",
  "freq", "power", "band", "xmin", "xmax", "ymin", "ymax", "sample", "value",
  "wave"))

# Resolve an rr argument that may be an rr data.frame or an ecgProcess object.
.ecg_rr_df <- function(rr) {
  if (inherits(rr, "PhysioExperiment")) {
    ecg <- S4Vectors::metadata(rr)$ecg
    if (is.null(ecg)) stop("PhysioExperiment has no ecg metadata; run ecgProcess()",
                           call. = FALSE)
    d <- if (!is.null(ecg$rr_corrected)) ecg$rr_corrected else ecg$rr
    return(d)
  }
  stopifnot(is.data.frame(rr), all(c("rr_ms", "time_sec") %in% names(rr)))
  if (is.null(rr$channel)) rr$channel <- 1L
  rr
}

#' Plot RR Tachogram
#'
#' Plots the RR interval series against time, marking ectopic beats when an
#' \code{is_ectopic} column is present.
#'
#' @param rr An rr data.frame (columns \code{rr_ms}, \code{time_sec} and
#'   optionally \code{channel}, \code{is_ectopic}) or an object returned by
#'   \code{\link{ecgProcess}}.
#' @return A ggplot object.
#' @seealso \code{\link{ecgRRintervals}}, \code{\link{ecgQualityCheck}}
#' @references Task Force (1996). Heart rate variability. \emph{Circulation},
#'   93(5), 1043-1065.
#' @export
#' @examples
#' pe <- make_ecg(20000, sr = 500)
#' rr <- ecgQualityCheck(ecgRRintervals(pe, ecgDetectRpeaks(pe)))
#' plotTachogram(rr)
plotTachogram <- function(rr) {
  d <- .ecg_rr_df(rr)
  d$channel <- factor(d$channel)
  p <- ggplot2::ggplot(d, ggplot2::aes(time_sec, rr_ms)) +
    ggplot2::geom_line(colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_point(colour = "#1b6ca8", size = 0.8)
  if (!is.null(d$is_ectopic) && any(d$is_ectopic, na.rm = TRUE)) {
    p <- p + ggplot2::geom_point(
      data = d[isTRUE(d$is_ectopic) | (!is.na(d$is_ectopic) & d$is_ectopic), ],
      colour = "#d1495b", size = 1.6, shape = 17)
  }
  if (nlevels(d$channel) > 1L) p <- p + ggplot2::facet_wrap(~ channel)
  p + ggplot2::labs(x = "Time (s)", y = "RR interval (ms)",
                    title = "RR tachogram") +
    PhysioCore::theme_physio()
}

.poincare_ellipse <- function(m, sd1, sd2, n = 120L) {
  th <- seq(0, 2 * pi, length.out = n)
  # SD2 axis along the identity line (1,1)/sqrt2; SD1 axis perpendicular
  u <- c(1, 1) / sqrt(2); v <- c(1, -1) / sqrt(2)
  x <- m + sd2 * cos(th) * u[1] + sd1 * sin(th) * v[1]
  y <- m + sd2 * cos(th) * u[2] + sd1 * sin(th) * v[2]
  data.frame(x = x, y = y)
}

#' Plot Poincare (Lorenz) Recurrence Map
#'
#' Scatters each RR interval against the next (\eqn{RR_n} vs \eqn{RR_{n+1}}) and
#' overlays the Poincare SD1/SD2 dispersion ellipse from
#' \code{\link{ecgHRVpoincare}}. The SD1 and SD2 used for the ellipse are
#' attached as the \code{"sd1"}/\code{"sd2"} attributes of the returned plot.
#'
#' @param rr An rr data.frame or an object from \code{\link{ecgProcess}}.
#' @return A ggplot object with \code{"sd1"}/\code{"sd2"} attributes (per
#'   channel).
#' @seealso \code{\link{ecgHRVpoincare}}, \code{\link{ecgHRVnonlinear}}
#' @references Brennan, M., Palaniswami, M. & Kamen, P. (2001). Do existing
#'   measures of Poincare plot geometry reflect nonlinear features of heart rate
#'   variability? \emph{IEEE TBME}, 48(11), 1342-1347.
#' @export
#' @examples
#' pe <- make_ecg(30000, sr = 500)
#' rr <- ecgRRintervals(pe, ecgDetectRpeaks(pe))
#' plotPoincare(rr)
plotPoincare <- function(rr) {
  d <- .ecg_rr_df(rr)
  pc <- ecgHRVpoincare(d)
  channels <- sort(unique(d$channel))

  pts <- do.call(rbind, lapply(channels, function(ch) {
    v <- d$rr_ms[d$channel == ch]
    if (length(v) < 2L) return(NULL)
    data.frame(channel = ch, rr_n = v[-length(v)], rr_n1 = v[-1])
  }))
  ell <- do.call(rbind, lapply(channels, function(ch) {
    row <- pc[pc$channel == ch, ]
    m <- mean(d$rr_ms[d$channel == ch])
    e <- .poincare_ellipse(m, row$sd1, row$sd2)
    e$channel <- ch; e
  }))
  pts$channel <- factor(pts$channel); ell$channel <- factor(ell$channel)

  p <- ggplot2::ggplot(pts, ggplot2::aes(rr_n, rr_n1)) +
    ggplot2::geom_point(colour = "#1b6ca8", alpha = 0.5, size = 0.8) +
    ggplot2::geom_path(data = ell, ggplot2::aes(x, y), inherit.aes = FALSE,
                       colour = "#d1495b", linewidth = 0.7) +
    ggplot2::coord_equal()
  if (nlevels(pts$channel) > 1L) p <- p + ggplot2::facet_wrap(~ channel)
  p <- p + ggplot2::labs(x = expression(RR[n] ~ "(ms)"),
                         y = expression(RR[n + 1] ~ "(ms)"),
                         title = "Poincare plot (SD1/SD2 ellipse)") +
    PhysioCore::theme_physio()
  attr(p, "sd1") <- stats::setNames(pc$sd1, pc$channel)
  attr(p, "sd2") <- stats::setNames(pc$sd2, pc$channel)
  p
}

#' Plot HRV Power Spectral Density with Frequency Bands
#'
#' Plots the RR-interval power spectrum with the VLF, LF and HF bands shaded and
#' their integrated powers annotated. The PSD and band integration reuse the same
#' internals as \code{\link{ecgHRVfreq}}, so the shaded band areas equal the
#' \code{ecgHRVfreq} band powers. The band powers are attached as the
#' \code{"band_power"} attribute.
#'
#' @param rr An rr data.frame or an object from \code{\link{ecgProcess}}.
#' @param method Spectral method ("welch" default, "ar" or "lomb"), matching
#'   \code{\link{ecgHRVfreq}}.
#' @param detrend Smoothness-priors detrending, matching \code{\link{ecgHRVfreq}}.
#' @param channel Channel to plot (default: first).
#' @param vlf_band,lf_band,hf_band Band edges (Hz).
#' @return A ggplot object with a \code{"band_power"} attribute.
#' @seealso \code{\link{ecgHRVfreq}}
#' @references Task Force (1996). Heart rate variability. \emph{Circulation},
#'   93(5), 1043-1065.
#' @export
#' @examples
#' pe <- make_ecg(40000, sr = 500)
#' rr <- ecgRRintervals(pe, ecgDetectRpeaks(pe))
#' plotHRVpsd(rr)
plotHRVpsd <- function(rr, method = c("welch", "ar", "lomb"), detrend = FALSE,
                       channel = NULL, vlf_band = c(0.003, 0.04),
                       lf_band = c(0.04, 0.15), hf_band = c(0.15, 0.4)) {
  method <- match.arg(method)
  d <- .ecg_rr_df(rr)
  if (is.null(channel)) channel <- sort(unique(d$channel))[1]
  dch <- d[d$channel == channel, ]
  ps <- switch(method,
    welch = .hrv_welch(dch$time_sec, dch$rr_ms, detrend = detrend),
    ar = .hrv_ar_burg(dch$time_sec, dch$rr_ms, detrend = detrend),
    lomb = .hrv_lomb(dch$time_sec, dch$rr_ms))

  df <- data.frame(freq = ps$freqs, power = ps$psd)
  df <- df[df$freq <= hf_band[2] * 1.2, , drop = FALSE]

  bp <- data.frame(
    band = c("VLF", "LF", "HF"),
    power = c(.band_power(ps$freqs, ps$psd, vlf_band),
              .band_power(ps$freqs, ps$psd, lf_band),
              .band_power(ps$freqs, ps$psd, hf_band)),
    stringsAsFactors = FALSE)
  bands <- data.frame(
    band = c("VLF", "LF", "HF"),
    xmin = c(vlf_band[1], lf_band[1], hf_band[1]),
    xmax = c(vlf_band[2], lf_band[2], hf_band[2]))

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bands, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf,
                   fill = band), alpha = 0.18) +
    ggplot2::geom_line(data = df, ggplot2::aes(freq, power),
                       colour = "grey20", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c(VLF = "#8d96a3", LF = "#66a182",
                                          HF = "#edae49"), name = NULL) +
    ggplot2::labs(x = "Frequency (Hz)", y = "PSD (ms^2/Hz)",
                  title = "HRV power spectrum",
                  subtitle = sprintf("VLF=%.0f  LF=%.0f  HF=%.0f (ms^2)",
                                     bp$power[1], bp$power[2], bp$power[3])) +
    PhysioCore::theme_physio()
  attr(p, "band_power") <- bp
  p
}

#' Plot Annotated ECG Waveform
#'
#' Plots the ECG waveform with P, QRS and T fiducials from
#' \code{\link{ecgDelineate}}. The detected R-peak sample indices are attached as
#' the \code{"peaks"} attribute (a beat-editing hook); the delineation is
#' attached as \code{"delineation"}.
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param peaks R-peak table from \code{\link{ecgDetectRpeaks}}.
#' @param channel Channel to plot (default: 1).
#' @param window_sec Seconds of signal to display from the start (default 10;
#'   NULL for the whole recording).
#' @param assay_name Input assay name (default: first assay).
#' @return A ggplot object with \code{"peaks"} and \code{"delineation"}
#'   attributes.
#' @seealso \code{\link{ecgDelineate}}, \code{\link{ecgDetectRpeaks}}
#' @references Pan, J. & Tompkins, W.J. (1985). A real-time QRS detection
#'   algorithm. \emph{IEEE TBME}, 32(3), 230-236.
#' @export
#' @examples
#' pe <- make_ecg(5000, sr = 500)
#' pk <- ecgDetectRpeaks(pe)
#' plotEcgWave(pe, pk)
plotEcgWave <- function(x, peaks, channel = 1, window_sec = 10,
                        assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  sr <- samplingRate(x)
  sig <- SummarizedExperiment::assay(x, assay_name)[, channel]
  n <- length(sig)
  last <- if (is.null(window_sec)) n else min(n, as.integer(window_sec * sr))
  idx <- seq_len(last)
  sig_df <- data.frame(time_sec = (idx - 1) / sr, value = sig[idx])

  del <- ecgDelineate(x, peaks[peaks$channel == channel, , drop = FALSE],
                      assay_name = assay_name)
  fid <- do.call(rbind, lapply(c("p_peak", "r_peak", "t_peak"), function(w) {
    s <- del[[w]]
    s <- s[is.finite(s) & s <= last]
    if (!length(s)) return(NULL)
    data.frame(time_sec = (s - 1) / sr, value = sig[s],
               wave = switch(w, p_peak = "P", r_peak = "R", t_peak = "T"))
  }))

  p <- ggplot2::ggplot(sig_df, ggplot2::aes(time_sec, value)) +
    ggplot2::geom_line(colour = "grey25", linewidth = 0.3)
  if (!is.null(fid)) {
    p <- p + ggplot2::geom_point(data = fid,
      ggplot2::aes(time_sec, value, colour = wave), size = 1.8) +
      ggplot2::scale_colour_manual(values = c(P = "#66a182", R = "#d1495b",
                                              `T` = "#edae49"), name = NULL)
  }
  p <- p + ggplot2::labs(x = "Time (s)", y = "Amplitude",
                         title = "Annotated ECG waveform") +
    PhysioCore::theme_physio()
  attr(p, "peaks") <- peaks$sample[peaks$channel == channel]
  attr(p, "delineation") <- del
  p
}
