#' Detect Atrial Fibrillation / Frequent Ectopy from RR Irregularity
#'
#' Classifies the cardiac rhythm of each channel as regular sinus, atrial
#' fibrillation (AF), or frequent ectopy from the irregularity of the RR
#' interval series. HRV metrics are physiologically undefined in AF, so this
#' gate is used to invalidate them (see \code{rhythm_check} in
#' \code{\link{ecgHRVtime}}, \code{\link{ecgHRVfreq}} and
#' \code{\link{ecgHRVnonlinear}}).
#'
#' Three descriptors are combined: the ectopic-beat fraction from
#' \code{\link{ecgQualityCheck}}, the RMSSD/mean-RR normalized irregularity
#' index, and the normalized Shannon entropy of the successive-difference (dRR)
#' histogram. AF requires \emph{both} large beat-to-beat variability
#' (\code{rmssd_ratio}) and a broad, disordered dRR distribution
#' (\code{drr_entropy}); a large but concentrated irregularity (e.g. an ectopic
#' couplet) falls through to \code{"frequent_ectopy"}.
#'
#' @param rr A data.frame with columns \code{channel} and \code{rr_ms}, as
#'   returned by \code{\link{ecgRRintervals}}.
#' @param window_beats If non-NULL, an integer window length (in beats): instead
#'   of a single per-channel verdict, a continuous AF probability is returned for
#'   each non-overlapping window of \code{window_beats} beats (see Value).
#' @param ectopic_frac_max Maximum ectopic fraction for a \code{"sinus"}
#'   verdict (default 0.05).
#' @param rmssd_ratio_af Minimum RMSSD/mean-RR ratio consistent with AF
#'   (default 0.15).
#' @param entropy_af Minimum normalized dRR Shannon entropy consistent with AF
#'   (default 0.7, range 0-1).
#' @param n_bins Number of histogram bins for the dRR entropy (default 16).
#' @param drr_range Numeric length-2 dRR window (ms) for the entropy histogram
#'   (default c(-600, 600)).
#' @param threshold_ms Deviation threshold (ms) passed to
#'   \code{\link{ecgQualityCheck}} for the ectopic fraction (default 300).
#' @param min_beats Minimum RR intervals required to classify; fewer yields
#'   \code{NA} (default 20).
#' @return With \code{window_beats = NULL} (default), a data.frame with one row
#'   per channel and columns \code{channel}, \code{rhythm} (\code{"sinus"},
#'   \code{"AF"} or \code{"frequent_ectopy"}), \code{ectopic_frac},
#'   \code{rmssd_ratio}, \code{drr_entropy}, \code{af_prob} (a continuous AF
#'   probability in 0-1) and \code{n_beats}. With \code{window_beats} set, a
#'   data.frame with one row per non-overlapping window and columns
#'   \code{channel}, \code{window}, \code{beat_start}, \code{n_beats},
#'   \code{af_prob} and its components (\code{rmssd_ratio}, \code{drr_entropy},
#'   \code{sample_entropy}, \code{sd1_sd2}).
#' @references Tateno, K. & Glass, L. (2001). Automatic detection of atrial
#'   fibrillation using the coefficient of variation and density histograms of
#'   RR and dRR intervals. \emph{Medical & Biological Engineering & Computing},
#'   39(6), 664-671. Sarkar, S., Ritscher, D. & Mehra, R. (2008). A detector
#'   for a chronic implantable atrial tachyarrhythmia monitor. \emph{IEEE
#'   Transactions on Biomedical Engineering}, 55(3), 1219-1224.
#' @seealso \code{\link{ecgQualityCheck}}, \code{\link{ecgHRVtime}}
#' @export
#' @examples
#' set.seed(1)
#' rr <- data.frame(channel = 1L, rr_ms = 800 + rnorm(100, sd = 20),
#'                  time_sec = cumsum(rep(0.8, 100)))
#' ecgDetectAF(rr)
ecgDetectAF <- function(rr, window_beats = NULL, ectopic_frac_max = 0.05,
                        rmssd_ratio_af = 0.15, entropy_af = 0.7, n_bins = 16L,
                        drr_range = c(-600, 600), threshold_ms = 300,
                        min_beats = 20L) {
  stopifnot(is.data.frame(rr))
  if (!all(c("channel", "rr_ms") %in% names(rr))) {
    stop("Missing required columns: channel, rr_ms")
  }
  if (!is.null(window_beats)) {
    return(.ecg_af_windows(rr, window_beats, n_bins, drr_range))
  }
  # ecgQualityCheck() requires a time_sec column (its values are unused by the
  # ectopic test); synthesize one when absent so the documented channel + rr_ms
  # contract holds.
  if (!"time_sec" %in% names(rr)) {
    rr$time_sec <- stats::ave(rr$rr_ms, rr$channel,
                              FUN = function(v) cumsum(v) / 1000)
  }
  qc <- ecgQualityCheck(rr, threshold_ms = threshold_ms)

  channels <- unique(rr$channel)
  out <- list()
  for (ch in channels) {
    v <- rr$rr_ms[rr$channel == ch]
    ect <- qc$is_ectopic[qc$channel == ch]
    n_beats <- length(v)
    mean_rr <- mean(v)
    drr <- diff(v)
    rmssd <- if (length(drr) >= 1) sqrt(mean(drr^2)) else NA_real_
    rmssd_ratio <- if (isTRUE(mean_rr > 0)) rmssd / mean_rr else NA_real_
    ectopic_frac <- if (length(ect)) mean(ect) else NA_real_
    drr_entropy <- if (length(drr) >= 2) {
      .drr_shannon_entropy(drr, n_bins, drr_range)
    } else NA_real_

    rhythm <- if (n_beats < min_beats) {
      NA_character_
    } else if (isTRUE(rmssd_ratio >= rmssd_ratio_af) &&
               isTRUE(drr_entropy >= entropy_af)) {
      "AF"
    } else if (isTRUE(ectopic_frac > ectopic_frac_max)) {
      "frequent_ectopy"
    } else {
      "sinus"
    }

    af_prob <- if (n_beats >= 5L) .af_probability(v, n_bins, drr_range)$af_prob
      else NA_real_

    out[[length(out) + 1]] <- data.frame(
      channel = as.integer(ch), rhythm = rhythm,
      ectopic_frac = ectopic_frac, rmssd_ratio = rmssd_ratio,
      drr_entropy = drr_entropy, af_prob = af_prob, n_beats = n_beats,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# Normalized Shannon entropy (0..1) of the dRR histogram on fixed edges, so it
# is comparable across channels. Regular sinus -> ~0; AF -> ~1.
.drr_shannon_entropy <- function(drr, n_bins = 16L, drr_range = c(-600, 600)) {
  if (n_bins < 2L) stop("n_bins must be >= 2 (entropy normalized by log(n_bins))")
  edges <- seq(drr_range[1], drr_range[2], length.out = n_bins + 1L)
  x <- pmax(pmin(drr, drr_range[2]), drr_range[1])   # clip to range
  counts <- graphics::hist(x, breaks = edges, plot = FALSE)$counts
  total <- sum(counts)
  if (total == 0) return(NA_real_)
  p <- counts / total
  p <- p[p > 0]
  -sum(p * log(p)) / log(n_bins)
}

# Gate an already-computed HRV metrics data.frame on the rhythm verdict:
# blank numeric metric columns for non-sinus channels, append rhythm/hrv_valid,
# and emit at most one warning. Sinus metric VALUES are left unchanged (the two
# descriptive columns rhythm/hrv_valid are appended in every case). A channel
# whose rhythm is undetermined (too few beats, rhythm = NA) keeps its metrics and
# is flagged hrv_valid = NA (unknown), never silently marked invalid.
.apply_rhythm_gate <- function(metrics, rr) {
  if (is.null(metrics) || nrow(metrics) == 0) return(metrics)
  af <- ecgDetectAF(rr)
  rhythm <- af$rhythm[match(metrics$channel, af$channel)]
  invalid <- !is.na(rhythm) & rhythm != "sinus"

  if (any(invalid)) {
    for (col in setdiff(names(metrics), "channel")) {
      if (is.numeric(metrics[[col]])) metrics[[col]][invalid] <- NA_real_
    }
    warning(sprintf(
      paste0("HRV metrics undefined in atrial fibrillation / frequent ectopy ",
             "(%s); values set to NA. Set rhythm_check = FALSE to override."),
      paste(sprintf("ch %s: %s", metrics$channel[invalid], rhythm[invalid]),
            collapse = "; ")), call. = FALSE)
  }

  metrics$rhythm <- rhythm
  metrics$hrv_valid <- ifelse(is.na(rhythm), NA, rhythm == "sinus")
  metrics
}
