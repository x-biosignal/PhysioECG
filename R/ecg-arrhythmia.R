# Logistic squash of a descriptor to [0, 1].
.logit_squash <- function(x, center, scale) {
  if (!is.finite(x)) return(NA_real_)
  1 / (1 + exp(-(x - center) / scale))
}

# Sample entropy of a numeric series (Chebyshev template matching, m and
# tolerance r = r_factor * SD). O(n^2); used on RR windows of tens of beats.
.ecg_sampen_vec <- function(v, m = 2L, r_factor = 0.2) {
  n <- length(v)
  if (n < m + 2L) return(NA_real_)
  r <- r_factor * stats::sd(v)
  if (!is.finite(r) || r <= 0) return(0)
  matches <- function(mm) {
    N <- n - mm + 1L
    total <- 0
    for (i in seq_len(N - 1L)) {
      ti <- v[i + 0:(mm - 1L)]
      for (j in (i + 1L):N) {
        if (max(abs(ti - v[j + 0:(mm - 1L)])) <= r) total <- total + 1
      }
    }
    total
  }
  B <- matches(m)
  A <- matches(m + 1L)
  if (B > 0 && A > 0) -log(A / B) else NA_real_
}

# Continuous AF probability of an RR (ms) window from four irregularity
# descriptors: normalized RMSSD, dRR Shannon entropy, RR sample entropy and the
# Poincare SD1/SD2 ratio. Returns the probability and its components.
.af_probability <- function(v, n_bins = 16L, drr_range = c(-600, 600)) {
  na_out <- list(af_prob = NA_real_, rmssd_ratio = NA_real_,
                 drr_entropy = NA_real_, sample_entropy = NA_real_,
                 sd1_sd2 = NA_real_)
  v <- v[is.finite(v)]
  if (length(v) < 5L) return(na_out)
  drr <- diff(v)
  mean_rr <- mean(v)
  rmssd_ratio <- if (mean_rr > 0) sqrt(mean(drr^2)) / mean_rr else NA_real_
  ent <- .drr_shannon_entropy(drr, n_bins, drr_range)
  sampen <- .ecg_sampen_vec(v)
  sd_drr <- stats::sd(drr)
  sdnn <- stats::sd(v)
  sd1 <- sqrt(0.5) * sd_drr
  sd2 <- sqrt(max(2 * sdnn^2 - 0.5 * sd_drr^2, 0))
  sd1_sd2 <- if (is.finite(sd2) && sd2 > 0) sd1 / sd2 else NA_real_

  s <- c(.logit_squash(rmssd_ratio, 0.10, 0.03),
         ent,
         .logit_squash(sampen, 1.0, 0.4),
         .logit_squash(sd1_sd2, 0.6, 0.12))
  list(af_prob = mean(s, na.rm = TRUE), rmssd_ratio = rmssd_ratio,
       drr_entropy = ent, sample_entropy = sampen, sd1_sd2 = sd1_sd2)
}

# Per-window AF probability for each channel (non-overlapping windows).
.ecg_af_windows <- function(rr, window_beats, n_bins, drr_range) {
  window_beats <- as.integer(window_beats)
  stopifnot(window_beats >= 5L)
  out <- list()
  for (ch in unique(rr$channel)) {
    v <- rr$rr_ms[rr$channel == ch]
    starts <- seq(1L, length(v) - window_beats + 1L, by = window_beats)
    for (k in seq_along(starts)) {
      idx <- starts[k]:(starts[k] + window_beats - 1L)
      p <- .af_probability(v[idx], n_bins, drr_range)
      out[[length(out) + 1]] <- data.frame(
        channel = as.integer(ch), window = k, beat_start = starts[k],
        n_beats = window_beats, af_prob = p$af_prob,
        rmssd_ratio = p$rmssd_ratio, drr_entropy = p$drr_entropy,
        sample_entropy = p$sample_entropy, sd1_sd2 = p$sd1_sd2,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) {
    return(data.frame(channel = integer(0), window = integer(0),
                      beat_start = integer(0), n_beats = integer(0),
                      af_prob = numeric(0), rmssd_ratio = numeric(0),
                      drr_entropy = numeric(0), sample_entropy = numeric(0),
                      sd1_sd2 = numeric(0)))
  }
  do.call(rbind, out)
}

#' Detect Premature Ventricular Contractions (PVCs)
#'
#' Flags ventricular ectopic beats from three classic features: an abnormal QRS
#' morphology (low correlation with an average-beat template), prematurity of the
#' preceding RR interval, and a full compensatory pause (the premature beat plus
#' its following pause together span about two normal cycles). The compensatory
#' pause is what distinguishes a PVC from a non-compensatory supraventricular
#' ectopic beat or sinus arrhythmia.
#'
#' @param x A PhysioExperiment object with ECG data.
#' @param peaks R-peak table from \code{\link{ecgDetectRpeaks}} (columns
#'   \code{channel}, \code{sample}).
#' @param rr RR-interval table from \code{\link{ecgRRintervals}} (columns
#'   \code{channel}, \code{rr_ms}).
#' @param qrs_ms QRS window half-width in ms for template matching (default 100).
#' @param corr_threshold Template-correlation below which a beat is
#'   morphologically abnormal (default 0.9).
#' @param prematurity_max Prematurity ratio (RR_pre / local median RR) below
#'   which a beat is premature (default 0.85).
#' @param compensatory_min Minimum (RR_pre + RR_post) / (2 * local median RR)
#'   for a full compensatory pause (default 0.90).
#' @param require_compensatory If TRUE, also require a compensatory pause for a
#'   PVC verdict (default FALSE).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per beat: \code{channel}, \code{beat},
#'   \code{sample}, \code{rr_pre_ms}, \code{rr_post_ms}, \code{prematurity},
#'   \code{compensatory}, \code{morph_corr} and the logical \code{is_pvc}.
#' @seealso \code{\link{ecgDetectRpeaks}}, \code{\link{ecgRRintervals}},
#'   \code{\link{ecgDetectAF}}
#' @references Petrenas, A., Marozas, V. & Sornmo, L. (2015). Low-complexity
#'   detection of atrial fibrillation in continuous long-term monitoring.
#'   \emph{Computers in Biology and Medicine}, 65, 184-191.
#' @export
#' @examples
#' pe <- make_ecg_pqrst(n_time = 5000, sr = 500)$pe
#' pk <- ecgDetectRpeaks(pe)
#' rr <- ecgRRintervals(pe, pk)
#' head(ecgDetectPVC(pe, pk, rr))
ecgDetectPVC <- function(x, peaks, rr, qrs_ms = 100, corr_threshold = 0.9,
                         prematurity_max = 0.85, compensatory_min = 0.90,
                         require_compensatory = FALSE, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  half <- max(1L, as.integer(round(qrs_ms / 1000 * sr / 2)))

  out <- list()
  for (ch in sort(unique(peaks$channel))) {
    pk <- sort(peaks$sample[peaks$channel == ch])
    rr_ch <- rr$rr_ms[rr$channel == ch]        # length length(pk) - 1
    npk <- length(pk)
    if (npk < 3L) next
    med_rr <- stats::median(rr_ch, na.rm = TRUE)

    # QRS template from the median across beats (robust to a few PVCs).
    ok <- pk - half >= 1L & pk + half <= n_time
    win <- vapply(pk[ok], function(p) data[(p - half):(p + half), ch],
                  numeric(2L * half + 1L))
    template <- apply(win, 1, stats::median)

    for (i in seq_len(npk)) {
      p <- pk[i]
      morph <- if (p - half >= 1L && p + half <= n_time) {
        suppressWarnings(stats::cor(data[(p - half):(p + half), ch], template))
      } else NA_real_
      rr_pre <- if (i >= 2L) rr_ch[i - 1L] else NA_real_
      rr_post <- if (i <= npk - 1L) rr_ch[i] else NA_real_
      prem <- if (is.finite(rr_pre) && med_rr > 0) rr_pre / med_rr else NA_real_
      comp <- if (is.finite(rr_pre) && is.finite(rr_post) && med_rr > 0) {
        (rr_pre + rr_post) / (2 * med_rr)
      } else NA_real_

      is_pvc <- isTRUE(morph < corr_threshold) && isTRUE(prem < prematurity_max)
      if (require_compensatory) is_pvc <- is_pvc && isTRUE(comp >= compensatory_min)

      out[[length(out) + 1]] <- data.frame(
        channel = as.integer(ch), beat = i, sample = p,
        rr_pre_ms = rr_pre, rr_post_ms = rr_post, prematurity = prem,
        compensatory = comp, morph_corr = morph, is_pvc = is_pvc,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) {
    return(data.frame(channel = integer(0), beat = integer(0),
                      sample = integer(0), rr_pre_ms = numeric(0),
                      rr_post_ms = numeric(0), prematurity = numeric(0),
                      compensatory = numeric(0), morph_corr = numeric(0),
                      is_pvc = logical(0)))
  }
  do.call(rbind, out)
}
