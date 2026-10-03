#' Heart-Rate Recovery (HRR)
#'
#' Computes post-exercise heart-rate recovery: the fall in heart rate from the
#' exercise peak to fixed offsets into recovery (Cole et al. 1999; Imai et al.
#' 1994), plus the recovery slope. Accepts either a heart-rate (bpm) or an RR
#' (ms) series.
#'
#' @param rr_or_hr A numeric vector of heart rate (bpm) or RR intervals (ms), or
#'   a data.frame with a \code{time_sec} column and an \code{hr} or \code{rr_ms}
#'   column.
#' @param time_sec Time of each sample in seconds (required for a vector input;
#'   taken from the \code{time_sec} column of a data.frame otherwise).
#' @param peak_time Time (seconds) of the exercise peak / start of recovery. If
#'   \code{NULL} (default), auto-detected as the time of maximum heart rate.
#' @param recovery_offsets_sec Offsets (seconds) after the peak at which to
#'   measure recovery (default \code{c(60, 120)}, i.e. HRR1 and HRR2).
#' @param input Interpretation of a vector input: "auto" (default; decided by
#'   magnitude -- values above 250 are treated as RR in ms), "hr" or "rr".
#' @param abnormal_bpm HRR at 1 minute below which recovery is flagged abnormal
#'   (default 12 bpm, per Cole et al. 1999).
#' @return A list with:
#'   \describe{
#'     \item{peak_time, peak_hr}{Recovery start time (s) and peak HR (bpm).}
#'     \item{hrr}{A data.frame with columns \code{offset_sec},
#'       \code{hr_at_offset} and \code{hrr_bpm} (peak HR minus HR at the offset).}
#'     \item{hrr1, hrr2}{Convenience scalars for the first two offsets.}
#'     \item{slope_bpm_per_min}{Slope of a linear fit of HR vs time over the
#'       recovery window (negative during recovery).}
#'     \item{hrr_1min}{HRR at exactly 60 s (used for the abnormal flag).}
#'     \item{abnormal}{TRUE if \code{hrr_1min < abnormal_bpm}.}
#'   }
#' @seealso \code{\link{ecgRRintervals}}, \code{\link{ecgHRVtime}}
#' @references Cole, C.R. et al. (1999). Heart-rate recovery immediately after
#'   exercise as a predictor of mortality. \emph{New England Journal of
#'   Medicine}, 341(18), 1351-1357. Imai, K. et al. (1994). Vagally mediated
#'   heart rate recovery after exercise is accelerated in athletes but blunted
#'   in patients with chronic heart failure. \emph{Journal of the American
#'   College of Cardiology}, 24(6), 1529-1535.
#' @export
#' @examples
#' # HR decays linearly from 170 to 120 bpm over 120 s of recovery
#' t <- seq(0, 120, by = 1)
#' hr <- 170 - (170 - 120) * t / 120
#' res <- ecgHRR(hr, time_sec = t, peak_time = 0)
#' res$hrr
#' res$abnormal
ecgHRR <- function(rr_or_hr, time_sec = NULL, peak_time = NULL,
                   recovery_offsets_sec = c(60, 120),
                   input = c("auto", "hr", "rr"), abnormal_bpm = 12) {
  input <- match.arg(input)

  if (is.data.frame(rr_or_hr)) {
    df <- rr_or_hr
    if (is.null(time_sec)) {
      if (!"time_sec" %in% names(df)) {
        stop("data.frame input requires a 'time_sec' column", call. = FALSE)
      }
      time_sec <- df$time_sec
    }
    if ("hr" %in% names(df)) {
      val <- df$hr
      if (input == "auto") input <- "hr"
    } else if ("rr_ms" %in% names(df)) {
      val <- df$rr_ms
      if (input == "auto") input <- "rr"
    } else {
      stop("data.frame input requires an 'hr' or 'rr_ms' column", call. = FALSE)
    }
  } else {
    val <- as.numeric(rr_or_hr)
    if (is.null(time_sec)) {
      stop("time_sec is required for a vector input", call. = FALSE)
    }
  }
  stopifnot(length(val) == length(time_sec))

  if (input == "auto") {
    input <- if (stats::median(val, na.rm = TRUE) > 250) "rr" else "hr"
  }
  hr <- if (input == "rr") 60000 / val else val          # bpm
  t <- as.numeric(time_sec)

  ok <- is.finite(hr) & is.finite(t)
  hr <- hr[ok]; t <- t[ok]
  if (length(hr) < 2L) stop("need at least two finite samples", call. = FALSE)
  o <- order(t); t <- t[o]; hr <- hr[o]

  if (is.null(peak_time)) {
    pk <- which.max(hr)
    peak_time <- t[pk]
    peak_hr <- hr[pk]
  } else {
    peak_hr <- stats::approx(t, hr, xout = peak_time, rule = 2)$y
  }

  hr_at <- function(offset) {
    stats::approx(t, hr, xout = peak_time + offset, rule = 1)$y
  }
  rows <- lapply(recovery_offsets_sec, function(off) {
    v <- hr_at(off)
    data.frame(offset_sec = off, hr_at_offset = v, hrr_bpm = peak_hr - v,
               stringsAsFactors = FALSE)
  })
  hrr <- do.call(rbind, rows)

  # recovery slope (bpm per minute) over [peak, peak + last offset]
  wmax <- peak_time + max(recovery_offsets_sec)
  sel <- t >= peak_time & t <= wmax
  slope <- if (sum(sel) >= 3L && stats::sd(t[sel]) > 0) {
    unname(stats::coef(stats::lm(hr[sel] ~ t[sel]))[2]) * 60
  } else {
    NA_real_
  }

  hr_1min <- hr_at(60)
  hrr_1min <- peak_hr - hr_1min
  abnormal <- if (is.finite(hrr_1min)) hrr_1min < abnormal_bpm else NA

  list(peak_time = peak_time, peak_hr = peak_hr, hrr = hrr,
       hrr1 = hrr$hrr_bpm[1],
       hrr2 = if (nrow(hrr) >= 2L) hrr$hrr_bpm[2] else NA_real_,
       slope_bpm_per_min = slope, hrr_1min = hrr_1min, abnormal = abnormal)
}
