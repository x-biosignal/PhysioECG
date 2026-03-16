#' HRV Time-Domain Metrics
#'
#' Compute standard heart rate variability (HRV) time-domain metrics from
#' RR interval data. Calculates SDNN, RMSSD, pNN50, mean RR interval, and
#' mean heart rate for each channel.
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @return A data.frame with one row per channel and the following columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{mean_rr}{Mean RR interval in milliseconds.}
#'     \item{sdnn}{Standard deviation of all RR intervals (ms), reflecting
#'       overall HRV.}
#'     \item{rmssd}{Root mean square of successive RR interval differences
#'       (ms), reflecting short-term vagal modulation.}
#'     \item{pnn50}{Percentage of successive RR intervals differing by more
#'       than 50 ms.}
#'     \item{mean_hr}{Mean heart rate in beats per minute (60000 / mean_rr).}
#'   }
#'
#' @references Task Force of the European Society of Cardiology and the North
#'   American Society of Pacing and Electrophysiology (1996). "Heart rate
#'   variability: Standards of measurement, physiological interpretation and
#'   clinical use." \emph{Circulation}, 93(5), 1043--1065.
#'
#' @seealso \code{\link{ecgRRintervals}} for computing RR intervals,
#'   \code{\link{ecgHRVfreq}} for frequency-domain HRV analysis,
#'   \code{\link{ecgHRVnonlinear}} for nonlinear HRV metrics,
#'   \code{\link{ecgQualityCheck}} for ectopic beat detection before analysis.
#'
#' @export
ecgHRVtime <- function(rr) {
  stopifnot(is.data.frame(rr))
  required <- c("channel", "rr_ms", "time_sec")
  missing_cols <- setdiff(required, names(rr))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  channels <- unique(rr$channel)
  results <- list()

  for (ch in channels) {
    ch_data <- rr[rr$channel == ch, ]
    rr_ms <- ch_data$rr_ms

    if (length(rr_ms) < 10) {
      warning(sprintf("Channel %s has fewer than 10 RR intervals; HRV metrics may be unreliable.", ch))
    }

    diffs <- diff(rr_ms)

    mean_rr <- mean(rr_ms)
    sdnn <- sd(rr_ms)
    rmssd <- sqrt(mean(diffs^2))
    pnn50 <- if (length(diffs) > 0) {
      100 * sum(abs(diffs) > 50) / length(diffs)
    } else {
      NA_real_
    }
    mean_hr <- 60000 / mean_rr

    results[[length(results) + 1]] <- data.frame(
      channel = ch,
      mean_rr = mean_rr,
      sdnn = sdnn,
      rmssd = rmssd,
      pnn50 = pnn50,
      mean_hr = mean_hr,
      stringsAsFactors = FALSE)
  }

  do.call(rbind, results)
}
