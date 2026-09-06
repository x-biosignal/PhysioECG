#' Time-Varying HRV (Sliding-Window Time and Frequency Metrics)
#'
#' Computes heart-rate-variability trajectories by sliding a time window over the
#' RR series and evaluating time-domain (SDNN, RMSSD) and frequency-domain
#' (LF, HF, LF/HF) metrics in each window, reusing \code{\link{ecgHRVtime}} and
#' \code{\link{ecgHRVfreq}} (Task Force 1996; Mainardi 2009).
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms} and
#'   \code{time_sec} (as returned by \code{\link{ecgRRintervals}}), assumed
#'   ordered in time per channel.
#' @param window_sec Analysis window length in seconds (default 300).
#' @param step_sec Step between successive windows in seconds (default 30).
#' @param freq_method Spectral method for \code{\link{ecgHRVfreq}}: "ar"
#'   (default), "welch" or "lomb".
#' @param detrend,detrend_lambda Passed to \code{\link{ecgHRVfreq}} for optional
#'   smoothness-priors detrending of the resampled tachogram.
#' @param min_beats Minimum RR intervals in a window to compute metrics; windows
#'   with fewer beats yield \code{NA} metrics (default 20).
#' @param rhythm_check Passed to the per-window HRV functions (default FALSE, so
#'   the trajectory is continuous and not gated by the AF detector).
#' @return A data.frame with one row per window per channel and columns
#'   \code{channel}, \code{window}, \code{time_start}, \code{time_center},
#'   \code{n_beats}, \code{sdnn}, \code{rmssd}, \code{lf}, \code{hf} and
#'   \code{lf_hf}.
#' @seealso \code{\link{ecgHRVtime}}, \code{\link{ecgHRVfreq}},
#'   \code{\link{ecgRRintervals}}
#' @references Task Force of the European Society of Cardiology and the North
#'   American Society of Pacing and Electrophysiology (1996). Heart rate
#'   variability: standards of measurement, physiological interpretation, and
#'   clinical use. \emph{Circulation}, 93(5), 1043-1065. Mainardi, L.T. (2009).
#'   On the quantification of heart rate variability spectral parameters using
#'   time-frequency and time-varying methods. \emph{Philosophical Transactions
#'   of the Royal Society A}, 367(1887), 255-275.
#' @export
#' @examples
#' set.seed(1)
#' n <- 600
#' rr <- data.frame(channel = 1L,
#'                  rr_ms = 800 + 25 * sin(2 * pi * 0.1 * cumsum(rep(0.8, n))) +
#'                    rnorm(n, sd = 10))
#' rr$time_sec <- cumsum(rr$rr_ms) / 1000
#' traj <- ecgHRVtimevarying(rr, window_sec = 120, step_sec = 30)
#' head(traj)
ecgHRVtimevarying <- function(rr, window_sec = 300, step_sec = 30,
                              freq_method = c("ar", "welch", "lomb"),
                              detrend = FALSE, detrend_lambda = 500,
                              min_beats = 20L, rhythm_check = FALSE) {
  freq_method <- match.arg(freq_method)
  stopifnot(is.data.frame(rr))
  req <- c("channel", "rr_ms", "time_sec")
  miss <- setdiff(req, names(rr))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  stopifnot(window_sec > 0, step_sec > 0)

  na_row <- function(ch, k, ws, nb) {
    data.frame(channel = as.integer(ch), window = k, time_start = ws,
               time_center = ws + window_sec / 2, n_beats = nb,
               sdnn = NA_real_, rmssd = NA_real_, lf = NA_real_, hf = NA_real_,
               lf_hf = NA_real_, stringsAsFactors = FALSE)
  }

  out <- list()
  for (ch in unique(rr$channel)) {
    d <- rr[rr$channel == ch, , drop = FALSE]
    d <- d[order(d$time_sec), , drop = FALSE]
    t0 <- min(d$time_sec); tn <- max(d$time_sec)
    starts <- if (tn - t0 >= window_sec) {
      seq(t0, tn - window_sec, by = step_sec)
    } else {
      t0
    }

    for (k in seq_along(starts)) {
      ws <- starts[k]
      win <- d[d$time_sec >= ws & d$time_sec < ws + window_sec, , drop = FALSE]
      nb <- nrow(win)
      if (nb < min_beats) {
        out[[length(out) + 1]] <- na_row(ch, k, ws, nb)
        next
      }
      tw <- suppressWarnings(ecgHRVtime(win, rhythm_check = rhythm_check))
      fw <- suppressWarnings(ecgHRVfreq(win, method = freq_method,
                                        detrend = detrend,
                                        detrend_lambda = detrend_lambda,
                                        rhythm_check = rhythm_check))
      out[[length(out) + 1]] <- data.frame(
        channel = as.integer(ch), window = k, time_start = ws,
        time_center = ws + window_sec / 2, n_beats = nb,
        sdnn = tw$sdnn, rmssd = tw$rmssd, lf = fw$lf, hf = fw$hf,
        lf_hf = fw$lf_hf_ratio, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}
