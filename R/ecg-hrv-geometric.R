#' Geometric HRV indices (triangular index and TINN)
#'
#' Computes the geometric heart-rate-variability measures defined by the
#' Task Force (1996) from a series of NN (normal-to-normal) intervals: the
#' HRV triangular index and the triangular interpolation of the NN interval
#' histogram (TINN).
#'
#' @param rr Numeric vector of NN/RR intervals in milliseconds, or a data frame
#'   with an \code{rr_ms} column (as produced by the RR-interval readers); a data
#'   frame is reduced to its \code{rr_ms} column.
#' @param bin_ms Histogram bin width in milliseconds. Defaults to
#'   \code{7.8125} ms (\code{1/128} s), the Task Force standard.
#' @return A list with:
#'   \describe{
#'     \item{HRV_triangular_index}{Total number of NN intervals divided by the
#'       height (modal count) of the NN histogram.}
#'     \item{TINN}{Baseline width, in milliseconds, of the triangle that best
#'       approximates the NN histogram in a least-squares sense.}
#'   }
#' @references Task Force of the ESC and NASPE (1996). Heart rate variability:
#'   standards of measurement. \emph{Circulation}, 93(5), 1043-1065.
#' @seealso \code{\link{ecgHRVsegmented}}, \code{\link{ecgHRVtime}}
#' @export
#' @examples
#' set.seed(1)
#' rr <- 800 + rnorm(500, sd = 40)  # NN intervals in ms
#' ecgHRVgeometric(rr)
ecgHRVgeometric <- function(rr, bin_ms = 7.8125) {
  if (is.data.frame(rr)) {                       # accept an RR data frame (channel, rr_ms, time_sec)
    if (!"rr_ms" %in% names(rr)) stop("`rr` data frame must have an `rr_ms` column.", call. = FALSE)
    rr <- rr$rr_ms
  }
  stopifnot(is.numeric(rr), is.numeric(bin_ms), length(bin_ms) == 1, bin_ms > 0)
  rr <- rr[is.finite(rr)]
  n <- length(rr)
  if (n < 2) {
    return(list(HRV_triangular_index = NA_real_, TINN = NA_real_))
  }

  # Regular-grid histogram at bin_ms resolution
  lo <- floor(min(rr) / bin_ms) * bin_ms
  hi <- ceiling(max(rr) / bin_ms) * bin_ms
  breaks <- seq(lo - bin_ms, hi + bin_ms, by = bin_ms)
  h <- graphics::hist(rr, breaks = breaks, plot = FALSE)
  counts <- h$counts
  centers <- h$mids

  Y <- max(counts)                       # modal height
  tri_index <- n / Y                     # HRV triangular index
  mode_i <- which.max(counts)
  X <- centers[mode_i]                   # modal location

  # TINN: base width of the least-squares triangular interpolation. The
  # triangle peaks at (X, Y) and falls linearly to zero at feet N (< X) and
  # M (> X); pick N, M minimising the squared error to the histogram.
  tri_sse <- function(N, M) {
    if (N >= X || M <= X) return(NA_real_)
    q <- numeric(length(centers))
    l <- centers >= N & centers <= X
    r <- centers > X & centers <= M
    q[l] <- Y * (centers[l] - N) / (X - N)
    q[r] <- Y * (M - centers[r]) / (M - X)
    sum((counts - q)^2)
  }

  best <- Inf
  TINN <- NA_real_
  for (li in seq_len(mode_i)) {
    for (ri in mode_i:length(centers)) {
      s <- tri_sse(centers[li], centers[ri])
      if (!is.na(s) && s < best) {
        best <- s
        TINN <- centers[ri] - centers[li]
      }
    }
  }
  if (!is.finite(best)) {                # mode at an edge: fall back to support
    nz <- which(counts > 0)
    TINN <- if (length(nz) >= 2) centers[max(nz)] - centers[min(nz)] else bin_ms
  }

  list(HRV_triangular_index = tri_index, TINN = TINN)
}

#' Segmented (long-term) HRV indices: SDANN and SDNN index
#'
#' Splits an NN-interval series into consecutive fixed-duration segments
#' (5 minutes by default) and returns the two Task Force (1996) long-term
#' variability measures.
#'
#' @param rr Numeric vector of NN/RR intervals in milliseconds.
#' @param segment_sec Segment length in seconds (default 300 = 5 minutes).
#' @return A list with:
#'   \describe{
#'     \item{SDANN}{Standard deviation of the per-segment mean NN intervals
#'       (ms) - long-term variability.}
#'     \item{SDNN_index}{Mean of the per-segment SDNN values (ms) - average
#'       short-term variability.}
#'     \item{n_segments}{Number of segments used.}
#'   }
#' @references Task Force of the ESC and NASPE (1996). \emph{Circulation},
#'   93(5), 1043-1065.
#' @seealso \code{\link{ecgHRVgeometric}}, \code{\link{ecgHRVtime}}
#' @export
#' @examples
#' set.seed(1)
#' rr <- 800 + rnorm(3000, sd = 40)
#' ecgHRVsegmented(rr, segment_sec = 300)
ecgHRVsegmented <- function(rr, segment_sec = 300) {
  stopifnot(is.numeric(rr), is.numeric(segment_sec), length(segment_sec) == 1,
            segment_sec > 0)
  rr <- rr[is.finite(rr)]
  n <- length(rr)
  if (n < 2) {
    return(list(SDANN = NA_real_, SDNN_index = NA_real_, n_segments = 0L))
  }

  # Cumulative time (s) at the end of each interval -> segment membership
  t_cum <- cumsum(rr) / 1000
  seg_id <- floor(t_cum / segment_sec)
  segs <- split(rr, seg_id)

  # Task Force: compute over COMPLETE segments only; drop a short trailing
  # segment so its transient mean does not inflate SDANN.
  seg_dur <- vapply(segs, function(v) sum(v) / 1000, numeric(1))
  segs <- segs[seg_dur >= 0.8 * segment_sec]
  if (length(segs) < 1) {
    return(list(SDANN = NA_real_, SDNN_index = NA_real_, n_segments = 0L))
  }

  seg_means <- vapply(segs, mean, numeric(1))
  seg_sds <- vapply(segs, function(v) if (length(v) >= 2) stats::sd(v) else NA_real_,
                    numeric(1))

  list(
    SDANN = if (length(seg_means) >= 2) stats::sd(seg_means) else NA_real_,
    SDNN_index = mean(seg_sds, na.rm = TRUE),
    n_segments = length(segs)
  )
}
