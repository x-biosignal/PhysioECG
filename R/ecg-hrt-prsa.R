.hrt_rr_vector <- function(rr) {
  if (is.data.frame(rr)) {
    stopifnot("rr_ms" %in% names(rr))
    rr$rr_ms
  } else {
    as.numeric(rr)
  }
}

# Turbulence slope: the steepest positive regression slope over any 5 successive
# points of the post-PVC sinus tachogram.
.hrt_turbulence_slope <- function(post) {
  n <- length(post)
  if (n < 5L) return(NA_real_)
  x <- seq_len(5L)
  best <- -Inf
  for (i in seq_len(n - 4L)) {
    y <- post[i:(i + 4L)]
    if (anyNA(y)) next
    slope <- stats::cov(x, y) / stats::var(x)
    if (slope > best) best <- slope
  }
  if (is.finite(best)) best else NA_real_
}

#' Heart-Rate Turbulence (HRT)
#'
#' Computes turbulence onset (TO) and turbulence slope (TS), the biphasic
#' sinus-rhythm response to a ventricular premature beat (Schmidt et al. 1999).
#' RR intervals are indexed by beat so that \code{rr[i]} is the interval between
#' beat \code{i} and beat \code{i+1}; for a PVC at beat \code{p} the coupling
#' interval is \code{rr[p-1]} and the compensatory pause \code{rr[p]}, which are
#' excluded. The local tachograms of all usable PVCs are averaged before TO/TS.
#'
#' @param rr RR intervals (ms): a numeric vector or a data.frame with an
#'   \code{rr_ms} column.
#' @param pvc_index Integer beat indices of the ventricular premature beats (e.g.
#'   the \code{beat} column of \code{\link{ecgDetectPVC}} where \code{is_pvc}).
#' @param n_post Number of post-PVC sinus intervals used for TS (default 15).
#' @return A list with \code{to} (turbulence onset, percent; negative is normal),
#'   \code{ts} (turbulence slope, ms per RR interval; positive is normal),
#'   \code{n_pvc} (PVCs averaged) and \code{tachogram} (the averaged post-PVC
#'   sinus RR series).
#' @seealso \code{\link{ecgDetectPVC}}, \code{\link{ecgPRSA}}
#' @references Schmidt, G. et al. (1999). Heart-rate turbulence after ventricular
#'   premature beats as a predictor of mortality after acute myocardial
#'   infarction. \emph{Lancet}, 353(9162), 1390-1396.
#' @export
#' @examples
#' # canonical post-PVC pattern: early acceleration then deceleration
#' rr <- c(rep(800, 5), 500, 1100, 780, 785, 790 + 5 * (0:12))
#' ecgHRT(rr, pvc_index = 7)
ecgHRT <- function(rr, pvc_index, n_post = 15L) {
  x <- .hrt_rr_vector(rr)
  n <- length(x)
  pvc_index <- as.integer(pvc_index)

  to_vals <- numeric(0)
  post_mat <- list()
  pre_mat <- list()
  for (p in pvc_index) {
    if (p - 3L < 1L || p + n_post > n) next     # need clean surrounds
    pre <- c(x[p - 3L], x[p - 2L])              # RR(-2), RR(-1)
    post <- x[(p + 1L):(p + n_post)]            # RR(1..n_post)
    to_vals <- c(to_vals,
                 100 * ((post[1] + post[2]) - sum(pre)) / sum(pre))
    pre_mat[[length(pre_mat) + 1]] <- pre
    post_mat[[length(post_mat) + 1]] <- post
  }
  if (!length(post_mat)) {
    return(list(to = NA_real_, ts = NA_real_, n_pvc = 0L,
                tachogram = rep(NA_real_, n_post)))
  }
  avg_post <- colMeans(do.call(rbind, post_mat))
  avg_pre <- colMeans(do.call(rbind, pre_mat))
  to_avg <- 100 * ((avg_post[1] + avg_post[2]) - sum(avg_pre)) / sum(avg_pre)

  list(to = to_avg, ts = .hrt_turbulence_slope(avg_post),
       n_pvc = length(post_mat), tachogram = avg_post,
       to_per_pvc = to_vals)
}

#' Phase-Rectified Signal Averaging (PRSA): DC and AC
#'
#' Computes the deceleration capacity (DC) and acceleration capacity (AC) of the
#' RR series by phase-rectified signal averaging (Bauer et al. 2006). Anchors are
#' points where the RR level increases (DC) or decreases (AC); segments centered
#' on the anchors are averaged into a PRSA curve, and the capacity is a
#' Haar-wavelet coefficient of that curve at scale \code{s}.
#'
#' @param rr RR intervals (ms): a numeric vector or a data.frame with an
#'   \code{rr_ms} column.
#' @param T Averaging window (beats) for anchor selection (default 1: an anchor
#'   is a single beat larger/smaller than the previous).
#' @param s Wavelet scale for the DC/AC coefficient (default 2).
#' @param L Half-length of the averaged PRSA segments (default 15; must be
#'   \eqn{\ge s}).
#' @return A list with \code{dc} (deceleration capacity, ms; positive with
#'   dominant decelerations), \code{ac} (acceleration capacity, ms; negative),
#'   the PRSA curves \code{prsa_dc}/\code{prsa_ac}, the anchor counts
#'   \code{n_dc}/\code{n_ac}, and the parameters \code{T}, \code{s}, \code{L}.
#' @seealso \code{\link{ecgHRT}}, \code{\link{ecgRRintervals}}
#' @references Bauer, A. et al. (2006). Deceleration capacity of heart rate as a
#'   predictor of mortality after myocardial infarction: cohort study.
#'   \emph{Lancet}, 367(9523), 1674-1681.
#' @export
#' @examples
#' set.seed(1)
#' rr <- 800 + cumsum(rnorm(600, sd = 3))       # RR random walk
#' prsa <- ecgPRSA(rr)
#' c(DC = prsa$dc, AC = prsa$ac)
ecgPRSA <- function(rr, T = 1L, s = 2L, L = 15L) {
  x <- .hrt_rr_vector(rr)
  x <- x[is.finite(x)]
  n <- length(x)
  T <- as.integer(T); s <- as.integer(s); L <- as.integer(L)
  stopifnot(T >= 1L, s >= 1L, L >= s, n > 2L * L + 2L * T)

  # anchor selection: compare the mean of the T beats ending at i with the T
  # beats before it.
  before_mean <- after_mean <- rep(NA_real_, n)
  for (i in (T + 1L):(n - T + 1L)) {
    after_mean[i] <- mean(x[i:(i + T - 1L)])
    before_mean[i] <- mean(x[(i - T):(i - 1L)])
  }
  is_dc <- which(after_mean > before_mean)
  is_ac <- which(after_mean < before_mean)

  prsa_curve <- function(anchors) {
    anchors <- anchors[anchors - L >= 1L & anchors + L <= n]
    if (!length(anchors)) return(list(curve = rep(NA_real_, 2L * L + 1L), n = 0L))
    seg <- vapply(anchors, function(a) x[(a - L):(a + L)], numeric(2L * L + 1L))
    list(curve = rowMeans(seg), n = length(anchors))
  }
  dc_p <- prsa_curve(is_dc)
  ac_p <- prsa_curve(is_ac)

  capacity <- function(X) {
    if (anyNA(X)) return(NA_real_)
    ctr <- L + 1L
    (sum(X[ctr:(ctr + s - 1L)]) - sum(X[(ctr - s):(ctr - 1L)])) / (2 * s)
  }

  list(dc = capacity(dc_p$curve), ac = capacity(ac_p$curve),
       prsa_dc = dc_p$curve, prsa_ac = ac_p$curve,
       n_dc = dc_p$n, n_ac = ac_p$n, T = T, s = s, L = L)
}
