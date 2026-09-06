#' Smoothness-priors detrending
#'
#' Removes slow trends from a uniformly sampled series (e.g. an interpolated
#' RR tachogram) using the smoothness-priors method of Tarvainen et al. (2002).
#' The trend estimate is \eqn{(I + \lambda^2 D_2^\top D_2)^{-1} z} where
#' \eqn{D_2} is the second-order difference operator, and the detrended
#' (stationary) component is \eqn{z_{stat} = (I - (I + \lambda^2
#' D_2^\top D_2)^{-1}) z}. Larger \code{lambda} lowers the cutoff frequency,
#' removing more low-frequency content.
#'
#' @param series Numeric vector, uniformly sampled.
#' @param lambda Smoothing parameter (regularisation). Higher values remove
#'   slower trends (lower cutoff). Default 500.
#' @param fs Sampling rate of \code{series} in Hz, used only to report the
#'   approximate half-power cutoff frequency. Default 4.
#' @return The detrended series (numeric vector, near-zero mean), carrying
#'   attributes \code{lambda} and \code{cutoff_hz} (the approximate half-power
#'   cutoff of the implied high-pass filter).
#' @references Tarvainen, M.P., Ranta-aho, P.O. & Karjalainen, P.A. (2002). An
#'   advanced detrending method with application to HRV analysis. \emph{IEEE
#'   Transactions on Biomedical Engineering}, 49(2), 172-175.
#' @seealso \code{\link{ecgHRVfreq}}
#' @export
#' @examples
#' t <- seq(0, 60, by = 0.25)
#' x <- 0.01 * t + sin(2 * pi * 0.25 * t)   # slow trend + HF oscillation
#' xd <- smoothnessPriorsDetrend(x, lambda = 500, fs = 4)
#' abs(mean(xd)) < 1e-6
smoothnessPriorsDetrend <- function(series, lambda = 500, fs = 4) {
  series <- as.numeric(series)
  n <- length(series)
  if (n < 3) {
    out <- series - mean(series)
    attr(out, "lambda") <- lambda
    attr(out, "cutoff_hz") <- NA_real_
    return(out)
  }

  # Second-order difference matrix D2: (n-2) x n
  D2 <- matrix(0, n - 2L, n)
  for (i in seq_len(n - 2L)) {
    D2[i, i]     <- 1
    D2[i, i + 1] <- -2
    D2[i, i + 2] <- 1
  }
  A <- diag(n) + lambda^2 * crossprod(D2)   # I + lambda^2 D2' D2
  trend <- solve(A, series)                 # (I + lambda^2 D2'D2)^-1 z
  stationary <- series - trend

  # True -3 dB (half-power) cutoff of the implied high-pass filter. The
  # low-pass (trend) amplitude gain is L(w) = 1/(1 + lambda^2 s^2) with
  # s = 2 - 2 cos(w); half power means |1 - L|^2 = 1/2, i.e.
  # s = 1 / (lambda * sqrt(sqrt(2) - 1)) and cos(w) = 1 - s/2.
  k <- sqrt(sqrt(2) - 1)
  arg <- 1 - 1 / (2 * lambda * k)
  cutoff_hz <- if (lambda > 0 && abs(arg) <= 1) {
    fs * acos(arg) / (2 * pi)
  } else {
    NA_real_
  }
  attr(stationary, "lambda") <- lambda
  attr(stationary, "cutoff_hz") <- cutoff_hz
  stationary
}
