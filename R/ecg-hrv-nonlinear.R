#' Poincare Plot Descriptors for HRV Analysis
#'
#' Computes Poincare plot descriptors (SD1, SD2, SD1/SD2 ratio) from
#' RR interval data. SD1 reflects short-term variability (perpendicular
#' to the identity line), while SD2 reflects long-term variability
#' (along the identity line).
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @return A data.frame with one row per channel and the following columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{sd1}{Standard deviation perpendicular to the identity line (ms),
#'       reflecting beat-to-beat (short-term) variability.}
#'     \item{sd2}{Standard deviation along the identity line (ms), reflecting
#'       long-term variability.}
#'     \item{sd1_sd2_ratio}{Ratio of SD1 to SD2, or \code{NA} if SD2 is
#'       zero.}
#'   }
#'
#' @references Task Force of the European Society of Cardiology and the North
#'   American Society of Pacing and Electrophysiology (1996). "Heart rate
#'   variability: Standards of measurement, physiological interpretation and
#'   clinical use." \emph{Circulation}, 93(5), 1043--1065.
#'
#' @seealso \code{\link{ecgHRVnonlinear}} for the combined nonlinear analysis
#'   wrapper, \code{\link{ecgSampleEntropy}} for sample entropy,
#'   \code{\link{ecgDFA}} for detrended fluctuation analysis.
#'
#' @export
ecgHRVpoincare <- function(rr) {
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

    if (length(rr_ms) < 2) {
      warning(sprintf("Channel %s has fewer than 2 RR intervals; Poincare descriptors cannot be computed.", ch))
      results[[length(results) + 1]] <- data.frame(
        channel = ch,
        sd1 = NA_real_,
        sd2 = NA_real_,
        sd1_sd2_ratio = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    sdsd_sq <- var(diff(rr_ms))
    sdnn_sq <- var(rr_ms)

    sd1 <- sqrt(0.5 * sdsd_sq)
    sd2_sq <- 2 * sdnn_sq - 0.5 * sdsd_sq
    sd2 <- if (sd2_sq > 0) sqrt(sd2_sq) else 0
    sd1_sd2_ratio <- if (sd2 > 0) sd1 / sd2 else NA_real_

    results[[length(results) + 1]] <- data.frame(
      channel = ch,
      sd1 = sd1,
      sd2 = sd2,
      sd1_sd2_ratio = sd1_sd2_ratio,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Heart-Rate Asymmetry (Poincare Plot Asymmetry) Descriptors
#'
#' Computes heart-rate asymmetry (HRA) descriptors from RR interval data.
#' The Poincare plot is asymmetric about the line of identity because heart-rate
#' \emph{decelerations} (prolongations of RR, points above the line, \eqn{RR_{n+1}
#' > RR_n}) and \emph{accelerations} (shortenings, points below the line)
#' contribute unequally to variability -- a signature of vagal activity and of the
#' temporal irreversibility of heart-rate dynamics. Decelerations typically
#' dominate short-term variability while accelerations dominate long-term
#' variability (Piskorski & Guzik, 2011).
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @return A data.frame with one row per channel and the columns:
#'   \describe{
#'     \item{channel}{Channel identifier.}
#'     \item{gi}{Guzik's Index (\%): the cumulative distance of the deceleration
#'       points from the line of identity divided by that of all points.}
#'     \item{si}{Slope Index (\%): the cumulative phase angle of the deceleration
#'       points divided by that of all points.}
#'     \item{ai}{Area Index (\%): the cumulative sector area of the deceleration
#'       points divided by that of all points.}
#'     \item{pi}{Porta's Index (\%): the number of acceleration points divided by
#'       the number of points off the line of identity.}
#'     \item{c1d, c1a}{Contributions of decelerations / accelerations to
#'       short-term HRV (\code{c1d + c1a = 1}).}
#'     \item{sd1d, sd1a}{Short-term variance components of decelerations /
#'       accelerations (ms).}
#'     \item{c2d, c2a}{Contributions of decelerations / accelerations to
#'       long-term HRV (\code{c2d + c2a = 1}).}
#'     \item{sd2d, sd2a}{Long-term variance components (ms).}
#'     \item{cd, ca}{Total contributions of decelerations / accelerations to HRV.}
#'     \item{sdnnd, sdnna}{Total variance components (ms).}
#'   }
#'   A value indicating asymmetry is \code{c1d > c1a} (decelerations dominate
#'   short-term variability). Columns are \code{NA} for channels with fewer than
#'   three RR intervals.
#'
#' @references
#' Guzik, P., et al. (2006). "Heart rate asymmetry by Poincare plots of RR
#'   intervals." \emph{Biomedizinische Technik}, 51(4), 272--275.
#'
#' Piskorski, J., & Guzik, P. (2011). "Asymmetric properties of long-term and
#'   total heart rate variability." \emph{Medical & Biological Engineering &
#'   Computing}, 49(11), 1289--1297. \doi{10.1007/s11517-011-0834-z}
#'
#' Porta, A., et al. (2008). "Temporal asymmetries of short-term heart period
#'   variability are linked to autonomic regulation." \emph{American Journal of
#'   Physiology}, 295(2), R550--R557.
#'
#' @seealso \code{\link{ecgHRVpoincare}} for the symmetric Poincare descriptors
#'   (SD1/SD2), \code{\link{ecgDFA}} for detrended fluctuation analysis.
#'
#' @export
ecgHRVasymmetry <- function(rr) {
  stopifnot(is.data.frame(rr))
  required <- c("channel", "rr_ms", "time_sec")
  missing_cols <- setdiff(required, names(rr))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  na_row <- function(ch) data.frame(
    channel = ch, gi = NA_real_, si = NA_real_, ai = NA_real_, pi = NA_real_,
    c1d = NA_real_, c1a = NA_real_, sd1d = NA_real_, sd1a = NA_real_,
    c2d = NA_real_, c2a = NA_real_, sd2d = NA_real_, sd2a = NA_real_,
    cd = NA_real_, ca = NA_real_, sdnnd = NA_real_, sdnna = NA_real_,
    stringsAsFactors = FALSE)

  channels <- unique(rr$channel)
  results <- list()

  for (ch in channels) {
    rr_ms <- rr$rr_ms[rr$channel == ch]
    n <- length(rr_ms)
    if (n < 3) {
      warning(sprintf("Channel %s has fewer than 3 RR intervals; HRA cannot be computed.", ch))
      results[[length(results) + 1]] <- na_row(ch)
      next
    }

    N <- n - 1L
    x <- rr_ms[-n]           # RR_n   (x-axis)
    y <- rr_ms[-1]           # RR_n+1 (y-axis)
    d <- y - x
    decel <- which(d > 0)    # points above the line of identity (prolongations)
    accel <- which(d < 0)    # points below the line of identity (shortenings)
    noch  <- which(d == 0)

    # distance to the line of identity, distance to the centroid line, phase angle,
    # and sector area (Piskorski & Guzik 2011; matches NeuroKit2's hrv_nonlinear)
    cx <- mean(x); cy <- mean(y)
    dist_l2 <- abs((x - cx) + (y - cy)) / sqrt(2)
    dist    <- abs(y - x) / sqrt(2)
    theta   <- abs(atan(1) - atan(y / x))
    S       <- 0.5 * theta * (x^2 + y^2)

    gi <- sum(dist[decel]) / sum(dist) * 100
    si <- sum(theta[decel]) / sum(theta) * 100
    ai <- sum(S[decel]) / sum(S) * 100
    pi_idx <- length(accel) / (N - length(noch)) * 100

    sd1d <- sqrt(sum(dist[decel]^2) / (N - 1))
    sd1a <- sqrt(sum(dist[accel]^2) / (N - 1))
    sd1I <- sqrt(sd1d^2 + sd1a^2)
    c1d <- (sd1d / sd1I)^2; c1a <- (sd1a / sd1I)^2

    ltd <- sum(dist_l2[decel]^2) / (N - 1)
    lta <- sum(dist_l2[accel]^2) / (N - 1)
    ltn <- sum(dist_l2[noch]^2) / (N - 1)
    sd2d <- sqrt(ltd + 0.5 * ltn); sd2a <- sqrt(lta + 0.5 * ltn)
    sd2I <- sqrt(sd2d^2 + sd2a^2)
    c2d <- (sd2d / sd2I)^2; c2a <- (sd2a / sd2I)^2

    sdnnd <- sqrt(0.5 * (sd1d^2 + sd2d^2))
    sdnna <- sqrt(0.5 * (sd1a^2 + sd2a^2))
    sdnn <- sqrt(sdnnd^2 + sdnna^2)
    cd <- (sdnnd / sdnn)^2; ca <- (sdnna / sdnn)^2

    results[[length(results) + 1]] <- data.frame(
      channel = ch, gi = gi, si = si, ai = ai, pi = pi_idx,
      c1d = c1d, c1a = c1a, sd1d = sd1d, sd1a = sd1a,
      c2d = c2d, c2a = c2a, sd2d = sd2d, sd2a = sd2a,
      cd = cd, ca = ca, sdnnd = sdnnd, sdnna = sdnna,
      stringsAsFactors = FALSE)
  }

  do.call(rbind, results)
}


#' Cardiac Autonomic Indices (Toichi CSI / CVI) from the Poincare Plot
#'
#' Computes the cardiac sympathetic and vagal indices (Toichi et al., 1997) from
#' RR interval data. The indices are read off the Poincare plot's \eqn{4\times}SD
#' bounding box: the long axis \eqn{L = 4\,\mathrm{SD2}} and the short axis
#' \eqn{T = 4\,\mathrm{SD1}}.
#' \describe{
#'   \item{\code{csi}}{Cardiac Sympathetic Index, \eqn{L/T} (= SD2/SD1); rises with
#'     sympathetic predominance.}
#'   \item{\code{cvi}}{Cardiac Vagal Index, \eqn{\log_{10}(L\times T)}; rises with
#'     parasympathetic (vagal) activity.}
#'   \item{\code{csi_modified}}{Modified CSI, \eqn{L^2/T}, a more sensitive
#'     sympathetic marker (Jeppesen et al., 2014) used in seizure detection.}
#' }
#' Here SD1 is the universal beat-to-beat term (\eqn{\mathrm{SDSD}/\sqrt2}) and SD2
#' is the geometric (paired-projection) long-term term, matching the convention in
#' which the Toichi indices are defined (as in NeuroKit2's \code{hrv_nonlinear});
#' this SD2 differs by a small \eqn{O(1/n)} amount from the analytical closed-form
#' SD2 reported by \code{\link{ecgHRVpoincare}}.
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @return A data.frame with one row per channel and the columns \code{channel},
#'   \code{csi}, \code{cvi}, and \code{csi_modified}. Columns are \code{NA} for
#'   channels with fewer than three RR intervals.
#'
#' @references
#' Toichi, M., Sugiura, T., Murai, T., & Sengoku, A. (1997). "A new method of
#'   assessing cardiac autonomic function and its comparison with spectral analysis
#'   and coefficient of variation of R-R interval." \emph{Journal of the Autonomic
#'   Nervous System}, 62(1-2), 79--84.
#'
#' Jeppesen, J., et al. (2014). "Using Lorenz plot and Cardiac Sympathetic Index of
#'   heart rate variability for detecting seizures for patients with epilepsy."
#'   \emph{IEEE EMBC}, 4563--4566.
#'
#' @seealso \code{\link{ecgHRVpoincare}} for SD1/SD2, \code{\link{ecgHRVasymmetry}}
#'   for heart-rate asymmetry.
#'
#' @export
ecgHRVautonomic <- function(rr) {
  stopifnot(is.data.frame(rr))
  required <- c("channel", "rr_ms", "time_sec")
  missing_cols <- setdiff(required, names(rr))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  channels <- unique(rr$channel)
  results <- list()

  for (ch in channels) {
    rr_ms <- rr$rr_ms[rr$channel == ch]
    n <- length(rr_ms)
    if (n < 3) {
      warning(sprintf("Channel %s has fewer than 3 RR intervals; CSI/CVI cannot be computed.", ch))
      results[[length(results) + 1]] <- data.frame(
        channel = ch, csi = NA_real_, cvi = NA_real_, csi_modified = NA_real_,
        stringsAsFactors = FALSE)
      next
    }

    a <- rr_ms[-n]; b <- rr_ms[-1]
    sd1 <- sqrt(0.5 * var(diff(rr_ms)))          # SDSD / sqrt(2) (universal short-term term)
    sd2 <- sd((a + b) / sqrt(2))                 # geometric paired-projection long-term term
    T <- 4 * sd1; L <- 4 * sd2                    # Poincare 4-SD bounding box
    results[[length(results) + 1]] <- data.frame(
      channel = ch,
      csi = L / T,
      cvi = log10(L * T),
      csi_modified = L^2 / T,
      stringsAsFactors = FALSE)
  }

  do.call(rbind, results)
}


#' Sample Entropy of RR Intervals
#'
#' Computes sample entropy (SampEn) from RR interval data. Sample entropy
#' measures the regularity or predictability of a time series. Lower values
#' indicate more regular (predictable) signals, while higher values indicate
#' more complex (irregular) signals.
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @param m Embedding dimension (default: 2). Length of template patterns
#'   to compare.
#' @param r_factor Tolerance factor (default: 0.2). The tolerance \code{r}
#'   is computed as \code{r_factor * sd(rr_ms)}.
#' @return A data.frame with one row per channel and the following columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{sample_entropy}{Sample entropy value (nats). Lower values
#'       indicate more regular signals; higher values indicate more complex
#'       signals. \code{NA} if the series is too short or constant.}
#'     \item{m}{Embedding dimension used.}
#'     \item{r}{Tolerance threshold (ms) computed as \code{r_factor * sd(rr_ms)}.}
#'   }
#'
#' @references Richman, J.S. & Moorman, J.R. (2000). "Physiological
#'   time-series analysis using approximate entropy and sample entropy."
#'   \emph{American Journal of Physiology-Heart and Circulatory Physiology},
#'   278(6), H2039--H2049. \doi{10.1152/ajpheart.2000.278.6.H2039}
#'
#' @seealso \code{\link{ecgHRVnonlinear}} for the combined nonlinear analysis
#'   wrapper, \code{\link{ecgHRVpoincare}} for Poincare plot descriptors,
#'   \code{\link{ecgDFA}} for detrended fluctuation analysis.
#'
#' @export
ecgSampleEntropy <- function(rr, m = 2L, r_factor = 0.2) {
  stopifnot(is.data.frame(rr))
  required <- c("channel", "rr_ms", "time_sec")
  missing_cols <- setdiff(required, names(rr))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  stopifnot(is.numeric(m) && length(m) == 1 && m >= 1)
  stopifnot(is.numeric(r_factor) && length(r_factor) == 1 && r_factor > 0)

  channels <- unique(rr$channel)
  results <- list()

  for (ch in channels) {
    ch_data <- rr[rr$channel == ch, ]
    rr_ms <- ch_data$rr_ms
    n <- length(rr_ms)

    if (n < m + 2) {
      warning(sprintf("Channel %s has too few RR intervals (%d) for SampEn with m=%d.", ch, n, m))
      results[[length(results) + 1]] <- data.frame(
        channel = ch,
        sample_entropy = NA_real_,
        m = as.integer(m),
        r = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    r <- r_factor * sd(rr_ms)

    if (r == 0) {
      results[[length(results) + 1]] <- data.frame(
        channel = ch,
        sample_entropy = NA_real_,
        m = as.integer(m),
        r = r,
        stringsAsFactors = FALSE
      )
      next
    }

    B <- .count_matches(rr_ms, m, r)
    A <- .count_matches(rr_ms, m + 1L, r)

    sampen <- if (B > 0 && A > 0) -log(A / B) else NA_real_

    results[[length(results) + 1]] <- data.frame(
      channel = ch,
      sample_entropy = sampen,
      m = as.integer(m),
      r = r,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Count template matches for sample entropy
#'
#' Counts the number of template vector pairs of length \code{dim} that are
#' within tolerance \code{r} (Chebyshev distance). Uses a vectorised approach
#' with an embedding matrix to avoid nested R loops for better performance on
#' long RR series.
#'
#' @param x Numeric vector (time series).
#' @param dim Embedding dimension.
#' @param r Tolerance.
#' @return Integer count of matching template pairs.
#' @keywords internal
.count_matches <- function(x, dim, r) {
  n <- length(x)
  n_templates <- n - dim
  if (n_templates < 1) return(0L)

  # Build embedding matrix: each row is a template of length dim
  emb <- matrix(NA_real_, nrow = n_templates, ncol = dim)
  for (d in seq_len(dim)) {
    emb[, d] <- x[d:(d + n_templates - 1L)]
  }

  # Vectorised pairwise Chebyshev distance check
  # Process by column: for each lag dimension, compute absolute differences
  # between all pairs and track running maximum
  count <- 0L

  # For moderate n_templates, use column-wise vectorisation with outer indexing
  # For very large n_templates, process in blocks to limit memory
  block_size <- 2000L
  n_blocks <- ceiling(n_templates / block_size)

  for (bi in seq_len(n_blocks)) {
    i_start <- (bi - 1L) * block_size + 1L
    i_end <- min(bi * block_size, n_templates)
    rows_i <- i_start:i_end

    for (bj in seq(bi, n_blocks)) {
      j_start <- if (bj == bi) i_start else (bj - 1L) * block_size + 1L
      j_end <- min(bj * block_size, n_templates)
      rows_j <- j_start:j_end

      # Compute max absolute difference across all dim columns
      max_diff <- abs(outer(emb[rows_i, 1], emb[rows_j, 1], `-`))
      if (dim > 1) {
        for (d in 2:dim) {
          max_diff <- pmax(max_diff, abs(outer(emb[rows_i, d], emb[rows_j, d], `-`)))
        }
      }

      # Count pairs where max_diff < r (upper triangle only to avoid double-counting)
      if (bj == bi) {
        # Same block: only count upper triangle (i < j)
        upper <- upper.tri(max_diff, diag = FALSE)
        count <- count + sum(max_diff[upper] < r)
      } else {
        # Different blocks: all pairs are unique
        count <- count + sum(max_diff < r)
      }
    }
  }

  count
}


#' Detrended Fluctuation Analysis of RR Intervals
#'
#' Performs detrended fluctuation analysis (DFA) on RR interval data to
#' characterize fractal scaling properties. Computes alpha1 (short-range
#' correlations, 4--16 beats) and alpha2 (long-range correlations, 16--64
#' beats).
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @param short_range Numeric vector of length 2 defining the scale range
#'   (in beats) for alpha1 (default: c(4, 16)).
#' @param long_range Numeric vector of length 2 defining the scale range
#'   (in beats) for alpha2 (default: c(16, 64)).
#' @return A data.frame with one row per channel and the following columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{alpha1}{Short-range scaling exponent. Values near 1.0 indicate
#'       fractal-like (healthy) correlations; values near 0.5 indicate
#'       uncorrelated (random) behavior; values near 1.5 suggest
#'       Brownian noise. \code{NA} if the series is too short.}
#'     \item{alpha2}{Long-range scaling exponent with the same interpretation
#'       as alpha1 but over larger time scales. \code{NA} if the series is
#'       too short.}
#'   }
#'
#' @references Peng, C.-K., et al. (1994). "Mosaic organization of DNA
#'   nucleotides." \emph{Physical Review E}, 49(2), 1685--1689.
#'   \doi{10.1103/PhysRevE.49.1685}
#'
#' @seealso \code{\link{ecgHRVnonlinear}} for the combined nonlinear analysis
#'   wrapper, \code{\link{ecgSampleEntropy}} for sample entropy,
#'   \code{\link{ecgHRVpoincare}} for Poincare plot descriptors.
#'
#' @export
ecgDFA <- function(rr,
                   short_range = c(4, 16),
                   long_range = c(16, 64)) {
  stopifnot(is.data.frame(rr))
  required <- c("channel", "rr_ms", "time_sec")
  missing_cols <- setdiff(required, names(rr))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  stopifnot(is.numeric(short_range) && length(short_range) == 2)
  stopifnot(is.numeric(long_range) && length(long_range) == 2)

  channels <- unique(rr$channel)
  results <- list()

  for (ch in channels) {
    ch_data <- rr[rr$channel == ch, ]
    rr_ms <- ch_data$rr_ms
    n <- length(rr_ms)

    # Integrate: cumulative sum of mean-subtracted series
    y <- cumsum(rr_ms - mean(rr_ms))

    # Generate scales: LOGARITHMICALLY spaced integers from 4 up to n/4. DFA
    # requires log-spaced scales so both the short range (alpha1, ~4-16) and the
    # long range (alpha2, ~16-64) are populated; linear spacing over a large n
    # leaves those ranges with < 2 scales, so the slope fits return NA.
    max_scale <- max(5L, as.integer(n %/% 4))
    scales <- unique(as.integer(round(
      exp(seq(log(4), log(max_scale), length.out = 30L)))))
    scales <- scales[scales >= 4 & scales <= max_scale]

    if (length(scales) < 2) {
      warning(sprintf("Channel %s has too few RR intervals (%d) for DFA.", ch, n))
      results[[length(results) + 1]] <- data.frame(
        channel = ch,
        alpha1 = NA_real_,
        alpha2 = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    fluct <- numeric(length(scales))
    for (si in seq_along(scales)) {
      s <- scales[si]
      n_windows <- n %/% s
      if (n_windows < 1) {
        fluct[si] <- NA_real_
        next
      }

      rms_vals <- numeric(n_windows)
      for (w in seq_len(n_windows)) {
        idx_start <- (w - 1L) * s + 1L
        idx_end <- w * s
        segment <- y[idx_start:idx_end]

        # Linear detrending
        t_seg <- seq_len(s)
        fit <- .fast_lm(t_seg, segment)
        trend <- fit[1] + fit[2] * t_seg
        residuals <- segment - trend
        rms_vals[w] <- sqrt(mean(residuals^2))
      }

      fluct[si] <- mean(rms_vals)
    }

    # Remove any NA values
    valid <- !is.na(fluct) & fluct > 0
    log_scales <- log(scales[valid])
    log_fluct <- log(fluct[valid])

    # Compute alpha1 (short-range) and alpha2 (long-range)
    alpha1 <- .dfa_slope(log_scales, log_fluct, scales[valid], short_range)
    alpha2 <- .dfa_slope(log_scales, log_fluct, scales[valid], long_range)

    results[[length(results) + 1]] <- data.frame(
      channel = ch,
      alpha1 = alpha1,
      alpha2 = alpha2,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


#' Compute DFA slope for a given scale range
#'
#' @param log_scales Log of scale values.
#' @param log_fluct Log of fluctuation values.
#' @param scales Original scale values.
#' @param range Numeric vector of length 2 (min_scale, max_scale).
#' @return Numeric slope (alpha exponent).
#' @keywords internal
.dfa_slope <- function(log_scales, log_fluct, scales, range) {
  idx <- which(scales >= range[1] & scales <= range[2])
  if (length(idx) < 2) return(NA_real_)
  fit <- .fast_lm(log_scales[idx], log_fluct[idx])
  fit[2]  # slope
}


#' Fast simple linear regression (intercept + slope)
#'
#' @param x Predictor vector.
#' @param y Response vector.
#' @return Numeric vector of length 2: c(intercept, slope).
#' @keywords internal
.fast_lm <- function(x, y) {
  n <- length(x)
  sx <- sum(x)
  sy <- sum(y)
  sxy <- sum(x * y)
  sxx <- sum(x^2)
  denom <- n * sxx - sx^2
  if (abs(denom) < .Machine$double.eps) return(c(mean(y), 0))
  slope <- (n * sxy - sx * sy) / denom
  intercept <- (sy - slope * sx) / n
  c(intercept, slope)
}


#' Nonlinear HRV Analysis (Convenience Wrapper)
#'
#' Computes all nonlinear HRV metrics by calling \code{\link{ecgHRVpoincare}},
#' \code{\link{ecgSampleEntropy}}, and \code{\link{ecgDFA}}, and merging
#' the results into a single data.frame.
#'
#' @param rr A data.frame with columns \code{channel}, \code{rr_ms}, and
#'   \code{time_sec}, as returned by \code{\link{ecgRRintervals}}.
#' @param m Embedding dimension for sample entropy (default: 2).
#' @param r_factor Tolerance factor for sample entropy (default: 0.2).
#' @param short_range Scale range for DFA alpha1 (default: c(4, 16)).
#' @param long_range Scale range for DFA alpha2 (default: c(16, 64)).
#' @param rhythm_check Logical; if \code{TRUE} (default), gate the result on
#'   \code{\link{ecgDetectAF}}: metrics for channels classified as AF or
#'   frequent ectopy are set to \code{NA} and \code{rhythm}/\code{hrv_valid}
#'   columns are appended (HRV is undefined in atrial fibrillation).
#' @return A data.frame with columns: channel, sd1, sd2, sd1_sd2_ratio,
#'   sample_entropy, m, r, alpha1, alpha2.
#'
#' @references Shaffer, F. & Ginsberg, J.P. (2017). "An overview of heart rate
#'   variability metrics and norms." \emph{Frontiers in Public Health}, 5, 258.
#'   \doi{10.3389/fpubh.2017.00258}
#'
#'   Task Force of the European Society of Cardiology and the North American
#'   Society of Pacing and Electrophysiology (1996). "Heart rate variability:
#'   Standards of measurement, physiological interpretation and clinical use."
#'   \emph{Circulation}, 93(5), 1043--1065.
#'
#' @seealso \code{\link{ecgHRVpoincare}} for Poincare plot descriptors,
#'   \code{\link{ecgSampleEntropy}} for sample entropy,
#'   \code{\link{ecgDFA}} for detrended fluctuation analysis,
#'   \code{\link{ecgHRVtime}} for time-domain HRV metrics,
#'   \code{\link{ecgHRVfreq}} for frequency-domain HRV analysis.
#'
#' @export
ecgHRVnonlinear <- function(rr,
                            m = 2L,
                            r_factor = 0.2,
                            short_range = c(4, 16),
                            long_range = c(16, 64),
                            rhythm_check = TRUE) {
  poincare <- ecgHRVpoincare(rr)
  sampen <- ecgSampleEntropy(rr, m = m, r_factor = r_factor)
  dfa <- ecgDFA(rr, short_range = short_range, long_range = long_range)

  result <- merge(poincare, sampen, by = "channel")
  result <- merge(result, dfa, by = "channel")
  if (isTRUE(rhythm_check)) result <- .apply_rhythm_gate(result, rr)
  result
}

# Short-range DFA scaling exponent (alpha1) from a raw RR vector.
.dfa_alpha1 <- function(rr_vals, short_range = c(4, 16)) {
  rr_vals <- rr_vals[is.finite(rr_vals)]
  n <- length(rr_vals)
  if (n < max(8L, as.integer(short_range[2]))) return(NA_real_)

  y <- cumsum(rr_vals - mean(rr_vals))
  max_scale <- max(5L, as.integer(n %/% 4))
  # LOG-spaced scales (see ecgDFA): linear spacing under-populates the short range.
  scales <- unique(as.integer(round(
    exp(seq(log(4), log(max_scale), length.out = 30L)))))
  scales <- scales[scales >= 4 & scales <= max_scale]
  if (length(scales) < 2) return(NA_real_)

  fluct <- numeric(length(scales))
  for (si in seq_along(scales)) {
    s <- scales[si]
    nw <- n %/% s
    if (nw < 1) { fluct[si] <- NA_real_; next }
    rms <- numeric(nw)
    for (w in seq_len(nw)) {
      seg <- y[((w - 1L) * s + 1L):(w * s)]
      t_seg <- seq_len(s)
      fit <- .fast_lm(t_seg, seg)
      rms[w] <- sqrt(mean((seg - (fit[1] + fit[2] * t_seg))^2))
    }
    fluct[si] <- mean(rms)
  }

  valid <- !is.na(fluct) & fluct > 0
  .dfa_slope(log(scales[valid]), log(fluct[valid]), scales[valid], short_range)
}

#' Rolling (time-resolved) DFA alpha1
#'
#' Computes the short-range detrended-fluctuation scaling exponent (alpha1)
#' over a sliding window of beats, yielding a time-resolved alpha1 trajectory.
#' This is the basis of the DFA-a1 exercise-intensity thresholds of Gronwald &
#' Rogers (2020).
#'
#' @param rr A data.frame with columns \code{channel} and \code{rr_ms}.
#' @param window_beats Window length in beats (default 300).
#' @param step_beats Step between successive windows in beats (default 30).
#' @param short_range Scale range (in beats) for alpha1 (default c(4, 16)).
#' @return A data.frame with one row per window and columns \code{channel},
#'   \code{beat_start}, \code{beat_end}, \code{beat_center} and \code{alpha1}.
#' @references Peng, C.K. et al. (1995). Quantification of scaling exponents.
#'   \emph{Chaos}, 5(1), 82-87. Gronwald, T. & Rogers, B. (2020). Fractal
#'   correlation properties of heart rate variability as a biomarker.
#'   \emph{Frontiers in Physiology}, 11, 550572.
#' @seealso \code{\link{ecgDFA}}, \code{\link{dfaThresholdCrossings}}
#' @export
#' @examples
#' set.seed(1)
#' rr <- data.frame(channel = 1L, rr_ms = 800 + cumsum(rnorm(1000, sd = 2)))
#' traj <- ecgDFArolling(rr, window_beats = 200, step_beats = 50)
ecgDFArolling <- function(rr, window_beats = 300, step_beats = 30,
                          short_range = c(4, 16)) {
  stopifnot(is.data.frame(rr))
  if (!"rr_ms" %in% names(rr)) stop("Missing required column: rr_ms")
  if (!"channel" %in% names(rr)) rr$channel <- 1L
  stopifnot(window_beats >= 8, step_beats >= 1)

  channels <- unique(rr$channel)
  out <- list()
  for (ch in channels) {
    v <- rr$rr_ms[rr$channel == ch]
    n <- length(v)
    starts <- if (n < window_beats) 1L else
      seq(1L, n - window_beats + 1L, by = step_beats)
    for (s in starts) {
      e <- min(s + window_beats - 1L, n)
      out[[length(out) + 1]] <- data.frame(
        channel = ch, beat_start = s, beat_end = e,
        beat_center = (s + e) / 2,
        alpha1 = .dfa_alpha1(v[s:e], short_range),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

#' Detect DFA-a1 threshold crossings
#'
#' Finds where a rolling alpha1 trajectory crosses the two exercise-intensity
#' markers of Gronwald & Rogers: AT1 (aerobic threshold, alpha1 = 0.75) and
#' AT2 (anaerobic threshold, alpha1 = 0.5).
#'
#' @param alpha1 Numeric vector of rolling alpha1 values (e.g. the
#'   \code{alpha1} column from \code{\link{ecgDFArolling}}).
#' @param at1,at2 Threshold values (defaults 0.75 and 0.5).
#' @return A list with \code{at1}, \code{at2} and integer vectors
#'   \code{at1_crossings} / \code{at2_crossings} giving the indices immediately
#'   before each crossing (either direction).
#' @seealso \code{\link{ecgDFArolling}}
#' @export
#' @examples
#' a1 <- seq(1.2, 0.3, length.out = 20)   # descending intensity ramp
#' dfaThresholdCrossings(a1)
dfaThresholdCrossings <- function(alpha1, at1 = 0.75, at2 = 0.5) {
  a <- as.numeric(alpha1)
  # Keep original indexing (do not compress out NA windows) so crossing indices
  # map back to ecgDFArolling() rows; skip only pairs with a non-finite endpoint.
  cross <- function(thr) {
    if (length(a) < 2) return(integer(0))
    lhs <- a[-length(a)]
    rhs <- a[-1]
    ok <- is.finite(lhs) & is.finite(rhs)
    which(ok & ((lhs > thr & rhs <= thr) | (lhs < thr & rhs >= thr)))
  }
  list(at1 = at1, at2 = at2,
       at1_crossings = cross(at1), at2_crossings = cross(at2))
}
