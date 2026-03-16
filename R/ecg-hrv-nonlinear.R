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

    # Generate scales: integer values from 4 up to n/4
    max_scale <- max(4L, as.integer(n %/% 4))
    scales <- unique(as.integer(seq(4, max_scale, length.out = min(30L, max_scale - 3L))))
    scales <- scales[scales >= 4]

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
                            long_range = c(16, 64)) {
  poincare <- ecgHRVpoincare(rr)
  sampen <- ecgSampleEntropy(rr, m = m, r_factor = r_factor)
  dfa <- ecgDFA(rr, short_range = short_range, long_range = long_range)

  result <- merge(poincare, sampen, by = "channel")
  result <- merge(result, dfa, by = "channel")
  result
}
