# HRV Frequency-Domain Analysis

Computes heart rate variability frequency-domain metrics from RR
interval data using Welch's method or Lomb-Scargle periodogram.

## Usage

``` r
ecgHRVfreq(
  rr,
  method = c("welch", "lomb", "ar"),
  vlf_band = c(0.003, 0.04),
  lf_band = c(0.04, 0.15),
  hf_band = c(0.15, 0.4),
  detrend = FALSE,
  detrend_lambda = 500,
  ar_order = NULL,
  rhythm_check = TRUE
)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

- method:

  Spectral estimation method: `"welch"` (default) for uniformly
  resampled FFT-based PSD, `"lomb"` for Lomb-Scargle periodogram on
  unevenly sampled data, or `"ar"` for a Burg autoregressive PSD.

- vlf_band:

  Numeric vector of length 2 defining VLF band in Hz (default: c(0.003,
  0.04)).

- lf_band:

  Numeric vector of length 2 defining LF band in Hz (default: c(0.04,
  0.15)).

- hf_band:

  Numeric vector of length 2 defining HF band in Hz (default: c(0.15,
  0.4)).

- detrend:

  Logical; if `TRUE`, apply smoothness-priors detrending
  ([`smoothnessPriorsDetrend`](https://x-biosignal.github.io/PhysioECG/reference/smoothnessPriorsDetrend.md))
  to the resampled tachogram before spectral estimation (Welch method).
  Default `FALSE` (mean removal).

- detrend_lambda:

  Smoothing parameter passed to
  [`smoothnessPriorsDetrend`](https://x-biosignal.github.io/PhysioECG/reference/smoothnessPriorsDetrend.md)
  when `detrend = TRUE` (default 500).

- ar_order:

  Integer AR model order for `method = "ar"`; `NULL` (default) selects
  the order automatically by AIC.

- rhythm_check:

  Logical; if `TRUE` (default), gate the result on
  [`ecgDetectAF`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectAF.md):
  metrics for channels classified as AF or frequent ectopy are set to
  `NA` and `rhythm`/`hrv_valid` columns are appended (HRV is undefined
  in atrial fibrillation).

## Value

A data.frame with one row per channel and the following columns:

- channel:

  Integer channel index.

- vlf:

  Absolute power in the very-low-frequency band (ms^2).

- lf:

  Absolute power in the low-frequency band (ms^2), associated with
  sympathetic and parasympathetic modulation.

- hf:

  Absolute power in the high-frequency band (ms^2), associated with
  parasympathetic (vagal) modulation.

- lf_hf_ratio:

  Ratio of LF to HF power, or `NA` if HF power is zero.

- total_power:

  Sum of VLF, LF, and HF power (ms^2).

- lf_nu, hf_nu:

  LF and HF power in normalized units, `100 * LF/(LF+HF)` and
  `100 * HF/(LF+HF)`; these sum to 100.

- lf_peak, hf_peak:

  Frequency (Hz) of the spectral peak within the LF and HF bands, or
  `NA` if the band is empty.

## References

Task Force of the European Society of Cardiology and the North American
Society of Pacing and Electrophysiology (1996). "Heart rate variability:
Standards of measurement, physiological interpretation and clinical
use." *Circulation*, 93(5), 1043–1065.

## See also

[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)
for computing RR intervals,
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)
for time-domain HRV metrics,
[`ecgHRVnonlinear`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVnonlinear.md)
for nonlinear HRV analysis,
[`ecgRRcorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRcorrect.md)
for ectopic beat correction before analysis.

## Examples

``` r
n <- 300
time_sec <- cumsum(rep(0.85, n))
rr_ms <- 850 + 30 * sin(2 * pi * 0.1 * time_sec)
rr <- data.frame(channel = rep(1L, n), rr_ms = rr_ms, time_sec = time_sec)
result <- ecgHRVfreq(rr, method = "welch")
```
