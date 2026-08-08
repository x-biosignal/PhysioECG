# Calibrate two-belt respiratory inductance plethysmography

Fits the Konno-Mead two-degree-of-freedom volume model from ribcage and
abdominal belt signals. Least-squares calibration uses a reference
volume trace; qualitative diagnostic calibration (QDC) equalizes the
variability of the two belt contributions.

## Usage

``` r
ripCalibrate(
  rc,
  ab,
  reference = NULL,
  sr = NULL,
  method = c("qdc", "lsq"),
  band = c(0.1, 0.6),
  window = NULL,
  fit_intercept = FALSE,
  known_volume = NULL
)
```

## Arguments

- rc:

  Ribcage signal as a numeric vector, or a `PhysioExperiment`.

- ab:

  Abdominal numeric vector, or a channel index/label when `rc` is a
  `PhysioExperiment`.

- reference:

  Optional reference volume signal, required for `"lsq"`.

- sr:

  Sampling rate in Hz for numeric signals. A `PhysioExperiment` supplies
  its own sampling rate.

- method:

  Calibration method, `"qdc"` or `"lsq"`.

- band:

  Two-element respiratory band in Hz.

- window:

  Optional inclusive calibration sample range `c(lo, hi)`.

- fit_intercept:

  Logical; include an intercept in least squares.

- known_volume:

  Optional known inspiratory volume used for absolute scaling.

## Value

A `rip_calibration` object containing coefficients, calibration
diagnostics, and the calibrated volume-proportional signal.

## References

Konno K, Mead J (1967). Measurement of the separate volume changes of
rib cage and abdomen during breathing. *Journal of Applied Physiology*,
22:407-422.

Sackner MA, Watson H, Belsito AS, et al. (1989). Calibration of
respiratory inductive plethysmograph during natural breathing. *Journal
of Applied Physiology*, 66:410-420.
[doi:10.1152/jappl.1989.66.1.410](https://doi.org/10.1152/jappl.1989.66.1.410)

## Examples

``` r
sr <- 10
time <- (0:599) / sr
rc <- 2 * sin(2 * pi * 0.25 * time)
ab <- sin(2 * pi * 0.25 * time)
ripCalibrate(rc, ab, sr = sr)
#> <rip_calibration> method=qdc, sr=10 Hz
#>   coefficients: RC=1, AB=2  K=2  scale=1
```
