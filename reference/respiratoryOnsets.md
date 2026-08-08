# Detect respiratory cycle onsets

Detects end-expiratory troughs in a zero-phase band-passed respiratory
signal. Consecutive troughs delimit complete breaths, and the
intervening maximum marks end inspiration.

## Usage

``` r
respiratoryOnsets(
  x,
  channel = NULL,
  band = c(0.1, 0.6),
  min_period_sec = 1,
  sr = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment`, numeric respiration vector, or
  `rip_calibration`.

- channel:

  Optional channel index or label.

- band:

  Two-element respiratory band in Hz.

- min_period_sec:

  Minimum breath period in seconds.

- sr:

  Sampling rate for a numeric vector.

- assay_name:

  Optional assay name for a `PhysioExperiment`.

## Value

A data frame with one row per complete breath, containing sample and
time coordinates for onset, peak, and offset plus inspiratory amplitude.

## See also

[`breathingRate()`](https://x-biosignal.github.io/PhysioECG/reference/breathingRate.md),
[`tidalMetrics()`](https://x-biosignal.github.io/PhysioECG/reference/tidalMetrics.md),
[`ripCalibrate()`](https://x-biosignal.github.io/PhysioECG/reference/ripCalibrate.md)

## Examples

``` r
sr <- 10
signal <- sin(2 * pi * 0.25 * (0:599) / sr)
head(respiratoryOnsets(signal, sr = sr))
#>   channel breath onset_sample peak_sample offset_sample onset_time peak_time
#> 1       1      1           31          51            71          3         5
#> 2       1      2           71          91           111          7         9
#> 3       1      3          111         131           151         11        13
#> 4       1      4          151         171           191         15        17
#> 5       1      5          191         211           231         19        21
#> 6       1      6          231         251           271         23        25
#>   offset_time insp_amplitude
#> 1           7              2
#> 2          11              2
#> 3          15              2
#> 4          19              2
#> 5          23              2
#> 6          27              2
```
