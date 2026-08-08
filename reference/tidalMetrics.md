# Calculate tidal breathing metrics

Derives inspiratory and expiratory timing, duty cycle, tidal amplitude,
mean inspiratory flow, breathing rate, and minute ventilation.

## Usage

``` r
tidalMetrics(
  x,
  onsets = NULL,
  channel = NULL,
  calibration = NULL,
  band = c(0.1, 0.6),
  min_period_sec = 1,
  sr = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A respiration signal source accepted by
  [`respiratoryOnsets()`](https://x-biosignal.github.io/PhysioECG/reference/respiratoryOnsets.md).

- onsets:

  Optional precomputed table from
  [`respiratoryOnsets()`](https://x-biosignal.github.io/PhysioECG/reference/respiratoryOnsets.md).

- channel, band, min_period_sec, sr, assay_name:

  As in
  [`respiratoryOnsets()`](https://x-biosignal.github.io/PhysioECG/reference/respiratoryOnsets.md).

- calibration:

  Optional `rip_calibration` or positive numeric amplitude-to-liters
  scale.

## Value

A named list containing per-breath `breaths` and per-channel `summary`
data frames.

## Examples

``` r
sr <- 10
signal <- sin(2 * pi * 0.25 * (0:599) / sr)
tidalMetrics(signal, sr = sr)
#> $breaths
#>    channel breath Ti Te Ttot duty_cycle tidal_volume mean_insp_flow rate_bpm
#> 1        1      1  2  2    4        0.5            2              1       15
#> 2        1      2  2  2    4        0.5            2              1       15
#> 3        1      3  2  2    4        0.5            2              1       15
#> 4        1      4  2  2    4        0.5            2              1       15
#> 5        1      5  2  2    4        0.5            2              1       15
#> 6        1      6  2  2    4        0.5            2              1       15
#> 7        1      7  2  2    4        0.5            2              1       15
#> 8        1      8  2  2    4        0.5            2              1       15
#> 9        1      9  2  2    4        0.5            2              1       15
#> 10       1     10  2  2    4        0.5            2              1       15
#> 11       1     11  2  2    4        0.5            2              1       15
#> 12       1     12  2  2    4        0.5            2              1       15
#> 13       1     13  2  2    4        0.5            2              1       15
#> 14       1     14  2  2    4        0.5            2              1       15
#> 
#> $summary
#>   channel n_breaths rate_bpm tidal_volume_mean minute_ventilation
#> 1       1        14       15                 2                 30
#>   duty_cycle_mean
#> 1             0.5
#> 
```
