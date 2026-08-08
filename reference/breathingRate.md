# Estimate breathing rate

Estimates rate from complete inter-trough breath periods or from the
dominant respiratory-band spectral peak.

## Usage

``` r
breathingRate(
  x,
  channel = NULL,
  method = c("peaks", "spectral"),
  band = c(0.1, 0.6),
  min_period_sec = 1,
  sr = NULL,
  assay_name = NULL
)
```

## Arguments

- x, channel, band, min_period_sec, sr, assay_name:

  As in
  [`respiratoryOnsets()`](https://x-biosignal.github.io/PhysioECG/reference/respiratoryOnsets.md).

- method:

  `"peaks"` for breath intervals or `"spectral"` for the dominant
  respiratory frequency.

## Value

A data frame with one row per selected channel and columns `channel`,
`method`, `rate_bpm`, `n_breaths`, and `sd_bpm`.

## Examples

``` r
sr <- 10
signal <- sin(2 * pi * 0.25 * (0:599) / sr)
breathingRate(signal, sr = sr)
#>   channel method rate_bpm n_breaths sd_bpm
#> 1       1  peaks       15        14      0
```
