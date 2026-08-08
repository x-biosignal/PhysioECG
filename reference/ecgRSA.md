# Respiratory Sinus Arrhythmia (RSA)

Quantifies respiratory sinus arrhythmia – the respiration-linked
oscillation of the RR interval – by the peak-valley method (Grossman) or
the Porges-Bohrer method (band-limited variance of the tachogram).

## Usage

``` r
ecgRSA(
  rr,
  resp = NULL,
  method = c("peak_valley", "porges_bohrer"),
  fs = 4,
  resp_band = c(0.12, 0.4)
)
```

## Arguments

- rr:

  A data.frame with columns `rr_ms` and `time_sec` (and optionally
  `channel`), or a numeric vector of RR intervals (ms), in which case
  beat times are taken as their cumulative sum.

- resp:

  Optional respiration surrogate (unused by the current methods,
  reserved for respiration-gated peak-valley).

- method:

  "peak_valley" (mean band-filtered peak-to-valley RR amplitude,
  default) or "porges_bohrer" (natural log of the band-filtered RR
  variance).

- fs:

  Resampling rate (Hz) of the tachogram (default 4).

- resp_band:

  Respiratory band (Hz) for filtering (default `c(0.12, 0.4)`).

## Value

A list with `method`, `rsa` (the RSA magnitude: mean peak-to-valley RR
in ms, or `ln(variance)` for Porges-Bohrer), and, for "peak_valley",
`n_cycles` (number of respiratory half-cycles used).

## References

Grossman, P., van Beek, J. & Wientjes, C. (1990). A comparison of three
quantification methods for estimation of respiratory sinus arrhythmia.
*Psychophysiology*, 27(6), 702-714. Lewis, G.F. et al. (2012).
Statistical strategies to quantify RSA: are commonly used metrics
equivalent? *Biological Psychology*, 89(2), 349-364.

## See also

[`ecgEDR`](https://x-biosignal.github.io/PhysioECG/reference/ecgEDR.md),
[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)

## Examples

``` r
set.seed(1)
t <- cumsum(rep(0.8, 300))
rr <- data.frame(rr_ms = 800 + 40 * sin(2 * pi * 0.25 * t) + rnorm(300, 5),
                 time_sec = t)
ecgRSA(rr)$rsa
#> [1] 70.13977
```
