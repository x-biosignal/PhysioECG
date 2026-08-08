# Sample Entropy of RR Intervals

Computes sample entropy (SampEn) from RR interval data. Sample entropy
measures the regularity or predictability of a time series. Lower values
indicate more regular (predictable) signals, while higher values
indicate more complex (irregular) signals.

## Usage

``` r
ecgSampleEntropy(rr, m = 2L, r_factor = 0.2)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

- m:

  Embedding dimension (default: 2). Length of template patterns to
  compare.

- r_factor:

  Tolerance factor (default: 0.2). The tolerance `r` is computed as
  `r_factor * sd(rr_ms)`.

## Value

A data.frame with one row per channel and the following columns:

- channel:

  Integer channel index.

- sample_entropy:

  Sample entropy value (nats). Lower values indicate more regular
  signals; higher values indicate more complex signals. `NA` if the
  series is too short or constant.

- m:

  Embedding dimension used.

- r:

  Tolerance threshold (ms) computed as `r_factor * sd(rr_ms)`.

## References

Richman, J.S. & Moorman, J.R. (2000). "Physiological time-series
analysis using approximate entropy and sample entropy." *American
Journal of Physiology-Heart and Circulatory Physiology*, 278(6),
H2039–H2049.
[doi:10.1152/ajpheart.2000.278.6.H2039](https://doi.org/10.1152/ajpheart.2000.278.6.H2039)

## See also

[`ecgHRVnonlinear`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVnonlinear.md)
for the combined nonlinear analysis wrapper,
[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md)
for Poincare plot descriptors,
[`ecgDFA`](https://x-biosignal.github.io/PhysioECG/reference/ecgDFA.md)
for detrended fluctuation analysis.
