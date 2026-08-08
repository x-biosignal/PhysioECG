# Detrended Fluctuation Analysis of RR Intervals

Performs detrended fluctuation analysis (DFA) on RR interval data to
characterize fractal scaling properties. Computes alpha1 (short-range
correlations, 4–16 beats) and alpha2 (long-range correlations, 16–64
beats).

## Usage

``` r
ecgDFA(rr, short_range = c(4, 16), long_range = c(16, 64))
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

- short_range:

  Numeric vector of length 2 defining the scale range (in beats) for
  alpha1 (default: c(4, 16)).

- long_range:

  Numeric vector of length 2 defining the scale range (in beats) for
  alpha2 (default: c(16, 64)).

## Value

A data.frame with one row per channel and the following columns:

- channel:

  Integer channel index.

- alpha1:

  Short-range scaling exponent. Values near 1.0 indicate fractal-like
  (healthy) correlations; values near 0.5 indicate uncorrelated (random)
  behavior; values near 1.5 suggest Brownian noise. `NA` if the series
  is too short.

- alpha2:

  Long-range scaling exponent with the same interpretation as alpha1 but
  over larger time scales. `NA` if the series is too short.

## References

Peng, C.-K., et al. (1994). "Mosaic organization of DNA nucleotides."
*Physical Review E*, 49(2), 1685–1689.
[doi:10.1103/PhysRevE.49.1685](https://doi.org/10.1103/PhysRevE.49.1685)

## See also

[`ecgHRVnonlinear`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVnonlinear.md)
for the combined nonlinear analysis wrapper,
[`ecgSampleEntropy`](https://x-biosignal.github.io/PhysioECG/reference/ecgSampleEntropy.md)
for sample entropy,
[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md)
for Poincare plot descriptors.
