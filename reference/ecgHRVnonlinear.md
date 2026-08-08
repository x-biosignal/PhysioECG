# Nonlinear HRV Analysis (Convenience Wrapper)

Computes all nonlinear HRV metrics by calling
[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md),
[`ecgSampleEntropy`](https://x-biosignal.github.io/PhysioECG/reference/ecgSampleEntropy.md),
and
[`ecgDFA`](https://x-biosignal.github.io/PhysioECG/reference/ecgDFA.md),
and merging the results into a single data.frame.

## Usage

``` r
ecgHRVnonlinear(
  rr,
  m = 2L,
  r_factor = 0.2,
  short_range = c(4, 16),
  long_range = c(16, 64),
  rhythm_check = TRUE
)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

- m:

  Embedding dimension for sample entropy (default: 2).

- r_factor:

  Tolerance factor for sample entropy (default: 0.2).

- short_range:

  Scale range for DFA alpha1 (default: c(4, 16)).

- long_range:

  Scale range for DFA alpha2 (default: c(16, 64)).

- rhythm_check:

  Logical; if `TRUE` (default), gate the result on
  [`ecgDetectAF`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectAF.md):
  metrics for channels classified as AF or frequent ectopy are set to
  `NA` and `rhythm`/`hrv_valid` columns are appended (HRV is undefined
  in atrial fibrillation).

## Value

A data.frame with columns: channel, sd1, sd2, sd1_sd2_ratio,
sample_entropy, m, r, alpha1, alpha2.

## References

Shaffer, F. & Ginsberg, J.P. (2017). "An overview of heart rate
variability metrics and norms." *Frontiers in Public Health*, 5, 258.
[doi:10.3389/fpubh.2017.00258](https://doi.org/10.3389/fpubh.2017.00258)

Task Force of the European Society of Cardiology and the North American
Society of Pacing and Electrophysiology (1996). "Heart rate variability:
Standards of measurement, physiological interpretation and clinical
use." *Circulation*, 93(5), 1043–1065.

## See also

[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md)
for Poincare plot descriptors,
[`ecgSampleEntropy`](https://x-biosignal.github.io/PhysioECG/reference/ecgSampleEntropy.md)
for sample entropy,
[`ecgDFA`](https://x-biosignal.github.io/PhysioECG/reference/ecgDFA.md)
for detrended fluctuation analysis,
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)
for time-domain HRV metrics,
[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md)
for frequency-domain HRV analysis.
