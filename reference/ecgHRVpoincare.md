# Poincare Plot Descriptors for HRV Analysis

Computes Poincare plot descriptors (SD1, SD2, SD1/SD2 ratio) from RR
interval data. SD1 reflects short-term variability (perpendicular to the
identity line), while SD2 reflects long-term variability (along the
identity line).

## Usage

``` r
ecgHRVpoincare(rr)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

## Value

A data.frame with one row per channel and the following columns:

- channel:

  Integer channel index.

- sd1:

  Standard deviation perpendicular to the identity line (ms), reflecting
  beat-to-beat (short-term) variability.

- sd2:

  Standard deviation along the identity line (ms), reflecting long-term
  variability.

- sd1_sd2_ratio:

  Ratio of SD1 to SD2, or `NA` if SD2 is zero.

## References

Task Force of the European Society of Cardiology and the North American
Society of Pacing and Electrophysiology (1996). "Heart rate variability:
Standards of measurement, physiological interpretation and clinical
use." *Circulation*, 93(5), 1043–1065.

## See also

[`ecgHRVnonlinear`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVnonlinear.md)
for the combined nonlinear analysis wrapper,
[`ecgSampleEntropy`](https://x-biosignal.github.io/PhysioECG/reference/ecgSampleEntropy.md)
for sample entropy,
[`ecgDFA`](https://x-biosignal.github.io/PhysioECG/reference/ecgDFA.md)
for detrended fluctuation analysis.
