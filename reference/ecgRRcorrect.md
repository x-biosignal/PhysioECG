# Correct Ectopic Beats in RR Interval Data

Replaces or removes ectopic beats identified by
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md).
The `"interpolate"` method uses linear interpolation from surrounding
non-ectopic intervals, `"cubic_spline"` uses natural cubic-spline
interpolation (C1-smooth, following Kubios-style automatic correction),
and `"remove"` simply drops ectopic rows.

## Usage

``` r
ecgRRcorrect(rr, method = c("interpolate", "cubic_spline", "remove"))
```

## Arguments

- rr:

  A data.frame with an `is_ectopic` logical column, as returned by
  [`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md).

- method:

  Correction method: `"interpolate"` (default) replaces ectopic values
  with linear interpolation; `"cubic_spline"` uses natural cubic-spline
  interpolation (falls back to linear when fewer than four valid beats
  are available in a channel); `"remove"` drops ectopic rows.

## Value

A data.frame with corrected RR intervals. The `is_ectopic` column is
removed from the output. The result carries the ectopic burden as
attributes `n_ectopic` (count) and `pct_ectopic` (percent of beats
flagged ectopic).

## References

Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006). *Advanced Methods
and Tools for ECG Data Analysis*. Artech House.

## See also

[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md)
for detecting ectopic beats,
[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)
for computing RR intervals,
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)
for time-domain HRV analysis.
