# HRV Time-Domain Metrics

Compute standard heart rate variability (HRV) time-domain metrics from
RR interval data. Calculates SDNN, RMSSD, pNN50, mean RR interval, and
mean heart rate for each channel.

## Usage

``` r
ecgHRVtime(rr, rhythm_check = TRUE)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

- rhythm_check:

  Logical; if `TRUE` (default), gate the result on
  [`ecgDetectAF`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectAF.md):
  metrics for channels classified as AF or frequent ectopy are set to
  `NA` and `rhythm`/`hrv_valid` columns are appended (HRV is undefined
  in atrial fibrillation).

## Value

A data.frame with one row per channel and the following columns:

- channel:

  Integer channel index.

- mean_rr:

  Mean RR interval in milliseconds.

- sdnn:

  Standard deviation of all RR intervals (ms), reflecting overall HRV.

- rmssd:

  Root mean square of successive RR interval differences (ms),
  reflecting short-term vagal modulation.

- pnn50:

  Percentage of successive RR intervals differing by more than 50 ms.

- mean_hr:

  Mean heart rate in beats per minute (60000 / mean_rr).

## References

Task Force of the European Society of Cardiology and the North American
Society of Pacing and Electrophysiology (1996). "Heart rate variability:
Standards of measurement, physiological interpretation and clinical
use." *Circulation*, 93(5), 1043–1065.

## See also

[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)
for computing RR intervals,
[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md)
for frequency-domain HRV analysis,
[`ecgHRVnonlinear`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVnonlinear.md)
for nonlinear HRV metrics,
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md)
for ectopic beat detection before analysis.
