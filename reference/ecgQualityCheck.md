# Detect Ectopic Beats in RR Interval Data

Identifies ectopic (abnormal) beats by comparing each RR interval to the
local median computed over a sliding window of 5 beats. Beats with
deviation exceeding the threshold are marked as ectopic.

## Usage

``` r
ecgQualityCheck(rr, threshold_ms = 300)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms`, and `time_sec`, as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

- threshold_ms:

  Maximum allowed deviation from local median in milliseconds (default:
  300).

## Value

The input data.frame with an additional logical column `is_ectopic`.

## References

Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006). *Advanced Methods
and Tools for ECG Data Analysis*. Artech House.

Task Force of the European Society of Cardiology and the North American
Society of Pacing and Electrophysiology (1996). "Heart rate variability:
Standards of measurement, physiological interpretation and clinical
use." *Circulation*, 93(5), 1043–1065.

## See also

[`ecgRRcorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRcorrect.md)
for correcting detected ectopic beats,
[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)
for computing RR intervals,
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)
for time-domain HRV analysis,
[`ecgSignalQuality`](https://x-biosignal.github.io/PhysioECG/reference/ecgSignalQuality.md)
for signal quality assessment.
