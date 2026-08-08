# Compute RR Intervals from Detected R-Peaks

Calculates the time intervals between consecutive R-peaks for each
channel. The resulting RR interval series is the standard input for all
HRV analysis functions in this package.

## Usage

``` r
ecgRRintervals(x, peaks)
```

## Arguments

- x:

  A PhysioExperiment object.

- peaks:

  A data.frame of detected peaks as returned by
  [`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md),
  with columns `channel`, `sample`, and `time_sec`.

## Value

A data.frame with one row per consecutive beat pair and the following
columns:

- channel:

  Integer channel index (1-based).

- rr_ms:

  RR interval in milliseconds (time between successive R-peaks).

- time_sec:

  Time of the first beat in each pair (seconds from signal onset).

Returns a zero-row data.frame with the same column structure if fewer
than two peaks are available.

## References

Pan, J. & Tompkins, W.J. (1985). "A real-time QRS detection algorithm."
*IEEE Transactions on Biomedical Engineering*, 32(3), 230–236.
[doi:10.1109/TBME.1985.325532](https://doi.org/10.1109/TBME.1985.325532)

## See also

[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
for R-peak detection,
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)
for time-domain HRV metrics,
[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md)
for frequency-domain HRV analysis,
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md)
for ectopic beat detection.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 60)
peaks <- ecgDetectRpeaks(pe)
rr <- ecgRRintervals(pe, peaks)
head(rr)
} # }
```
