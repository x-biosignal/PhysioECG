# Plot RR Tachogram

Plots the RR interval series against time, marking ectopic beats when an
`is_ectopic` column is present.

## Usage

``` r
plotTachogram(rr)
```

## Arguments

- rr:

  An rr data.frame (columns `rr_ms`, `time_sec` and optionally
  `channel`, `is_ectopic`) or an object returned by
  [`ecgProcess`](https://x-biosignal.github.io/PhysioECG/reference/ecgProcess.md).

## Value

A ggplot object.

## References

Task Force (1996). Heart rate variability. *Circulation*, 93(5),
1043-1065.

## See also

[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md),
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md)

## Examples

``` r
pe <- make_ecg(20000, sr = 500)
rr <- ecgQualityCheck(ecgRRintervals(pe, ecgDetectRpeaks(pe)))
plotTachogram(rr)
```
