# Plot Poincare (Lorenz) Recurrence Map

Scatters each RR interval against the next (\\RR_n\\ vs \\RR\_{n+1}\\)
and overlays the Poincare SD1/SD2 dispersion ellipse from
[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md).
The SD1 and SD2 used for the ellipse are attached as the `"sd1"`/`"sd2"`
attributes of the returned plot.

## Usage

``` r
plotPoincare(rr)
```

## Arguments

- rr:

  An rr data.frame or an object from
  [`ecgProcess`](https://x-biosignal.github.io/PhysioECG/reference/ecgProcess.md).

## Value

A ggplot object with `"sd1"`/`"sd2"` attributes (per channel).

## References

Brennan, M., Palaniswami, M. & Kamen, P. (2001). Do existing measures of
Poincare plot geometry reflect nonlinear features of heart rate
variability? *IEEE TBME*, 48(11), 1342-1347.

## See also

[`ecgHRVpoincare`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVpoincare.md),
[`ecgHRVnonlinear`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVnonlinear.md)

## Examples

``` r
pe <- make_ecg(30000, sr = 500)
rr <- ecgRRintervals(pe, ecgDetectRpeaks(pe))
plotPoincare(rr)
```
