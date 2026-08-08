# Smoothness-priors detrending

Removes slow trends from a uniformly sampled series (e.g. an
interpolated RR tachogram) using the smoothness-priors method of
Tarvainen et al. (2002). The trend estimate is \\(I + \lambda^2 D_2^\top
D_2)^{-1} z\\ where \\D_2\\ is the second-order difference operator, and
the detrended (stationary) component is \\z\_{stat} = (I - (I +
\lambda^2 D_2^\top D_2)^{-1}) z\\. Larger `lambda` lowers the cutoff
frequency, removing more low-frequency content.

## Usage

``` r
smoothnessPriorsDetrend(series, lambda = 500, fs = 4)
```

## Arguments

- series:

  Numeric vector, uniformly sampled.

- lambda:

  Smoothing parameter (regularisation). Higher values remove slower
  trends (lower cutoff). Default 500.

- fs:

  Sampling rate of `series` in Hz, used only to report the approximate
  half-power cutoff frequency. Default 4.

## Value

The detrended series (numeric vector, near-zero mean), carrying
attributes `lambda` and `cutoff_hz` (the approximate half-power cutoff
of the implied high-pass filter).

## References

Tarvainen, M.P., Ranta-aho, P.O. & Karjalainen, P.A. (2002). An advanced
detrending method with application to HRV analysis. *IEEE Transactions
on Biomedical Engineering*, 49(2), 172-175.

## See also

[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md)

## Examples

``` r
t <- seq(0, 60, by = 0.25)
x <- 0.01 * t + sin(2 * pi * 0.25 * t)   # slow trend + HF oscillation
xd <- smoothnessPriorsDetrend(x, lambda = 500, fs = 4)
abs(mean(xd)) < 1e-6
#> [1] TRUE
```
