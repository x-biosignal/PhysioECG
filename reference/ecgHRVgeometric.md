# Geometric HRV indices (triangular index and TINN)

Computes the geometric heart-rate-variability measures defined by the
Task Force (1996) from a series of NN (normal-to-normal) intervals: the
HRV triangular index and the triangular interpolation of the NN interval
histogram (TINN).

## Usage

``` r
ecgHRVgeometric(rr, bin_ms = 7.8125)
```

## Arguments

- rr:

  Numeric vector of NN/RR intervals in milliseconds, or a data frame
  with an `rr_ms` column (as produced by the RR-interval readers); a
  data frame is reduced to its `rr_ms` column.

- bin_ms:

  Histogram bin width in milliseconds. Defaults to `7.8125` ms (`1/128`
  s), the Task Force standard.

## Value

A list with:

- HRV_triangular_index:

  Total number of NN intervals divided by the height (modal count) of
  the NN histogram.

- TINN:

  Baseline width, in milliseconds, of the triangle that best
  approximates the NN histogram in a least-squares sense.

## References

Task Force of the ESC and NASPE (1996). Heart rate variability:
standards of measurement. *Circulation*, 93(5), 1043-1065.

## See also

[`ecgHRVsegmented`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVsegmented.md),
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)

## Examples

``` r
set.seed(1)
rr <- 800 + rnorm(500, sd = 40)  # NN intervals in ms
ecgHRVgeometric(rr)
#> $HRV_triangular_index
#> [1] 9.803922
#> 
#> $TINN
#> [1] 156.25
#> 
```
