# Segmented (long-term) HRV indices: SDANN and SDNN index

Splits an NN-interval series into consecutive fixed-duration segments (5
minutes by default) and returns the two Task Force (1996) long-term
variability measures.

## Usage

``` r
ecgHRVsegmented(rr, segment_sec = 300)
```

## Arguments

- rr:

  Numeric vector of NN/RR intervals in milliseconds.

- segment_sec:

  Segment length in seconds (default 300 = 5 minutes).

## Value

A list with:

- SDANN:

  Standard deviation of the per-segment mean NN intervals (ms) -
  long-term variability.

- SDNN_index:

  Mean of the per-segment SDNN values (ms) - average short-term
  variability.

- n_segments:

  Number of segments used.

## References

Task Force of the ESC and NASPE (1996). *Circulation*, 93(5), 1043-1065.

## See also

[`ecgHRVgeometric`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVgeometric.md),
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)

## Examples

``` r
set.seed(1)
rr <- 800 + rnorm(3000, sd = 40)
ecgHRVsegmented(rr, segment_sec = 300)
#> $SDANN
#> [1] 1.7121
#> 
#> $SDNN_index
#> [1] 41.39023
#> 
#> $n_segments
#> [1] 8
#> 
```
