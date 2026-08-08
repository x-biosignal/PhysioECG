# Ectopic Beat Burden

Convenience summary of the ectopic-beat burden in an RR-interval
data.frame that has been passed through
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md).

## Usage

``` r
ecgEctopicBurden(rr)
```

## Arguments

- rr:

  A data.frame with an `is_ectopic` logical column.

## Value

A list with `n_ectopic` (count of ectopic beats), `n_total` (total
beats) and `pct_ectopic` (percent ectopic).

## See also

[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md),
[`ecgRRcorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRcorrect.md)

## Examples

``` r
rr <- data.frame(channel = 1L, rr_ms = c(800, 810, 400, 805, 815),
                 is_ectopic = c(FALSE, FALSE, TRUE, FALSE, FALSE))
ecgEctopicBurden(rr)
#> $n_ectopic
#> [1] 1
#> 
#> $n_total
#> [1] 5
#> 
#> $pct_ectopic
#> [1] 20
#> 
```
