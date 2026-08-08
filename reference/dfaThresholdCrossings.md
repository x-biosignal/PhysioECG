# Detect DFA-a1 threshold crossings

Finds where a rolling alpha1 trajectory crosses the two
exercise-intensity markers of Gronwald & Rogers: AT1 (aerobic threshold,
alpha1 = 0.75) and AT2 (anaerobic threshold, alpha1 = 0.5).

## Usage

``` r
dfaThresholdCrossings(alpha1, at1 = 0.75, at2 = 0.5)
```

## Arguments

- alpha1:

  Numeric vector of rolling alpha1 values (e.g. the `alpha1` column from
  [`ecgDFArolling`](https://x-biosignal.github.io/PhysioECG/reference/ecgDFArolling.md)).

- at1, at2:

  Threshold values (defaults 0.75 and 0.5).

## Value

A list with `at1`, `at2` and integer vectors `at1_crossings` /
`at2_crossings` giving the indices immediately before each crossing
(either direction).

## See also

[`ecgDFArolling`](https://x-biosignal.github.io/PhysioECG/reference/ecgDFArolling.md)

## Examples

``` r
a1 <- seq(1.2, 0.3, length.out = 20)   # descending intensity ramp
dfaThresholdCrossings(a1)
#> $at1
#> [1] 0.75
#> 
#> $at2
#> [1] 0.5
#> 
#> $at1_crossings
#> [1] 10
#> 
#> $at2_crossings
#> [1] 15
#> 
```
