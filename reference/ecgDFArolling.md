# Rolling (time-resolved) DFA alpha1

Computes the short-range detrended-fluctuation scaling exponent (alpha1)
over a sliding window of beats, yielding a time-resolved alpha1
trajectory. This is the basis of the DFA-a1 exercise-intensity
thresholds of Gronwald & Rogers (2020).

## Usage

``` r
ecgDFArolling(rr, window_beats = 300, step_beats = 30, short_range = c(4, 16))
```

## Arguments

- rr:

  A data.frame with columns `channel` and `rr_ms`.

- window_beats:

  Window length in beats (default 300).

- step_beats:

  Step between successive windows in beats (default 30).

- short_range:

  Scale range (in beats) for alpha1 (default c(4, 16)).

## Value

A data.frame with one row per window and columns `channel`,
`beat_start`, `beat_end`, `beat_center` and `alpha1`.

## References

Peng, C.K. et al. (1995). Quantification of scaling exponents. *Chaos*,
5(1), 82-87. Gronwald, T. & Rogers, B. (2020). Fractal correlation
properties of heart rate variability as a biomarker. *Frontiers in
Physiology*, 11, 550572.

## See also

[`ecgDFA`](https://x-biosignal.github.io/PhysioECG/reference/ecgDFA.md),
[`dfaThresholdCrossings`](https://x-biosignal.github.io/PhysioECG/reference/dfaThresholdCrossings.md)

## Examples

``` r
set.seed(1)
rr <- data.frame(channel = 1L, rr_ms = 800 + cumsum(rnorm(1000, sd = 2)))
traj <- ecgDFArolling(rr, window_beats = 200, step_beats = 50)
```
