# Find T-wave end by searching forward from T-peak for baseline return

Uses a tangent-intercept method: draws a tangent at the steepest
downslope of the T wave and finds where it crosses the baseline level.
Falls back to amplitude threshold if tangent method fails.

## Usage

``` r
.find_t_end(sig, t_peak, sr, n_time)
```

## Arguments

- sig:

  Signal vector.

- t_peak:

  Sample index of the T-wave peak.

- sr:

  Sampling rate in Hz.

- n_time:

  Total number of samples.

## Value

Sample index of T-wave end, or NA_integer\_.
