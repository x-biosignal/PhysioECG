# Find QRS onset by searching backward from Q-wave trough

Strategy: first locate the Q-wave trough (local minimum before R), then
estimate baseline from far before Q, and walk backward from Q until the
signal amplitude is within a small fraction of the Q deflection from
baseline.

## Usage

``` r
.find_qrs_onset(sig, r_pos, sr, n_time)
```
