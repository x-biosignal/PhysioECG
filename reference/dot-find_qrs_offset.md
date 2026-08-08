# Find QRS offset by searching forward from S-wave trough

Strategy: first locate the S-wave trough (local minimum after R), then
estimate baseline from far after S, and walk forward from S until the
signal returns to near baseline.

## Usage

``` r
.find_qrs_offset(sig, r_pos, sr, n_time)
```
