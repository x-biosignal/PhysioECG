# Heart-Rate Recovery (HRR)

Computes post-exercise heart-rate recovery: the fall in heart rate from
the exercise peak to fixed offsets into recovery (Cole et al. 1999; Imai
et al. 1994), plus the recovery slope. Accepts either a heart-rate (bpm)
or an RR (ms) series.

## Usage

``` r
ecgHRR(
  rr_or_hr,
  time_sec = NULL,
  peak_time = NULL,
  recovery_offsets_sec = c(60, 120),
  input = c("auto", "hr", "rr"),
  abnormal_bpm = 12
)
```

## Arguments

- rr_or_hr:

  A numeric vector of heart rate (bpm) or RR intervals (ms), or a
  data.frame with a `time_sec` column and an `hr` or `rr_ms` column.

- time_sec:

  Time of each sample in seconds (required for a vector input; taken
  from the `time_sec` column of a data.frame otherwise).

- peak_time:

  Time (seconds) of the exercise peak / start of recovery. If `NULL`
  (default), auto-detected as the time of maximum heart rate.

- recovery_offsets_sec:

  Offsets (seconds) after the peak at which to measure recovery (default
  `c(60, 120)`, i.e. HRR1 and HRR2).

- input:

  Interpretation of a vector input: "auto" (default; decided by
  magnitude – values above 250 are treated as RR in ms), "hr" or "rr".

- abnormal_bpm:

  HRR at 1 minute below which recovery is flagged abnormal (default 12
  bpm, per Cole et al. 1999).

## Value

A list with:

- peak_time, peak_hr:

  Recovery start time (s) and peak HR (bpm).

- hrr:

  A data.frame with columns `offset_sec`, `hr_at_offset` and `hrr_bpm`
  (peak HR minus HR at the offset).

- hrr1, hrr2:

  Convenience scalars for the first two offsets.

- slope_bpm_per_min:

  Slope of a linear fit of HR vs time over the recovery window (negative
  during recovery).

- hrr_1min:

  HRR at exactly 60 s (used for the abnormal flag).

- abnormal:

  TRUE if `hrr_1min < abnormal_bpm`.

## References

Cole, C.R. et al. (1999). Heart-rate recovery immediately after exercise
as a predictor of mortality. *New England Journal of Medicine*, 341(18),
1351-1357. Imai, K. et al. (1994). Vagally mediated heart rate recovery
after exercise is accelerated in athletes but blunted in patients with
chronic heart failure. *Journal of the American College of Cardiology*,
24(6), 1529-1535.

## See also

[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md),
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)

## Examples

``` r
# HR decays linearly from 170 to 120 bpm over 120 s of recovery
t <- seq(0, 120, by = 1)
hr <- 170 - (170 - 120) * t / 120
res <- ecgHRR(hr, time_sec = t, peak_time = 0)
res$hrr
#>   offset_sec hr_at_offset hrr_bpm
#> 1         60          145      25
#> 2        120          120      50
res$abnormal
#> [1] FALSE
```
