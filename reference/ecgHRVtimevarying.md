# Time-Varying HRV (Sliding-Window Time and Frequency Metrics)

Computes heart-rate-variability trajectories by sliding a time window
over the RR series and evaluating time-domain (SDNN, RMSSD) and
frequency-domain (LF, HF, LF/HF) metrics in each window, reusing
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)
and
[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md)
(Task Force 1996; Mainardi 2009).

## Usage

``` r
ecgHRVtimevarying(
  rr,
  window_sec = 300,
  step_sec = 30,
  freq_method = c("ar", "welch", "lomb"),
  detrend = FALSE,
  detrend_lambda = 500,
  min_beats = 20L,
  rhythm_check = FALSE
)
```

## Arguments

- rr:

  A data.frame with columns `channel`, `rr_ms` and `time_sec` (as
  returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)),
  assumed ordered in time per channel.

- window_sec:

  Analysis window length in seconds (default 300).

- step_sec:

  Step between successive windows in seconds (default 30).

- freq_method:

  Spectral method for
  [`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md):
  "ar" (default), "welch" or "lomb".

- detrend, detrend_lambda:

  Passed to
  [`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md)
  for optional smoothness-priors detrending of the resampled tachogram.

- min_beats:

  Minimum RR intervals in a window to compute metrics; windows with
  fewer beats yield `NA` metrics (default 20).

- rhythm_check:

  Passed to the per-window HRV functions (default FALSE, so the
  trajectory is continuous and not gated by the AF detector).

## Value

A data.frame with one row per window per channel and columns `channel`,
`window`, `time_start`, `time_center`, `n_beats`, `sdnn`, `rmssd`, `lf`,
`hf` and `lf_hf`.

## References

Task Force of the European Society of Cardiology and the North American
Society of Pacing and Electrophysiology (1996). Heart rate variability:
standards of measurement, physiological interpretation, and clinical
use. *Circulation*, 93(5), 1043-1065. Mainardi, L.T. (2009). On the
quantification of heart rate variability spectral parameters using
time-frequency and time-varying methods. *Philosophical Transactions of
the Royal Society A*, 367(1887), 255-275.

## See also

[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md),
[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md),
[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)

## Examples

``` r
set.seed(1)
n <- 600
rr <- data.frame(channel = 1L,
                 rr_ms = 800 + 25 * sin(2 * pi * 0.1 * cumsum(rep(0.8, n))) +
                   rnorm(n, sd = 10))
rr$time_sec <- cumsum(rr$rr_ms) / 1000
traj <- ecgHRVtimevarying(rr, window_sec = 120, step_sec = 30)
head(traj)
#>   channel window  time_start time_center n_beats     sdnn    rmssd       lf
#> 1       1      1   0.8057793    60.80578     150 19.04766 15.23051 286.0625
#> 2       1      2  30.8057793    90.80578     150 19.60577 15.93059 296.2325
#> 3       1      3  60.8057793   120.80578     150 20.10344 16.30058 326.4678
#> 4       1      4  90.8057793   150.80578     150 20.81888 17.25596 314.6790
#> 5       1      5 120.8057793   180.80578     150 21.21406 17.74164 342.8496
#> 6       1      6 150.8057793   210.80578     150 20.86647 17.52086 318.8135
#>         hf    lf_hf
#> 1 20.42544 14.00521
#> 2 22.10450 13.40146
#> 3 24.45791 13.34815
#> 4 27.57396 11.41218
#> 5 27.30200 12.55768
#> 6 25.95065 12.28538
```
