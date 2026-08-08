# Detect Atrial Fibrillation / Frequent Ectopy from RR Irregularity

Classifies the cardiac rhythm of each channel as regular sinus, atrial
fibrillation (AF), or frequent ectopy from the irregularity of the RR
interval series. HRV metrics are physiologically undefined in AF, so
this gate is used to invalidate them (see `rhythm_check` in
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md),
[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md)
and
[`ecgHRVnonlinear`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVnonlinear.md)).

## Usage

``` r
ecgDetectAF(
  rr,
  window_beats = NULL,
  ectopic_frac_max = 0.05,
  rmssd_ratio_af = 0.15,
  entropy_af = 0.7,
  n_bins = 16L,
  drr_range = c(-600, 600),
  threshold_ms = 300,
  min_beats = 20L
)
```

## Arguments

- rr:

  A data.frame with columns `channel` and `rr_ms`, as returned by
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md).

- window_beats:

  If non-NULL, an integer window length (in beats): instead of a single
  per-channel verdict, a continuous AF probability is returned for each
  non-overlapping window of `window_beats` beats (see Value).

- ectopic_frac_max:

  Maximum ectopic fraction for a `"sinus"` verdict (default 0.05).

- rmssd_ratio_af:

  Minimum RMSSD/mean-RR ratio consistent with AF (default 0.15).

- entropy_af:

  Minimum normalized dRR Shannon entropy consistent with AF (default
  0.7, range 0-1).

- n_bins:

  Number of histogram bins for the dRR entropy (default 16).

- drr_range:

  Numeric length-2 dRR window (ms) for the entropy histogram (default
  c(-600, 600)).

- threshold_ms:

  Deviation threshold (ms) passed to
  [`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md)
  for the ectopic fraction (default 300).

- min_beats:

  Minimum RR intervals required to classify; fewer yields `NA` (default
  20).

## Value

With `window_beats = NULL` (default), a data.frame with one row per
channel and columns `channel`, `rhythm` (`"sinus"`, `"AF"` or
`"frequent_ectopy"`), `ectopic_frac`, `rmssd_ratio`, `drr_entropy`,
`af_prob` (a continuous AF probability in 0-1) and `n_beats`. With
`window_beats` set, a data.frame with one row per non-overlapping window
and columns `channel`, `window`, `beat_start`, `n_beats`, `af_prob` and
its components (`rmssd_ratio`, `drr_entropy`, `sample_entropy`,
`sd1_sd2`).

## Details

Three descriptors are combined: the ectopic-beat fraction from
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md),
the RMSSD/mean-RR normalized irregularity index, and the normalized
Shannon entropy of the successive-difference (dRR) histogram. AF
requires *both* large beat-to-beat variability (`rmssd_ratio`) and a
broad, disordered dRR distribution (`drr_entropy`); a large but
concentrated irregularity (e.g. an ectopic couplet) falls through to
`"frequent_ectopy"`.

## References

Tateno, K. & Glass, L. (2001). Automatic detection of atrial
fibrillation using the coefficient of variation and density histograms
of RR and dRR intervals. *Medical & Biological Engineering & Computing*,
39(6), 664-671. Sarkar, S., Ritscher, D. & Mehra, R. (2008). A detector
for a chronic implantable atrial tachyarrhythmia monitor. *IEEE
Transactions on Biomedical Engineering*, 55(3), 1219-1224.

## See also

[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md),
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md)

## Examples

``` r
set.seed(1)
rr <- data.frame(channel = 1L, rr_ms = 800 + rnorm(100, sd = 20),
                 time_sec = cumsum(rep(0.8, 100)))
ecgDetectAF(rr)
#>   channel rhythm ectopic_frac rmssd_ratio drr_entropy   af_prob n_beats
#> 1       1  sinus            0  0.03164052   0.2485076 0.5644174     100
```
