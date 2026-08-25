# One-Call ECG Processing Pipeline

Runs the standard ECG preprocessing chain end to end – baseline
correction, R-peak detection, RR extraction, ectopic-beat correction and
(optionally) a bSQI signal-quality gate – and stores the peaks, RR
series and provenance in the object's `metadata()` (after the neurokit2
`ecg_process()` pattern).

## Usage

``` r
ecgProcess(
  x,
  detector = "pan_tompkins",
  correct = "cubic_spline",
  sqi_gate = TRUE,
  baseline = TRUE,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- detector:

  R-peak detector passed to
  [`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
  (default "pan_tompkins").

- correct:

  Ectopic-beat correction passed to
  [`ecgRRcorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRcorrect.md):
  "cubic_spline" (default), "interpolate" or "remove".

- sqi_gate:

  If TRUE (default), also compute the bSQI quality gate
  ([`ecgBeatSQIgate`](https://x-biosignal.github.io/PhysioECG/reference/ecgBeatSQIgate.md)).

- baseline:

  If TRUE (default), baseline-correct with
  [`ecgBaselineCorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgBaselineCorrect.md)
  and detect on the corrected assay.

- assay_name:

  Input assay name (default: first assay).

## Value

The (baseline-corrected) PhysioExperiment with an `ecg` entry in
`metadata()` holding `peaks`, `rr` (with ectopic flags), `rr_corrected`,
`bsqi` and the detect-assay name, plus a provenance record appended via
[`withProvenance`](https://x-biosignal.github.io/PhysioCore//reference/withProvenance.html).

## See also

[`ecgAnalyze`](https://x-biosignal.github.io/PhysioECG/reference/ecgAnalyze.md),
[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md),
[`ecgRRcorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRcorrect.md),
[`ecgBeatSQIgate`](https://x-biosignal.github.io/PhysioECG/reference/ecgBeatSQIgate.md)

## Examples

``` r
pe <- make_ecg(n_time = 20000, sr = 500)
proc <- ecgProcess(pe)
names(S4Vectors::metadata(proc)$ecg)
#> [1] "detect_assay" "detector"     "correct"      "peaks"        "rr"          
#> [6] "rr_corrected" "bsqi"        
```
