# Detect Premature Ventricular Contractions (PVCs)

Flags ventricular ectopic beats from three classic features: an abnormal
QRS morphology (low correlation with an average-beat template),
prematurity of the preceding RR interval, and a full compensatory pause
(the premature beat plus its following pause together span about two
normal cycles). The compensatory pause is what distinguishes a PVC from
a non-compensatory supraventricular ectopic beat or sinus arrhythmia.

## Usage

``` r
ecgDetectPVC(
  x,
  peaks,
  rr,
  qrs_ms = 100,
  corr_threshold = 0.9,
  prematurity_max = 0.85,
  compensatory_min = 0.9,
  require_compensatory = FALSE,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- peaks:

  R-peak table from
  [`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
  (columns `channel`, `sample`).

- rr:

  RR-interval table from
  [`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)
  (columns `channel`, `rr_ms`).

- qrs_ms:

  QRS window half-width in ms for template matching (default 100).

- corr_threshold:

  Template-correlation below which a beat is morphologically abnormal
  (default 0.9).

- prematurity_max:

  Prematurity ratio (RR_pre / local median RR) below which a beat is
  premature (default 0.85).

- compensatory_min:

  Minimum (RR_pre + RR_post) / (2 \* local median RR) for a full
  compensatory pause (default 0.90).

- require_compensatory:

  If TRUE, also require a compensatory pause for a PVC verdict (default
  FALSE).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with one row per beat: `channel`, `beat`, `sample`,
`rr_pre_ms`, `rr_post_ms`, `prematurity`, `compensatory`, `morph_corr`
and the logical `is_pvc`.

## References

Petrenas, A., Marozas, V. & Sornmo, L. (2015). Low-complexity detection
of atrial fibrillation in continuous long-term monitoring. *Computers in
Biology and Medicine*, 65, 184-191.

## See also

[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md),
[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md),
[`ecgDetectAF`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectAF.md)

## Examples

``` r
pe <- make_ecg_pqrst(n_time = 5000, sr = 500)$pe
pk <- ecgDetectRpeaks(pe)
rr <- ecgRRintervals(pe, pk)
head(ecgDetectPVC(pe, pk, rr))
#>   channel beat sample rr_pre_ms rr_post_ms prematurity compensatory morph_corr
#> 1       1    1    567        NA        834          NA           NA  0.9992236
#> 2       1    2    984       834        834           1            1  0.9992244
#> 3       1    3   1401       834        834           1            1  0.9993503
#> 4       1    4   1818       834        834           1            1  0.9992215
#> 5       1    5   2235       834        834           1            1  0.9993522
#> 6       1    6   2652       834        834           1            1  0.9993327
#>   is_pvc
#> 1  FALSE
#> 2  FALSE
#> 3  FALSE
#> 4  FALSE
#> 5  FALSE
#> 6  FALSE
```
