# Beat-level Signal Quality Index (bSQI)

Quantifies ECG quality per channel by the agreement between two
independent QRS detectors. A clean signal makes both detectors agree on
nearly every beat (bSQI near 1); noise makes them disagree (bSQI drops).

## Usage

``` r
ecgBeatSQI(
  x,
  detector_a = "pan_tompkins",
  detector_b = "elgendi",
  match_window_ms = 150,
  threshold_factor = 0.6,
  refractory_ms = 200,
  beta = 0.08,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- detector_a, detector_b:

  Detector methods to compare (see
  [`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md));
  defaults `"pan_tompkins"` and `"elgendi"`.

- match_window_ms:

  Two beats match when within this window (default 150).

- threshold_factor, refractory_ms, beta:

  Passed to the detectors.

- assay_name:

  Assay to use (default: the default assay).

## Value

A data.frame with one row per channel: `n_a`, `n_b` (beat counts),
`n_matched`, and `bSQI = matched / (n_a + n_b - matched)` (Jaccard
agreement in `[0, 1]`; `NA` for a silent channel).

## References

Li, Q., Mark, R.G. & Clifford, G.D. (2008). Robust heart rate estimation
from multiple asynchronous noisy sources. *Physiological Measurement*,
29(1), 15-32.

## See also

[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md),
[`ecgBeatSQIgate`](https://x-biosignal.github.io/PhysioECG/reference/ecgBeatSQIgate.md)

## Examples

``` r
pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
ecgBeatSQI(pe)
#>   channel n_a n_b n_matched bSQI
#> 1       1  11  11        11    1
```
