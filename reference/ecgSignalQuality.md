# Assess ECG Signal Quality Per Channel

Computes per-channel signal quality metrics for ECG data stored in a
PhysioExperiment object. When detected R-peak locations are provided the
signal-to-noise ratio is estimated from QRS vs.\\ baseline power;
otherwise a variance-based estimate is used.

## Usage

``` r
ecgSignalQuality(x, peaks = NULL, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- peaks:

  Optional data.frame of detected R-peaks as returned by
  [`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md),
  with at least columns `channel` and `sample`.

- assay_name:

  Name of the assay to use. If `NULL` the default assay is used.

## Value

A data.frame with one row per channel and columns:

- channel:

  Integer channel index.

- snr_db:

  Signal-to-noise ratio in decibels.

- baseline_wander:

  RMS amplitude of the low-frequency drift.

- saturation_ratio:

  Fraction of samples within 1\\ min or max.

- quality_score:

  Composite quality score in the range \[0, 1\], where 1 indicates
  excellent quality.

## References

Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006). *Advanced Methods
and Tools for ECG Data Analysis*. Artech House.

Shaffer, F. & Ginsberg, J.P. (2017). "An overview of heart rate
variability metrics and norms." *Frontiers in Public Health*, 5, 258.
[doi:10.3389/fpubh.2017.00258](https://doi.org/10.3389/fpubh.2017.00258)

## See also

[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
for R-peak detection,
[`ecgBaselineCorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgBaselineCorrect.md)
for baseline wander correction,
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md)
for ectopic beat detection.
