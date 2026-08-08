# Correct Baseline Wander in ECG Signals

Removes low-frequency baseline drift from ECG data using either a
high-pass moving-average subtraction or a running-median subtraction.

## Usage

``` r
ecgBaselineCorrect(
  x,
  method = c("highpass", "median"),
  cutoff = 0.5,
  assay_name = NULL,
  output_assay = "baseline_corrected"
)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- method:

  Correction method: `"highpass"` (default) subtracts a moving average;
  `"median"` subtracts a running median.

- cutoff:

  Approximate cutoff frequency in Hz (default: 0.5). Controls the window
  size of the moving average or median filter.

- assay_name:

  Name of the input assay. If `NULL` the default assay is used.

- output_assay:

  Name of the assay in which to store the corrected signal (default:
  `"baseline_corrected"`).

## Value

A PhysioExperiment with the corrected signal stored in `output_assay`.

## References

Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006). *Advanced Methods
and Tools for ECG Data Analysis*. Artech House.

## See also

[`ecgSignalQuality`](https://x-biosignal.github.io/PhysioECG/reference/ecgSignalQuality.md)
for signal quality assessment,
[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
for R-peak detection,
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md)
for ectopic beat detection.
