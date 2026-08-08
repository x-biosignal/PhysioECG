# Delineate ECG Waveform Morphology

For each detected R-peak, identifies QRS complex boundaries (onset and
offset) and P-wave and T-wave peaks. QRS boundaries are detected using a
gradient-based search from the R-peak, while P and T waves are found as
local maxima in physiologically plausible time windows.

## Usage

``` r
ecgDelineate(x, peaks, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- peaks:

  A data.frame of detected R-peaks as returned by
  [`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md),
  with columns `channel` and `sample`.

- assay_name:

  Name of the assay to use. If `NULL`, the default assay is used.

## Value

A data.frame with one row per beat and the following columns:

- channel:

  Integer channel index (1-based).

- beat:

  Integer beat number within the channel (1-based).

- r_peak:

  Sample index of the R-peak.

- qrs_onset:

  Sample index of QRS complex onset.

- qrs_offset:

  Sample index of QRS complex offset (J-point).

- qrs_duration_ms:

  QRS complex duration in milliseconds.

- p_peak:

  Sample index of P-wave peak, or `NA` if not found in the search window
  (300–80 ms before R-peak).

- t_peak:

  Sample index of T-wave peak, or `NA` if not found in the search window
  (80–500 ms after R-peak).

- t_end:

  Sample index of T-wave end estimated by tangent-intercept method, or
  `NA` if T-wave not found.

Returns a zero-row data.frame with the same column structure if no beats
are delineated.

## References

Goldberger, A.L., et al. (2000). "PhysioBank, PhysioToolkit, and
PhysioNet: Components of a new research resource for complex physiologic
signals." *Circulation*, 101(23), e215–e220.
[doi:10.1161/01.CIR.101.23.e215](https://doi.org/10.1161/01.CIR.101.23.e215)

## See also

[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
for R-peak detection,
[`ecgIntervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgIntervals.md)
for computing clinical ECG intervals from delineation results,
[`ecgSignalQuality`](https://x-biosignal.github.io/PhysioECG/reference/ecgSignalQuality.md)
for signal quality assessment.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
peaks <- ecgDetectRpeaks(pe)
delin <- ecgDelineate(pe, peaks)
head(delin)
} # }
```
