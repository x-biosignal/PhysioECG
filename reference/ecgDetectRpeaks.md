# Detect R-Peaks in ECG Signal

Identifies R-peaks in ECG data using one of four QRS detectors. All
share a bandpass front end and automatic inverted-lead handling but
differ in their feature transform and threshold rule: `"pan_tompkins"`
(5-15 Hz, derivative-square-integrate, adaptive dual threshold),
`"hamilton"` (8-16 Hz, rectified derivative, 80 ms integration, running
mean of the last eight QRS/noise peaks), `"elgendi"` (8-20 Hz, squared,
two moving averages with an event-block threshold), and `"christov"`
(9-30 Hz, rectified derivative smoothed over 30 ms, adaptive steep-slope
M-threshold with decay).

## Usage

``` r
ecgDetectRpeaks(
  x,
  method = c("pan_tompkins", "hamilton", "elgendi", "christov"),
  threshold_factor = 0.6,
  refractory_ms = 200,
  assay_name = NULL,
  beta = 0.08
)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- method:

  Detection method: one of `"pan_tompkins"` (default), `"hamilton"`,
  `"elgendi"` or `"christov"`.

- threshold_factor:

  Retained for backward compatibility; the Pan-Tompkins path uses its
  built-in fractional threshold (default: 0.6).

- refractory_ms:

  Refractory period in milliseconds. No two peaks can be closer than
  this (default: 200).

- assay_name:

  Name of the assay to use. If `NULL`, the default assay is used.

- beta:

  Offset factor for the Elgendi event threshold
  `MA_beat + beta * mean(squared signal)` (default: 0.08); ignored by
  the other methods.

## Value

A data.frame with one row per detected R-peak and the following columns:

- channel:

  Integer channel index (1-based).

- sample:

  Integer sample index of the R-peak within the assay matrix.

- time_sec:

  Time of the R-peak in seconds from signal onset.

- amplitude:

  Amplitude of the raw signal at the R-peak location (in original units,
  not inverted).

Returns a zero-row data.frame with the same column structure if no peaks
are detected.

## References

Pan, J. & Tompkins, W.J. (1985). "A real-time QRS detection algorithm."
*IEEE Transactions on Biomedical Engineering*, 32(3), 230–236.
[doi:10.1109/TBME.1985.325532](https://doi.org/10.1109/TBME.1985.325532)

Hamilton, P. (2002). "Open source ECG analysis." *Computers in
Cardiology*, 29, 101–104.

Elgendi, M. (2013). "Fast QRS detection with an optimized
knowledge-based method." *PLoS ONE*, 8(9), e73557.
[doi:10.1371/journal.pone.0073557](https://doi.org/10.1371/journal.pone.0073557)

Christov, I. (2004). "Real time electrocardiogram QRS detection using
combined adaptive threshold." *BioMedical Engineering OnLine*, 3, 28.
[doi:10.1186/1475-925X-3-28](https://doi.org/10.1186/1475-925X-3-28)

## See also

[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)
for computing RR intervals from detected peaks,
[`ecgBeatSQI`](https://x-biosignal.github.io/PhysioECG/reference/ecgBeatSQI.md)
for detector-agreement quality,
[`ecgDelineate`](https://x-biosignal.github.io/PhysioECG/reference/ecgDelineate.md)
for full waveform morphology analysis,
[`ecgSignalQuality`](https://x-biosignal.github.io/PhysioECG/reference/ecgSignalQuality.md)
for signal quality assessment.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
peaks <- ecgDetectRpeaks(pe)
head(peaks)
} # }
```
