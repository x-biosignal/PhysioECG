# Compute ECG Intervals from Delineation

Calculates standard clinical ECG intervals from the waveform delineation
produced by
[`ecgDelineate`](https://x-biosignal.github.io/PhysioECG/reference/ecgDelineate.md):
PR interval, QT interval, QTc (Bazett correction), QRS duration, and RR
interval.

## Usage

``` r
ecgIntervals(delineation, sr)
```

## Arguments

- delineation:

  A data.frame as returned by
  [`ecgDelineate`](https://x-biosignal.github.io/PhysioECG/reference/ecgDelineate.md),
  with columns `channel`, `beat`, `r_peak`, `qrs_onset`, `qrs_offset`,
  `p_peak`, `t_peak`, `t_end`.

- sr:

  Sampling rate in Hz.

## Value

A data.frame with one row per beat and the following columns:

- channel:

  Integer channel index (1-based).

- beat:

  Integer beat number within the channel.

- pr_ms:

  PR interval in milliseconds (P-wave peak to QRS onset), or `NA` if the
  P wave was not detected.

- qt_ms:

  QT interval in milliseconds (QRS onset to T-wave end), or `NA` if the
  T wave was not detected.

- qtc_ms:

  Corrected QT interval using Bazett's formula (`QT / sqrt(RR_sec)`), or
  `NA` if QT or RR is unavailable. Kept as a backward-compatible alias
  of `qtc_bazett`.

- qtc_bazett, qtc_fridericia, qtc_framingham, qtc_hodges:

  QT corrected for heart rate by, respectively, Bazett
  (`QT / sqrt(RR)`), Fridericia (`QT / RR^(1/3)`), Framingham
  (`QT + 154 * (1 - RR)`) and Hodges (`QT + 1.75 * (HR - 60)`); QT in
  ms, RR in seconds. All four equal `qt_ms` at `RR = 1 s`.

- qrs_ms:

  QRS complex duration in milliseconds.

- rr_ms:

  RR interval in milliseconds to the next beat, or `NA` for the last
  beat in each channel.

Returns a zero-row data.frame with the same column structure if no beats
are present.

## References

Goldberger, A.L., et al. (2000). "PhysioBank, PhysioToolkit, and
PhysioNet: Components of a new research resource for complex physiologic
signals." *Circulation*, 101(23), e215–e220.
[doi:10.1161/01.CIR.101.23.e215](https://doi.org/10.1161/01.CIR.101.23.e215)

## See also

[`ecgDelineate`](https://x-biosignal.github.io/PhysioECG/reference/ecgDelineate.md)
for waveform delineation,
[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
for R-peak detection,
[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)
for RR interval computation.

## Examples

``` r
if (FALSE) { # \dontrun{
pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
peaks <- ecgDetectRpeaks(pe)
delin <- ecgDelineate(pe, peaks)
intervals <- ecgIntervals(delin, samplingRate(pe))
head(intervals)
} # }
```
