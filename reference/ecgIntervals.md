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
set.seed(1)
pe <- make_ecg_pqrst(n_time = 5000, sr = 500, heart_rate = 72)$pe
peaks <- ecgDetectRpeaks(pe)
delin <- ecgDelineate(pe, peaks)
intervals <- ecgIntervals(delin, samplingRate(pe))
head(intervals)
#>   channel beat pr_ms qt_ms   qtc_ms qtc_bazett qtc_fridericia qtc_framingham
#> 1       1    1   118   414 453.3330   453.3330       439.8234        439.564
#> 2       1    2   118   476 521.2234   521.2234       505.6907        501.564
#> 3       1    3   118   410 448.9530   448.9530       435.5739        435.564
#> 4       1    4   104   498 545.3136   545.3136       529.0629        523.564
#> 5       1    5   106   504 551.8836   551.8836       535.4372        529.564
#> 6       1    6   120   566 619.7741   619.7741       601.3044        591.564
#>   qtc_hodges qrs_ms rr_ms
#> 1   434.8993     92   834
#> 2   496.8993     86   834
#> 3   430.8993     88   834
#> 4   518.8993    104   834
#> 5   524.8993    104   834
#> 6   586.8993     96   834
```
