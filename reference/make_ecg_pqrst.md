# Create Simulated ECG with PQRST Morphology

Generates a PhysioExperiment object containing synthetic ECG data with
physiologically realistic P, Q, R, S, and T wave components. Returns
both the PhysioExperiment and a data.frame of known fiducial points for
validation testing of waveform delineation algorithms.

## Usage

``` r
make_ecg_pqrst(
  n_time = 10000,
  n_channels = 1,
  sr = 500,
  heart_rate = 72,
  noise_sd = 0.02
)
```

## Arguments

- n_time:

  Number of time points (default: 10000, i.e., 20 seconds at 500 Hz).

- n_channels:

  Number of ECG channels (default: 1).

- sr:

  Sampling rate in Hz (default: 500).

- heart_rate:

  Heart rate in beats per minute (default: 72).

- noise_sd:

  Standard deviation of baseline noise (default: 0.02).

## Value

A list with two components:

- pe:

  PhysioExperiment object with the simulated ECG signal.

- fiducials:

  A data.frame with columns: `beat`, `r_peak`, `p_peak`, `q_point`,
  `s_point`, `t_peak`, `qrs_onset`, `qrs_offset` (all in sample
  indices).

## References

Goldberger, A.L., et al. (2000). "PhysioBank, PhysioToolkit, and
PhysioNet: Components of a new research resource for complex physiologic
signals." *Circulation*, 101(23), e215–e220.
[doi:10.1161/01.CIR.101.23.e215](https://doi.org/10.1161/01.CIR.101.23.e215)

Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006). *Advanced Methods
and Tools for ECG Data Analysis*. Artech House.

## See also

[`make_ecg`](https://x-biosignal.github.io/PhysioECG/reference/make_ecg.md)
for basic ECG data,
[`ecgDelineate`](https://x-biosignal.github.io/PhysioECG/reference/ecgDelineate.md)
for PQRST waveform delineation,
[`ecgIntervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgIntervals.md)
for computing clinical ECG intervals,
[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
for R-peak detection.

## Examples

``` r
result <- make_ecg_pqrst(n_time = 10000, sr = 500, heart_rate = 72)
pe <- result$pe
fiducials <- result$fiducials
head(fiducials)
#>   beat p_peak q_peak r_peak s_peak t_peak q_point s_point qrs_onset qrs_offset
#> 1    1    487    552    567    582    717     552     582       544        590
#> 2    2    904    969    984    999   1134     969     999       961       1007
#> 3    3   1321   1386   1401   1416   1551    1386    1416      1378       1424
#> 4    4   1738   1803   1818   1833   1968    1803    1833      1795       1841
#> 5    5   2155   2220   2235   2250   2385    2220    2250      2212       2258
#> 6    6   2572   2637   2652   2667   2802    2637    2667      2629       2675
```
