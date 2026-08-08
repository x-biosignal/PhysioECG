# Create Simulated ECG with Irregular R-R Intervals

Generates a PhysioExperiment object containing synthetic ECG data with
irregular beat timing. Every 5th beat is premature (60\\ interval),
followed by a compensatory pause (140\\ Useful for testing ectopic beat
detection and RR interval correction.

## Usage

``` r
make_ecg_irregular(n_time = 5000, sr = 500, heart_rate = 72)
```

## Arguments

- n_time:

  Number of time points (default: 5000).

- sr:

  Sampling rate in Hz (default: 500).

- heart_rate:

  Base heart rate in beats per minute (default: 72).

## Value

A PhysioExperiment object with a single `"raw"` assay containing the
simulated ECG signal with irregular beats.

## References

Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006). *Advanced Methods
and Tools for ECG Data Analysis*. Artech House.

Task Force of the European Society of Cardiology and the North American
Society of Pacing and Electrophysiology (1996). "Heart rate variability:
Standards of measurement, physiological interpretation and clinical
use." *Circulation*, 93(5), 1043–1065.

## See also

[`make_ecg`](https://x-biosignal.github.io/PhysioECG/reference/make_ecg.md)
for regular ECG data,
[`ecgQualityCheck`](https://x-biosignal.github.io/PhysioECG/reference/ecgQualityCheck.md)
for ectopic beat detection,
[`ecgRRcorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRcorrect.md)
for ectopic beat correction,
[`ecgRRintervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgRRintervals.md)
for RR interval computation.

## Examples

``` r
pe <- make_ecg_irregular(n_time = 5000, sr = 500, heart_rate = 72)
pe
#> class: PhysioExperiment
#> dim: 5000 x 1 
#> assays(1): raw
#> samplingRate: 500 Hz
#> channels(1): ECG1
#> colData names(2): label, type
```
