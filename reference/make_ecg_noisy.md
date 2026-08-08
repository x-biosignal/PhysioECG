# Create Simulated ECG with Noise Artifacts

Generates a PhysioExperiment object containing synthetic ECG data
contaminated with multiple noise sources: baseline wander (0.3 Hz
sinusoidal drift), powerline interference (50 Hz), and broadband
Gaussian noise. Useful for testing signal quality assessment, baseline
correction, and filtering pipelines.

## Usage

``` r
make_ecg_noisy(
  n_time = 5000,
  n_channels = 1,
  sr = 500,
  heart_rate = 72,
  baseline_amp = 0.3,
  powerline_amp = 0.1,
  noise_sd = 0.15
)
```

## Arguments

- n_time:

  Number of time points (default: 5000).

- n_channels:

  Number of ECG channels (default: 1).

- sr:

  Sampling rate in Hz (default: 500).

- heart_rate:

  Heart rate in beats per minute (default: 72).

- baseline_amp:

  Amplitude of baseline wander in arbitrary units (default: 0.3).

- powerline_amp:

  Amplitude of 50 Hz powerline noise (default: 0.1).

- noise_sd:

  Standard deviation of broadband Gaussian noise (default: 0.15).

## Value

A PhysioExperiment object with a single `"raw"` assay containing the
noisy ECG signal.

## References

Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006). *Advanced Methods
and Tools for ECG Data Analysis*. Artech House.

Shaffer, F. & Ginsberg, J.P. (2017). "An overview of heart rate
variability metrics and norms." *Frontiers in Public Health*, 5, 258.
[doi:10.3389/fpubh.2017.00258](https://doi.org/10.3389/fpubh.2017.00258)

## See also

[`make_ecg`](https://x-biosignal.github.io/PhysioECG/reference/make_ecg.md)
for clean ECG data,
[`ecgSignalQuality`](https://x-biosignal.github.io/PhysioECG/reference/ecgSignalQuality.md)
for signal quality assessment,
[`ecgBaselineCorrect`](https://x-biosignal.github.io/PhysioECG/reference/ecgBaselineCorrect.md)
for baseline wander removal,
[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
for R-peak detection.

## Examples

``` r
pe <- make_ecg_noisy(n_time = 5000, sr = 500)
pe
#> class: PhysioExperiment
#> dim: 5000 x 1 
#> assays(1): raw
#> samplingRate: 500 Hz
#> channels(1): ECG1
#> colData names(2): label, type
```
