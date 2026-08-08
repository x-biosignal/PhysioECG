# Create Simulated ECG PhysioExperiment

Generates a PhysioExperiment object containing synthetic ECG data with
Gaussian-shaped R-peaks at regular intervals. The resulting object is
suitable for testing R-peak detection, RR interval computation, and HRV
analysis pipelines.

## Usage

``` r
make_ecg(n_time = 5000, n_channels = 1, sr = 500, heart_rate = 72)
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

## Value

A PhysioExperiment object with a single `"raw"` assay containing the
simulated ECG signal.

## References

Pan, J. & Tompkins, W.J. (1985). "A real-time QRS detection algorithm."
*IEEE Transactions on Biomedical Engineering*, 32(3), 230–236.
[doi:10.1109/TBME.1985.325532](https://doi.org/10.1109/TBME.1985.325532)

Clifford, G.D., Azuaje, F. & McSharry, P.E. (2006). *Advanced Methods
and Tools for ECG Data Analysis*. Artech House.

## See also

[`make_ecg_irregular`](https://x-biosignal.github.io/PhysioECG/reference/make_ecg_irregular.md)
for ECG with ectopic beats,
[`make_ecg_pqrst`](https://x-biosignal.github.io/PhysioECG/reference/make_ecg_pqrst.md)
for ECG with full PQRST morphology,
[`make_ecg_noisy`](https://x-biosignal.github.io/PhysioECG/reference/make_ecg_noisy.md)
for ECG with noise artifacts,
[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
for R-peak detection.

## Examples

``` r
pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)
pe
#> class: PhysioExperiment
#> dim: 5000 x 1 
#> assays(1): raw
#> samplingRate: 500 Hz
#> channels(1): ECG1
#> colData names(2): label, type
```
