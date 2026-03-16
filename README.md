# PhysioECG <img src="man/figures/logo.png" align="right" height="139" alt="PhysioECG logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/x-biosignal/PhysioECG/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/x-biosignal/PhysioECG/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/PhysioECG)](https://CRAN.R-project.org/package=PhysioECG)
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioECG)](https://x-biosignal.r-universe.dev/PhysioECG)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**ECG Analysis Functions for PhysioExperiment Objects**

PhysioECG provides electrocardiography (ECG) analysis functions for the PhysioExperiment ecosystem. With 18 exported functions, it delivers a complete ECG analysis pipeline -- from R-peak detection using an adaptive Pan-Tompkins algorithm, through RR interval computation and ectopic beat correction, to comprehensive heart rate variability (HRV) analysis across time-domain, frequency-domain, and nonlinear methods. It also includes ECG morphology analysis (PQRST waveform delineation, clinical interval measurement) and signal quality assessment.

## Installation

You can install PhysioECG from [r-universe](https://x-biosignal.r-universe.dev):

```r
install.packages("PhysioECG",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioECG")
```

## Quick Start

```r
library(PhysioECG)

# Generate a simulated ECG signal (10 seconds at 500 Hz, 72 bpm)
pe <- make_ecg(n_time = 5000, sr = 500, heart_rate = 72)

# Detect R-peaks using the adaptive Pan-Tompkins algorithm
peaks <- ecgDetectRpeaks(pe)
head(peaks)
#>   channel sample time_sec amplitude
#> 1       1    347    0.692      1.50
#> 2       1    764    1.526      1.50
#> ...

# Compute RR intervals
rr <- ecgRRintervals(pe, peaks)

# Time-domain HRV metrics
hrv_time <- ecgHRVtime(rr)
hrv_time
#>   channel mean_rr   sdnn rmssd pnn50 mean_hr
#> 1       1  833.33   1.23  0.87     0   72.00

# Frequency-domain HRV analysis (Welch PSD)
hrv_freq <- ecgHRVfreq(rr, method = "welch")

# Nonlinear HRV (Poincare, Sample Entropy, DFA)
hrv_nl <- ecgHRVnonlinear(rr)

# ECG morphology: delineate P, QRS, T waves
delin <- ecgDelineate(pe, peaks)
intervals <- ecgIntervals(delin, sr = 500)
head(intervals)
#>   channel beat  pr_ms  qt_ms  qtc_ms qrs_ms  rr_ms
#> 1       1    1  126.0  360.0   394.1   80.0  833.3
```

## Features

### R-Peak Detection

Adaptive dual-threshold Pan-Tompkins detector with automatic inverted-signal handling:

- `ecgDetectRpeaks()` -- detect R-peaks using bandpass filtering (5--15 Hz), differentiation, squaring, moving-window integration, and adaptive thresholding

### RR Interval Analysis

RR interval computation and ectopic beat correction:

- `ecgRRintervals()` -- compute RR intervals from detected R-peaks
- `ecgQualityCheck()` -- detect ectopic beats using local median comparison
- `ecgRRcorrect()` -- correct ectopic beats by interpolation or removal

### HRV Time-Domain Analysis

Standard time-domain HRV metrics per ESC/NASPE Task Force guidelines:

- `ecgHRVtime()` -- SDNN, RMSSD, pNN50, mean RR, mean heart rate

### HRV Frequency-Domain Analysis

Power spectral density estimation with standard frequency band integration:

- `ecgHRVfreq()` -- VLF (0.003--0.04 Hz), LF (0.04--0.15 Hz), HF (0.15--0.4 Hz) power, LF/HF ratio, total power
- Supports Welch's method (resampled FFT) and Lomb-Scargle periodogram

### HRV Nonlinear Analysis

Nonlinear dynamics and complexity measures:

- `ecgHRVnonlinear()` -- convenience wrapper combining all nonlinear metrics
- `ecgHRVpoincare()` -- Poincare plot descriptors (SD1, SD2, SD1/SD2 ratio)
- `ecgSampleEntropy()` -- sample entropy for signal complexity/regularity
- `ecgDFA()` -- detrended fluctuation analysis (alpha1 short-range, alpha2 long-range)

### ECG Morphology Analysis

Waveform delineation and clinical interval measurement:

- `ecgDelineate()` -- identify QRS complex boundaries, P-wave and T-wave peaks, and T-wave end for each beat
- `ecgIntervals()` -- compute PR interval, QRS duration, QT interval, and QTc (Bazett correction)

### Signal Quality Assessment

Signal quality evaluation and baseline correction:

- `ecgSignalQuality()` -- per-channel SNR, baseline wander, saturation ratio, and composite quality score
- `ecgBaselineCorrect()` -- remove baseline wander using highpass moving-average or running-median subtraction

### Simulated ECG Data

Synthetic ECG generators for testing and demonstration:

- `make_ecg()` -- regular ECG with Gaussian-shaped R-peaks at fixed heart rate
- `make_ecg_pqrst()` -- realistic PQRST morphology with known fiducial points for validation
- `make_ecg_irregular()` -- ECG with simulated ectopic (premature) beats and compensatory pauses
- `make_ecg_noisy()` -- ECG contaminated with baseline wander, powerline interference, and broadband noise

## Dependencies

- **R** (>= 4.2)
- **[PhysioCore](https://github.com/x-biosignal/PhysioCore)**
- **SummarizedExperiment**
- **S4Vectors**
- **stats**

## PhysioExperiment Ecosystem

PhysioECG is the ECG analysis module of the PhysioExperiment ecosystem, a suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|---------|-------------|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioFilters](https://github.com/x-biosignal/PhysioFilters) | Signal filtering and preprocessing |
| [PhysioEpoch](https://github.com/x-biosignal/PhysioEpoch) | Epoching and segmentation |
| **PhysioECG** | ECG analysis and HRV |
| [PhysioEDA](https://github.com/x-biosignal/PhysioEDA) | Electrodermal activity analysis |
| [PhysioStats](https://github.com/x-biosignal/PhysioStats) | Statistical analysis |
| [PhysioVis](https://github.com/x-biosignal/PhysioVis) | Visualization |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to browse all available packages.

## References

- Pan, J. & Tompkins, W.J. (1985). "A real-time QRS detection algorithm." *IEEE Transactions on Biomedical Engineering*, 32(3), 230--236.
- Task Force of the ESC and NASPE (1996). "Heart rate variability: Standards of measurement, physiological interpretation and clinical use." *Circulation*, 93(5), 1043--1065.
- Shaffer, F. & Ginsberg, J.P. (2017). "An overview of heart rate variability metrics and norms." *Frontiers in Public Health*, 5, 258.
- Richman, J.S. & Moorman, J.R. (2000). "Physiological time-series analysis using approximate entropy and sample entropy." *American Journal of Physiology*, 278(6), H2039--H2049.
- Peng, C.-K., et al. (1994). "Mosaic organization of DNA nucleotides." *Physical Review E*, 49(2), 1685--1689.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Author

Yusuke Matsui
