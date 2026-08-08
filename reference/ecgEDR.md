# ECG-Derived Respiration (EDR)

Estimates a respiration surrogate from beat-to-beat modulation of the
QRS complex (Moody et al. 1985). Respiration shifts the cardiac
electrical axis, modulating R-wave amplitude, QRS area or upslope; the
per-beat feature series is resampled to a uniform grid to form the EDR
signal, and the respiratory rate is read from its spectral peak.

## Usage

``` r
ecgEDR(
  x,
  peaks,
  method = c("amplitude", "area", "slope"),
  fs = 4,
  resp_band = c(0.1, 0.5),
  qrs_ms = 50,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with ECG data.

- peaks:

  R-peak table from
  [`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md)
  (columns `channel`, `sample`, `time_sec`; `amplitude` is used by the
  "amplitude" method when present).

- method:

  Beat feature: "amplitude" (R-wave amplitude, default), "area" (QRS
  area) or "slope" (maximum QRS upslope).

- fs:

  Resampling rate (Hz) of the EDR signal (default 4).

- resp_band:

  Respiratory frequency band (Hz) for the rate estimate (default
  `c(0.1, 0.5)`, i.e. 6-30 breaths/min).

- qrs_ms:

  QRS half-window (ms) for the "area"/"slope" features (default 50).

- assay_name:

  Input assay name (default: first assay).

## Value

A list with:

- method:

  The feature used.

- beats:

  A data.frame (`channel`, `time_sec`, `feature`) of the per-beat
  feature.

- edr:

  A data.frame (`channel`, `time_sec`, `edr`) of the uniform-resampled
  respiration surrogate.

- resp_rate:

  A data.frame (`channel`, `resp_rate_hz`, `resp_rate_bpm`) of the
  estimated respiratory rate.

## References

Moody, G.B. et al. (1985). Derivation of respiratory signals from
multi-lead ECGs. *Computers in Cardiology*, 12, 113-116.

## See also

[`ecgDetectRpeaks`](https://x-biosignal.github.io/PhysioECG/reference/ecgDetectRpeaks.md),
[`ecgRSA`](https://x-biosignal.github.io/PhysioECG/reference/ecgRSA.md)

## Examples

``` r
pe <- make_ecg(n_time = 10000, sr = 250)
pk <- ecgDetectRpeaks(pe)
edr <- ecgEDR(pe, pk, method = "area")
edr$resp_rate
#>   channel resp_rate_hz resp_rate_bpm
#> 1       1    0.2012139      12.07284
```
