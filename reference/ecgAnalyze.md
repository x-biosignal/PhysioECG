# Aggregate ECG HRV and Morphology Metrics

Aggregates the full HRV panel – time-domain, frequency-domain (AR
spectrum with normalized units), geometric and nonlinear – together with
a QTc morphology summary into a single per-channel data.frame, from an
object produced by
[`ecgProcess`](https://x-biosignal.github.io/PhysioECG/reference/ecgProcess.md)
(after the neurokit2 `ecg_analyze()` pattern).

## Usage

``` r
ecgAnalyze(processed)
```

## Arguments

- processed:

  A PhysioExperiment returned by
  [`ecgProcess`](https://x-biosignal.github.io/PhysioECG/reference/ecgProcess.md).

## Value

A data.frame with one row per channel and columns for the time-domain
(`sdnn`, `rmssd`, `pnn50`, `mean_hr`), frequency-domain (`lf`, `hf`,
`lf_hf_ratio`, `lf_nu`, `hf_nu`, `total_power`), geometric
(`HRV_triangular_index`, `TINN`), nonlinear (`sd1`, `sd2`,
`sample_entropy`, `alpha1`, `alpha2`) and morphology (`qt_ms`,
`qtc_bazett`, `qrs_ms`) metrics. Metrics that cannot be computed (e.g.
DFA on a short recording) are `NA`.

## See also

[`ecgProcess`](https://x-biosignal.github.io/PhysioECG/reference/ecgProcess.md),
[`ecgHRVtime`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVtime.md),
[`ecgHRVfreq`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVfreq.md),
[`ecgHRVgeometric`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVgeometric.md),
[`ecgHRVnonlinear`](https://x-biosignal.github.io/PhysioECG/reference/ecgHRVnonlinear.md),
[`ecgIntervals`](https://x-biosignal.github.io/PhysioECG/reference/ecgIntervals.md)

## Examples

``` r
pe <- make_ecg(n_time = 40000, sr = 500)
ecgAnalyze(ecgProcess(pe))
#>   channel     sdnn    rmssd pnn50  mean_hr        lf       hf lf_hf_ratio
#> 1       1 2.623303 4.369986     0 71.94245 0.1751817 1.943026  0.09015925
#>      lf_nu    hf_nu total_power HRV_triangular_index    TINN      sd1      sd2
#> 1 8.270283 91.72972     2.15185              1.62069 23.4375 3.106644 2.027856
#>   sample_entropy   alpha1     alpha2 qt_ms qtc_bazett qrs_ms
#> 1       1.642228 0.266313 0.06310743   474   525.1299    128
```
